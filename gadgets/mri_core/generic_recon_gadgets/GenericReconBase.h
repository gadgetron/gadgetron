/** \file   GenericReconBase.h
    \brief  This serves an optional base class gadget for the generic chain.
            Some common functionalities are implemented here and can be reused in specific recon gadgets.
            This gadget is instantiated for IsmrmrdReconData and IsmrmrdImageArray
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "gadgetron_mricore_export.h"
#include "Gadget.h"
#include "GadgetronTimer.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"
#include "ismrmrd/serialization.h"
#include "ismrmrd/serialization_iostream.h"

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"

#include "ImageIOAnalyze.h"

#include "gadgetron_sha1.h"

#include "GenericReconStreamDef.h"

namespace Gadgetron {

    template <typename T> 
    class EXPORTGADGETSMRICORE GenericReconBase : public Gadget1<T>
    {
    public:
        GADGET_DECLARE(GenericReconBase);

        typedef Gadget1<T> BaseClass;

        GenericReconBase();
        ~GenericReconBase();

        /// ------------------------------------------------------------------------------------
        /// debug and timing
        GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);
        GADGET_PROPERTY(debug_folder, std::string, "If set, the debug output will be written out", "");
        GADGET_PROPERTY(perform_timing, bool, "Whether to perform timing on some computational steps", false);

        /// ms for every time tick
        GADGET_PROPERTY(time_tick, float, "Time tick in ms", 2.5);

    protected:

        // number of encoding spaces in the protocol
        size_t num_encoding_spaces_;

        // number of times the process function is called
        size_t process_called_times_;

        // --------------------------------------------------
        // variables for debug and timing
        // --------------------------------------------------

        // clock for timing
        Gadgetron::GadgetronTimer gt_timer_local_;
        Gadgetron::GadgetronTimer gt_timer_;

        // debug folder
        std::string debug_folder_full_path_;

        // exporter
        Gadgetron::ImageIOAnalyze gt_exporter_;

        // store buffer names
        std::map<std::string, std::pair<std::string, std::shared_ptr<std::ofstream> > > buffer_names_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(GadgetContainerMessage<T>* m1);
        virtual int close(unsigned long flags);

        // --------------------------------------------------
        // data stream functions
        // --------------------------------------------------

        void initialize_stream_name_buffer(const std::string& name);
        void close_stream_buffer();

        void stream_ismrmrd_header(const ISMRMRD::IsmrmrdHeader& hdr);

        // stream of ND array buffer
        template <typename DataType> 
        void stream_to_array_buffer(const std::string& name, const hoNDArray<DataType>& data)
        {
            if (this->buffer_names_.find(name)!=this->buffer_names_.end())
            {
                std::string buf_name = this->buffer_names_[name].first;

                if (!this->buffer_names_[name].second)
                {
                    GDEBUG_STREAM("Generic recon, create the stream for the first time - " << buf_name);
                    this->buffer_names_[name].second = std::make_shared<std::ofstream>(std::ofstream(buf_name, std::ios::out | std::ios::binary | std::ios::app));
                }

                std::ofstream& os = *this->buffer_names_[name].second;
                if (os.is_open())
                {
                    GDEBUG_STREAM("Generic recon, continue streaming the data to the buffer " << buf_name);
                    Gadgetron::Core::IO::write(os, data);
                    os.flush();
                }
                else
                {
                    GERROR_STREAM("Generic recon, already created stream is not in open status - " << buf_name << " ... ");
                }
            }
            else
            {
                GWARN_CONDITION_STREAM(this->verbose.value(), "Generic reconstruction chain, the pre-set buffer names do not include " << name << "; the data will not be saved into the buffer ...");
            }
        }

        // stream to ismrmrd image stream, e.g. for coil maps, g-maps and reconed images

        template <typename DataType> 
        void stream_images(std::ofstream& os, const std::vector< ISMRMRD::Image<DataType> >& ismrmrd_images)
        {
            ISMRMRD::OStreamView ws(os);
            ISMRMRD::ProtocolSerializer serializer(ws);

            for (auto im : ismrmrd_images)
            {
                serializer.serialize(im);
            }
            os.flush();
        }

        template <typename DataType> 
        void stream_to_ismrmrd_image_buffer(const std::string& name, const hoNDArray<DataType>& img, const hoNDArray< ISMRMRD::ImageHeader >& headers, const std::vector< ISMRMRD::MetaContainer >& meta)
        {
            if (this->buffer_names_.find(name)!=this->buffer_names_.end())
            {
                std::string buf_name = this->buffer_names_[name].first;

                // convert images to one or more ismrmrd images
                std::vector< ISMRMRD::Image<DataType> > ismrmrd_images;

                size_t RO = img.get_size(0);
                size_t E1 = img.get_size(1);
                size_t E2 = img.get_size(2);
                size_t CHA = img.get_size(3);

                size_t N = img.get_size(4);
                size_t S = img.get_size(5);
                size_t SLC = img.get_size(6);

                GDEBUG_STREAM("Generic recon, convert recon images to ismrmd images for " << name << " [RO E1 E2 CHA N S SLC] = [" << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC << "] - " << buf_name);

                ismrmrd_images.resize(N*S*SLC);

                for (auto slc=0; slc<SLC; slc++)
                {
                    for (auto s=0; s<S; s++)
                    {
                        for (auto n=0; n<N; n++)
                        {
                            size_t ind = n+s*N+slc*N*S;

                            ISMRMRD::Image<DataType>& a_img = ismrmrd_images[ind];

                            a_img.resize(RO, E1, E2, CHA);
                            memcpy(a_img.getDataPtr(), &img(0, 0, 0, 0, n, s, slc), sizeof(DataType)*RO*E1*E2*CHA);

                            ISMRMRD::ImageHeader hd = headers(n, s, slc);
                            hd.data_type = Gadgetron::Core::IO::ismrmrd_data_type<DataType>();
                            a_img.setHead(hd);

                            std::ostringstream str;
                            ISMRMRD::serialize(meta[ind], str);
                            a_img.setAttributeString(str.str());
                        }
                    }
                }

                if (!this->buffer_names_[name].second)
                {
                    GDEBUG_STREAM("Generic recon, create the ismrmrd image stream for the first time - " << buf_name);
                    this->buffer_names_[name].second = std::make_shared<std::ofstream>(std::ofstream(buf_name, std::ios::out | std::ios::binary | std::ios::app));
                }

                if (this->buffer_names_[name].second->is_open())
                {
                    GDEBUG_STREAM("Generic recon, continue streaming the data to the ismrmrd image buffer " << buf_name);
                    stream_images(*this->buffer_names_[name].second, ismrmrd_images);
                }
                else
                {
                    GERROR_STREAM("Generic recon, already created ismrmrd image stream is not in open status - " << buf_name << " ... ");
                }
            }
            else
            {
                GWARN_CONDITION_STREAM(this->verbose.value(), "Generic reconstruction chain, the pre-set buffer names do not include " << name << "; the data will not be saved into the buffer ...");
            }
        }
    };

    class EXPORTGADGETSMRICORE GenericReconKSpaceReadoutBase :public GenericReconBase < ISMRMRD::AcquisitionHeader >
    {
    public:
        GADGET_DECLARE(GenericReconKSpaceReadoutBase);

        typedef GenericReconBase < ISMRMRD::AcquisitionHeader > BaseClass;

        GenericReconKSpaceReadoutBase();
        virtual ~GenericReconKSpaceReadoutBase();
        virtual int close(unsigned long flags) { return BaseClass::close(flags); }
    };

    class EXPORTGADGETSMRICORE GenericReconDataBase :public GenericReconBase < IsmrmrdReconData >
    {
    public:
        GADGET_DECLARE(GenericReconDataBase);

        typedef GenericReconBase < IsmrmrdReconData > BaseClass;

        GenericReconDataBase();
        virtual ~GenericReconDataBase();
        virtual int close(unsigned long flags) { return BaseClass::close(flags); }
    };

    class EXPORTGADGETSMRICORE GenericReconImageBase :public GenericReconBase < IsmrmrdImageArray >
    {
    public:
        GADGET_DECLARE(GenericReconImageBase);

        typedef GenericReconBase < IsmrmrdImageArray > BaseClass;

        GenericReconImageBase();
        virtual ~GenericReconImageBase();
        virtual int close(unsigned long flags) { return BaseClass::close(flags); }
    };

    class EXPORTGADGETSMRICORE GenericReconImageHeaderBase :public GenericReconBase < ISMRMRD::ImageHeader >
    {
    public:
        GADGET_DECLARE(GenericReconImageHeaderBase);

        typedef GenericReconBase < ISMRMRD::ImageHeader > BaseClass;

        GenericReconImageHeaderBase();
        virtual ~GenericReconImageHeaderBase();
        virtual int close(unsigned long flags) { return BaseClass::close(flags); }
    };
}
