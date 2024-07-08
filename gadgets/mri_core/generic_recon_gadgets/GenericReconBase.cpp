
#include "GenericReconBase.h"
#include <boost/filesystem.hpp>

namespace Gadgetron {

    template <typename T> 
    GenericReconBase<T>::GenericReconBase() : num_encoding_spaces_(1), process_called_times_(0)
    {
        gt_timer_.set_timing_in_destruction(false);
        gt_timer_local_.set_timing_in_destruction(false);
    }

    template <typename T> 
    GenericReconBase<T>::~GenericReconBase()
    {
    }

    template <typename T> 
    int GenericReconBase<T>::process_config(ACE_Message_Block* mb)
    {
        if (!debug_folder.value().empty())
        {
            Gadgetron::get_debug_folder_path(debug_folder.value(), debug_folder_full_path_);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is " << debug_folder_full_path_);

            // Create debug folder if necessary
            boost::filesystem::path boost_folder_path(debug_folder_full_path_);
            try
            {
                boost::filesystem::create_directories(boost_folder_path);
            }
            catch (...)
            {
                GERROR("Error creating the debug folder.\n");
                return false;
            }
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is not set ... ");
        }

        // find the buffer names if they are set
        this->initialize_stream_name_buffer(GENERIC_RECON_ISMRMRD_HEADER);
        this->initialize_stream_name_buffer(GENERIC_RECON_UNDERSAMPLED_KSPACE);
        this->initialize_stream_name_buffer(GENERIC_RECON_REF_KSPACE);
        this->initialize_stream_name_buffer(GENERIC_RECON_REF_KSPACE_FOR_COILMAP);
        this->initialize_stream_name_buffer(GENERIC_RECON_COILMAP);
        this->initialize_stream_name_buffer(GENERIC_RECON_GFACTOR_MAP);
        this->initialize_stream_name_buffer(GENERIC_RECON_RECONED_KSPACE);
        this->initialize_stream_name_buffer(GENERIC_RECON_RECONED_COMPLEX_IMAGE);
        this->initialize_stream_name_buffer(GENERIC_RECON_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING);

        return GADGET_OK;
    }

    template <typename T> 
    void GenericReconBase<T>::initialize_stream_name_buffer(const std::string& name)
    {
        if (this->context.parameters.find(name) != this->context.parameters.end())
        {
            buffer_names_[name].first = this->context.parameters.at(name);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store " << name << " is " << buffer_names_[name].first);
        }
    }

    template <typename T> 
    void GenericReconBase<T>::close_stream_buffer()
    {
        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GenericReconBase<T>::close_stream_buffer ");

        // put the close message to all opened streams for the images into the ismrmrd format
        // the hoNDArray stream is not in ismrmrd format for now, which may be changed in the future
        for (auto const& x : this->buffer_names_)
        {
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "GenericReconBase<T>::close_stream_buffer, stream is for " << x.first << " - " << this->buffer_names_[x.first].first);

            if ( (x.first == GENERIC_RECON_ISMRMRD_HEADER)  
                | (x.first == GENERIC_RECON_COILMAP) 
                | (x.first == GENERIC_RECON_GFACTOR_MAP)
                | (x.first == GENERIC_RECON_RECONED_COMPLEX_IMAGE)
                | (x.first == GENERIC_RECON_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING)
            )
            {
                if(this->buffer_names_[x.first].second)
                {
                    std::ofstream& os = *this->buffer_names_[x.first].second;
                    if (os.is_open())
                    {
                        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GenericReconBase<T>::close_stream_buffer, stream is open for " << x.first << "; put in the close message ... ");
                        ISMRMRD::OStreamView ws(os);
                        ISMRMRD::ProtocolSerializer serializer(ws);
                        serializer.close();
                        os.flush();
                    }
                    else
                    {
                        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GenericReconBase<T>::close_stream_buffer, stream is not open for " << x.first << " ... ");
                    }
                }
            }
        }
    }

    template <typename T>
    int GenericReconBase<T>::process(GadgetContainerMessage<T>* m1)
    {
        return GADGET_OK;
    }

    template <typename T> 
    int GenericReconBase<T>::close(unsigned long flags)
    {
        GDEBUG_CONDITION_STREAM(this->verbose.value(), "GenericReconBase<T> - close(flags) : " << flags);
        if (BaseClass::close(flags) != GADGET_OK) return GADGET_FAIL;
        return GADGET_OK;
    }

    template <typename T> 
    void GenericReconBase<T>::stream_ismrmrd_header(const ISMRMRD::IsmrmrdHeader& hdr)
    {
        if (this->buffer_names_.find(GENERIC_RECON_ISMRMRD_HEADER)!=this->buffer_names_.end())
        {
            std::string buf_name = this->buffer_names_[GENERIC_RECON_ISMRMRD_HEADER].first;
            if (!this->buffer_names_[GENERIC_RECON_ISMRMRD_HEADER].second)
            {
                this->buffer_names_[GENERIC_RECON_ISMRMRD_HEADER].second = std::make_shared<std::ofstream>(std::ofstream(buf_name, std::ios::out ));
            }

            std::ofstream& os = *this->buffer_names_[GENERIC_RECON_ISMRMRD_HEADER].second;
            if (os.is_open())
            {
                GDEBUG_STREAM("Generic recon, stream the ismrmrd header to the array buffer " << buf_name);
                ISMRMRD::serialize(hdr, os);
                os.flush();
            }
            else
            {
                GERROR_STREAM("Generic recon, unable to open the ismrmrd header buffer " << buf_name << " ... ");
            }
        }
        else
        {
            GWARN_CONDITION_STREAM(this->verbose.value(), "Generic reconstruction chain, the pre-set buffer names do not include " << GENERIC_RECON_ISMRMRD_HEADER << "; the header will not be saved into the buffer ...");
        }
    }

    template class EXPORTGADGETSMRICORE GenericReconBase<IsmrmrdReconData>;
    template class EXPORTGADGETSMRICORE GenericReconBase<IsmrmrdImageArray>;
    template class EXPORTGADGETSMRICORE GenericReconBase<ISMRMRD::ImageHeader>;

    GenericReconKSpaceReadoutBase::GenericReconKSpaceReadoutBase() : BaseClass()
    {
    }

    GenericReconKSpaceReadoutBase::~GenericReconKSpaceReadoutBase()
    {
    }

    GenericReconDataBase::GenericReconDataBase() : BaseClass()
    {
    }

    GenericReconDataBase::~GenericReconDataBase()
    {
    }

    GenericReconImageBase::GenericReconImageBase() : BaseClass()
    {
    }

    GenericReconImageBase::~GenericReconImageBase()
    {
    }

    GenericReconImageHeaderBase::GenericReconImageHeaderBase() : BaseClass()
    {
    }

    GenericReconImageHeaderBase::~GenericReconImageHeaderBase()
    {
    }

    GADGET_FACTORY_DECLARE(GenericReconKSpaceReadoutBase)
    GADGET_FACTORY_DECLARE(GenericReconDataBase)
    GADGET_FACTORY_DECLARE(GenericReconImageBase)
    GADGET_FACTORY_DECLARE(GenericReconImageHeaderBase)
}
