
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
        if (this->context.parameters.find(GENERIC_RECON_ISMRMRD_HEADER) != this->context.parameters.end())
        {
            buffer_names_[GENERIC_RECON_ISMRMRD_HEADER] = this->context.parameters.at(GENERIC_RECON_ISMRMRD_HEADER);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store the ismrmrd header is " << buffer_names_[GENERIC_RECON_ISMRMRD_HEADER]);
        }

        if (this->context.parameters.find(GENERIC_RECON_UNDERSAMPLED_KSPACE) != this->context.parameters.end())
        {
            buffer_names_[GENERIC_RECON_UNDERSAMPLED_KSPACE] = this->context.parameters.at(GENERIC_RECON_UNDERSAMPLED_KSPACE);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store the undersampled kspace is " << buffer_names_[GENERIC_RECON_UNDERSAMPLED_KSPACE]);
        }

        if (this->context.parameters.find(GENERIC_RECON_REF_KSPACE) != this->context.parameters.end())
        {
            buffer_names_[GENERIC_RECON_REF_KSPACE] = this->context.parameters.at(GENERIC_RECON_REF_KSPACE);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store the reference kspace is " << buffer_names_[GENERIC_RECON_REF_KSPACE]);
        }

        if (this->context.parameters.find(GENERIC_RECON_REF_KSPACE_FOR_COILMAP) != this->context.parameters.end())
        {
            buffer_names_[GENERIC_RECON_REF_KSPACE_FOR_COILMAP] = this->context.parameters.at(GENERIC_RECON_REF_KSPACE_FOR_COILMAP);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store the prepared reference kspace used for the coil map estimation is " << buffer_names_[GENERIC_RECON_REF_KSPACE_FOR_COILMAP]);
        }

        if (this->context.parameters.find(GENERIC_RECON_COILMAP) != this->context.parameters.end())
        {
            buffer_names_[GENERIC_RECON_COILMAP] = this->context.parameters.at(GENERIC_RECON_COILMAP);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store the coil maps is " << buffer_names_[GENERIC_RECON_COILMAP]);
        }

        if (this->context.parameters.find(GENERIC_RECON_GFACTOR_MAP) != this->context.parameters.end())
        {
            buffer_names_[GENERIC_RECON_GFACTOR_MAP] = this->context.parameters.at(GENERIC_RECON_GFACTOR_MAP);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store the gfactor maps is " << buffer_names_[GENERIC_RECON_GFACTOR_MAP]);
        }

        if (this->context.parameters.find(GENERIC_RECON_RECONED_KSPACE) != this->context.parameters.end())
        {
            buffer_names_[GENERIC_RECON_RECONED_KSPACE] = this->context.parameters.at(GENERIC_RECON_RECONED_KSPACE);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store the reconstructed kspace is " << buffer_names_[GENERIC_RECON_RECONED_KSPACE]);
        }

        if (this->context.parameters.find(GENERIC_RECON_RECONED_COMPLEX_IMAGE) != this->context.parameters.end())
        {
            buffer_names_[GENERIC_RECON_RECONED_COMPLEX_IMAGE] = this->context.parameters.at(GENERIC_RECON_RECONED_COMPLEX_IMAGE);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store the complex images after recon is " << buffer_names_[GENERIC_RECON_RECONED_COMPLEX_IMAGE]);
        }

        if (this->context.parameters.find(GENERIC_RECON_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING) != this->context.parameters.end())
        {
            buffer_names_[GENERIC_RECON_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING] = this->context.parameters.at(GENERIC_RECON_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING);
            GDEBUG_CONDITION_STREAM(this->verbose.value(), "Buffer to store the complex images after recon and further post-processing is " << buffer_names_[GENERIC_RECON_RECONED_COMPLEX_IMAGE_AFTER_POSTPROCESSING]);
        }

        return GADGET_OK;
    }

    template <typename T>
    int GenericReconBase<T>::process(GadgetContainerMessage<T>* m1)
    {
        return GADGET_OK;
    }

    template <typename T> 
    void GenericReconBase<T>::stream_ismrmrd_header(const ISMRMRD::IsmrmrdHeader& hdr)
    {
        if (this->buffer_names_.find(GENERIC_RECON_ISMRMRD_HEADER)!=this->buffer_names_.end())
        {
            std::string buf_name = this->buffer_names_[GENERIC_RECON_ISMRMRD_HEADER];
            std::ofstream os(buf_name, std::ios::out );
            if (os.is_open())
            {
                GDEBUG_STREAM("Generic recon, stream the ismrmrd header to the array buffer " << buf_name);
                ISMRMRD::serialize(hdr, os);
                os.flush();
                os.close();
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
