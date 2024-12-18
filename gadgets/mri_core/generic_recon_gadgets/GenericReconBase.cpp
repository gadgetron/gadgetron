
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
    int GenericReconBase<T>::process_config(const mrd::Header& header)
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
        this->gt_streamer_.initialize_stream_name_buffer(this->context.parameters);
        this->gt_streamer_.verbose_ = this->verbose.value();

        return GADGET_OK;
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

    template class GenericReconBase<mrd::ReconData>;
    template class GenericReconBase<mrd::ImageArray>;
    template class GenericReconBase<mrd::ImageHeader>;

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
