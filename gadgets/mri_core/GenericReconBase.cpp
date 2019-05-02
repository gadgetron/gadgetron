
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

        return GADGET_OK;
    }

    template <typename T>
    int GenericReconBase<T>::process(GadgetContainerMessage<T>* m1)
    {
        return GADGET_OK;
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
