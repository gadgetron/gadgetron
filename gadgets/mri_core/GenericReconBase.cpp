
#include "GenericReconBase.h"

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

    GADGET_FACTORY_DECLARE(GenericReconDataBase)
    GADGET_FACTORY_DECLARE(GenericReconImageBase)
}
