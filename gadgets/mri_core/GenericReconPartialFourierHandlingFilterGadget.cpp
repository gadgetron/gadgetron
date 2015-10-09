
#include "GenericReconPartialFourierHandlingFilterGadget.h"
#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron { 

    GenericReconPartialFourierHandlingFilterGadget::GenericReconPartialFourierHandlingFilterGadget() : BaseClass()
    {
    }

    GenericReconPartialFourierHandlingFilterGadget::~GenericReconPartialFourierHandlingFilterGadget()
    {
    }

    int GenericReconPartialFourierHandlingFilterGadget::perform_partial_fourier_handling()
    {
        if (perform_timing.value()) { gt_timer_.start("GenericReconPartialFourierHandlingFilterGadget, partial_fourier_filter"); }

        GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::partial_fourier_filter(kspace_buf_,
            startRO_, endRO_, startE1_, endE1_, startE2_, endE2_,
            partial_fourier_filter_RO_width.value(), partial_fourier_filter_E1_width.value(),
            partial_fourier_filter_E2_width.value(), partial_fourier_filter_densityComp.value(),
            filter_pf_RO_, filter_pf_E1_, filter_pf_E2_, pf_res_), GADGET_FAIL);

        if (perform_timing.value()) { gt_timer_.stop(); }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(GenericReconPartialFourierHandlingFilterGadget)

}
