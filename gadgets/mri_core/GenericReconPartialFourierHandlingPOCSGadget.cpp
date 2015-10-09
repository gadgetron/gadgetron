
#include "GenericReconPartialFourierHandlingPOCSGadget.h"
#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron {

    GenericReconPartialFourierHandlingPOCSGadget::GenericReconPartialFourierHandlingPOCSGadget() : BaseClass()
    {
    }

    GenericReconPartialFourierHandlingPOCSGadget::~GenericReconPartialFourierHandlingPOCSGadget()
    {
    }

    int GenericReconPartialFourierHandlingPOCSGadget::perform_partial_fourier_handling()
    {
        if (perform_timing.value()) { gt_timer_.start("GenericReconPartialFourierHandlingPOCSGadget, partial_fourier_filter"); }

        GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::partial_fourier_POCS(kspace_buf_,
            startRO_, endRO_, startE1_, endE1_, startE2_, endE2_,
            partial_fourier_POCS_transitBand.value(), partial_fourier_POCS_transitBand.value(),
            partial_fourier_POCS_transitBand_E2.value(), partial_fourier_POCS_iters.value(),
            partial_fourier_POCS_thres.value(), pf_res_), GADGET_FAIL);

        if (perform_timing.value()) { gt_timer_.stop(); }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(GenericReconPartialFourierHandlingPOCSGadget)

}
