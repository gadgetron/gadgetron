
#include <GadgetronTimer.h>
#include "GenericReconPartialFourierHandlingPOCSGadget.h"
#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron {


hoNDArray<std::complex<float>> GenericReconPartialFourierHandlingPOCSGadget::perform_partial_fourier_handling(const hoNDArray<std::complex<float>> & kspace_buffer, size_t start_RO, size_t end_RO, size_t start_E1, size_t end_E1, size_t start_E2, size_t end_E2) const{

        std::optional<GadgetronTimer> gt_timer;
        if (perform_timing) { gt_timer = GadgetronTimer("GenericReconPartialFourierHandlingFilterGadget, partial_fourier_filter"); }

        hoNDArray<std::complex<float>> pf_res;

        Gadgetron::partial_fourier_POCS(kspace_buffer,
            start_RO,end_RO,start_E1,end_E1,start_E2,end_E2,
            partial_fourier_POCS_transitBand, partial_fourier_POCS_transitBand,
            partial_fourier_POCS_transitBand_E2, partial_fourier_POCS_iters,
            partial_fourier_POCS_thres, pf_res);

        return std::move(pf_res);

}

GADGETRON_GADGET_EXPORT(GenericReconPartialFourierHandlingPOCSGadget)
}
