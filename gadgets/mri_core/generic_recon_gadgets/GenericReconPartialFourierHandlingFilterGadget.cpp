
#include <GadgetronTimer.h>
#include "GenericReconPartialFourierHandlingFilterGadget.h"
#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron { 



    hoNDArray <std::complex<float>> GenericReconPartialFourierHandlingFilterGadget::perform_partial_fourier_handling(const hoNDArray <std::complex<float>> &kspace_buf, size_t start_RO, size_t end_RO,
                                     size_t start_E1, size_t end_E1, size_t start_E2, size_t end_E2) const
    {

        Core::optional<GadgetronTimer> gt_timer;
        if (perform_timing) { gt_timer = GadgetronTimer("GenericReconPartialFourierHandlingFilterGadget, partial_fourier_filter"); }


        std::lock_guard<std::mutex> guard(filter_mutex);

        hoNDArray<std::complex<float>> pf_res;

        Gadgetron::partial_fourier_filter(kspace_buf,
            start_RO, end_RO, start_E1, end_E1, start_E2, end_E2,
            partial_fourier_filter_RO_width, partial_fourier_filter_E1_width,
            partial_fourier_filter_E2_width, partial_fourier_filter_densityComp,
            filter_pf_RO_, filter_pf_E1_, filter_pf_E2_, pf_res);

        return pf_res;
    }

    // ----------------------------------------------------------------------------------------

    GADGETRON_GADGET_EXPORT(GenericReconPartialFourierHandlingFilterGadget)

}
