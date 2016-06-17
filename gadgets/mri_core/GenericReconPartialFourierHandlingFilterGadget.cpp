
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

        // this PF sampling related SNR scaling should be done in the recon step
        // the PF filter is only a noise level preserving filter
        /*long lenRO = endRO_ - startRO_ + 1;
        long lenE1 = endE1_ - startE1_ + 1;
        long lenE2 = endE2_ - startE2_ + 1;

        size_t RO = kspace_buf_.get_size(0);
        size_t E1 = kspace_buf_.get_size(1);
        size_t E2 = kspace_buf_.get_size(2);

        real_value_type partialFourierCompensationFactor = 1;

        if (lenRO < RO)
        {
            partialFourierCompensationFactor *= (real_value_type)(RO) / (real_value_type)(lenRO);
        }

        if (lenE1 < E1)
        {
            partialFourierCompensationFactor *= (real_value_type)(E1) / (real_value_type)(lenE1);
        }

        if (E2>1 && lenE2 < E2)
        {
            partialFourierCompensationFactor *= (real_value_type)(E2) / (real_value_type)(lenE2);
        }

        partialFourierCompensationFactor = std::sqrt(partialFourierCompensationFactor);
        if (verbose.value())
        {
            GDEBUG_STREAM("GenericReconPartialFourierHandlingFilterGadget, partial fourier scaling factor : " << partialFourierCompensationFactor);
        }

        if (partialFourierCompensationFactor>1)
        {
            Gadgetron::scal(partialFourierCompensationFactor, kspace_buf_);
        }*/

        GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::partial_fourier_filter(kspace_buf_,
            startRO_, endRO_, startE1_, endE1_, startE2_, endE2_,
            partial_fourier_filter_RO_width.value(), partial_fourier_filter_E1_width.value(),
            partial_fourier_filter_E2_width.value(), partial_fourier_filter_densityComp.value(),
            filter_pf_RO_, filter_pf_E1_, filter_pf_E2_, pf_res_), GADGET_FAIL);

        if (perform_timing.value()) { gt_timer_.stop(); }

        /*if (!debug_folder_full_path_.empty())
        {
            if (filter_pf_RO_.get_number_of_elements()>0) gt_exporter_.export_array_complex(filter_pf_RO_, debug_folder_full_path_ + "filter_pf_RO");
            if (filter_pf_E1_.get_number_of_elements()>0) gt_exporter_.export_array_complex(filter_pf_E1_, debug_folder_full_path_ + "filter_pf_E1");
            if (filter_pf_E2_.get_number_of_elements()>0) gt_exporter_.export_array_complex(filter_pf_E2_, debug_folder_full_path_ + "filter_pf_E2");
        }*/

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    GADGET_FACTORY_DECLARE(GenericReconPartialFourierHandlingFilterGadget)

}
