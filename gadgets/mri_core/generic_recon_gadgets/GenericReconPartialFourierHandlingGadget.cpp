
#include "GenericReconPartialFourierHandlingGadget.h"
#include <iomanip>
#include <GadgetronTimer.h>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "mri_core_kspace_filter.h"
#include "mri_core_utility.h"

namespace Gadgetron {
    GenericReconPartialFourierHandlingGadget::GenericReconPartialFourierHandlingGadget(
        const Core::Context& context, const Core::GadgetProperties& props) : BaseClass(context,props)
    {
        const auto& h = context.header;

        size_t NE = h.encoding.size();

        num_encoding_spaces = NE;

        GDEBUG_CONDITION_STREAM(verbose, "Number of encoding spaces: " << NE);

        acceFactorE1_.resize(NE, 1);
        acceFactorE2_.resize(NE, 1);

        size_t e;
        for (e = 0; e < h.encoding.size(); e++) {
            if (!h.encoding[e].parallel_imaging) {
                GDEBUG_STREAM("Parallel Imaging section not found in header for encoding " << e);
                acceFactorE1_[e] = 1;
                acceFactorE2_[e] = 1;
            } else {
                auto p_imaging = *h.encoding[0].parallel_imaging;

                acceFactorE1_[e] = p_imaging.acceleration_factor.kspace_encoding_step_1;
                acceFactorE2_[e] = p_imaging.acceleration_factor.kspace_encoding_step_2;
                GDEBUG_CONDITION_STREAM(verbose, "acceFactorE1 is " << acceFactorE1_[e]);
                GDEBUG_CONDITION_STREAM(verbose, "acceFactorE2 is " << acceFactorE2_[e]);
            }
        }
    }

    mrd::ImageArray GenericReconPartialFourierHandlingGadget::process_function(mrd::ImageArray recon_res) const {
        std::optional<GadgetronTimer> gt_timer;
        if (perform_timing) {
            gt_timer = GadgetronTimer("GenericReconPartialFourierHandlingGadget::process");
        }

        GDEBUG_CONDITION_STREAM(verbose, "GenericReconPartialFourierHandlingGadget::process(...) starts ... ");

        // some images do not need partial fourier handling processing
        if (recon_res.meta[0].count(skip_processing_meta_field) && recon_res.meta[0][skip_processing_meta_field].size() > 0) {
            return std::move(recon_res);
        }

        // call the partial foureir

        size_t encoding = (size_t)std::get<long>(recon_res.meta[0]["encoding"].front());
        if (encoding > num_encoding_spaces) throw std::runtime_error("Illegal number of encoding spaces provided");

        // std::string dataRole = std::string(recon_res.meta[0].as_str(GADGETRON_DATA_ROLE));

        // std::stringstream os;
        // os << "encoding_" << encoding << "_" << dataRole;
        // std::string str = os.str();

        size_t RO  = recon_res.data.get_size(0);
        size_t E1  = recon_res.data.get_size(1);
        size_t E2  = recon_res.data.get_size(2);
        size_t CHA = recon_res.data.get_size(3);
        size_t N   = recon_res.data.get_size(4);
        size_t S   = recon_res.data.get_size(5);
        size_t SLC = recon_res.data.get_size(6);

        // perform SNR unit scaling
        mrd::SamplingLimits sampling_limits;

        if (recon_res.meta[0].count("sampling_limits_RO")) {
            auto& sl = recon_res.meta[0]["sampling_limits_RO"];
            sampling_limits.kspace_encoding_step_0.minimum = (uint32_t)std::get<long>(sl[0]);
            sampling_limits.kspace_encoding_step_0.center = (uint32_t)std::get<long>(sl[1]);
            sampling_limits.kspace_encoding_step_0.maximum = (uint32_t)std::get<long>(sl[2]);
        }

        if (!((sampling_limits.kspace_encoding_step_0.minimum >= 0) && (sampling_limits.kspace_encoding_step_0.maximum < RO)
                && (sampling_limits.kspace_encoding_step_0.minimum <= sampling_limits.kspace_encoding_step_0.maximum))) {
            sampling_limits.kspace_encoding_step_0.minimum    = 0;
            sampling_limits.kspace_encoding_step_0.center = RO / 2;
            sampling_limits.kspace_encoding_step_0.maximum    = RO - 1;
        }

        if (recon_res.meta[0].count("sampling_limits_E1")) {
            auto& sl = recon_res.meta[0]["sampling_limits_E1"];
            sampling_limits.kspace_encoding_step_1.minimum = (uint32_t)std::get<long>(sl[0]);
            sampling_limits.kspace_encoding_step_1.center = (uint32_t)std::get<long>(sl[1]);
            sampling_limits.kspace_encoding_step_1.maximum = (uint32_t)std::get<long>(sl[2]);
        }

        if (!((sampling_limits.kspace_encoding_step_1.minimum >= 0) && (sampling_limits.kspace_encoding_step_1.maximum < E1)
                && (sampling_limits.kspace_encoding_step_1.minimum <= sampling_limits.kspace_encoding_step_1.maximum))) {
            sampling_limits.kspace_encoding_step_1.minimum    = 0;
            sampling_limits.kspace_encoding_step_1.center = E1 / 2;
            sampling_limits.kspace_encoding_step_1.maximum    = E1 - 1;
        }

        if (recon_res.meta[0].count("sampling_limits_E2")) {
            auto& sl = recon_res.meta[0]["sampling_limits_E2"];
            sampling_limits.kspace_encoding_step_2.minimum = (uint32_t)std::get<long>(sl[0]);
            sampling_limits.kspace_encoding_step_2.center = (uint32_t)std::get<long>(sl[1]);
            sampling_limits.kspace_encoding_step_2.maximum = (uint32_t)std::get<long>(sl[2]);
        }

        if (!((sampling_limits.kspace_encoding_step_2.minimum >= 0) && (sampling_limits.kspace_encoding_step_2.maximum < E2)
                && (sampling_limits.kspace_encoding_step_2.minimum <= sampling_limits.kspace_encoding_step_2.maximum))) {
            sampling_limits.kspace_encoding_step_2.minimum    = 0;
            sampling_limits.kspace_encoding_step_2.center = E2 / 2;
            sampling_limits.kspace_encoding_step_2.maximum    = E2 - 1;
        }

        // ----------------------------------------------------------
        // pf kspace sampling range
        // ----------------------------------------------------------
        // if image padding is performed, those dimension may not need partial fourier handling

        size_t startRO_ = sampling_limits.kspace_encoding_step_0.minimum;
        size_t endRO_   = sampling_limits.kspace_encoding_step_0.maximum;

        size_t startE1_ = 0;
        size_t endE1_   = E1 - 1;

        size_t startE2_ = 0;
        size_t endE2_   = E2 - 1;

        if (std::abs((double)(sampling_limits.kspace_encoding_step_1.maximum - E1 / 2) - (double)(E1 / 2 - sampling_limits.kspace_encoding_step_1.minimum))
            > acceFactorE1_[encoding]) {
            startE1_ = sampling_limits.kspace_encoding_step_1.minimum;
            endE1_   = sampling_limits.kspace_encoding_step_1.maximum;
        }

        if ((E2 > 1)
            && (std::abs((double)(sampling_limits.kspace_encoding_step_2.maximum - E2 / 2) - (double)(E2 / 2 - sampling_limits.kspace_encoding_step_2.minimum))
                > acceFactorE2_[encoding])) {
            startE2_ = sampling_limits.kspace_encoding_step_2.minimum;
            endE2_   = sampling_limits.kspace_encoding_step_2.maximum;
        }

        long lenRO = endRO_ - startRO_ + 1;
        long lenE1 = endE1_ - startE1_ + 1;
        long lenE2 = endE2_ - startE2_ + 1;

        if (lenRO == RO && lenE1 == E1 && lenE2 == E2) {
            GDEBUG_CONDITION_STREAM(verbose, "lenRO == RO && lenE1 == E1 && lenE2 == E2");
            return recon_res;
        }

        // ----------------------------------------------------------
        // go to kspace
        // ----------------------------------------------------------
        hoNDArray<std::complex<float>> kspace_buf;

        if (E2 > 1) {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(recon_res.data, kspace_buf);
        } else {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(recon_res.data, kspace_buf);
        }

        // ----------------------------------------------------------
        // pf handling
        // ----------------------------------------------------------
        auto pf_res = this->perform_partial_fourier_handling(kspace_buf, startRO_, endRO_, startE1_, endE1_, startE2_, endE2_);

        // ----------------------------------------------------------
        // go back to image domain
        // ----------------------------------------------------------
        if (E2 > 1) {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(pf_res, recon_res.data);
        } else {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(pf_res, recon_res.data);
        }

        GDEBUG_CONDITION_STREAM(verbose, "GenericReconPartialFourierHandlingGadget::process(...) ends ... ");
        return std::move(recon_res);

    }

    // ----------------------------------------------------------------------------------------
}
