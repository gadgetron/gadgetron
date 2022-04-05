#include "PartialFourierAdjustROGadget.h"

namespace {
int addPrePostZeros(size_t centre_column, size_t samples) {
    // 1 : pre zeros
    // 2 : post zeros
    // 0 : no zeros
    if (2 * centre_column == samples) {
        return 0;
    }
    if (2 * centre_column < samples) {
        return 1;
    }
    if (2 * centre_column > samples) {
        return 2;
    }
    return 0;
}
} // namespace

namespace Gadgetron {

PartialFourierAdjustROGadget::PartialFourierAdjustROGadget(const Core::Context& context,
                                                           const Core::GadgetProperties& props)
    : Core::ChannelGadget<Core::Acquisition>(context, props) {

    auto h = (context.header);

    if (h.encoding.size() != 1) {
        GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
        GERROR("This partial fourier gadget only supports one encoding space\n");
        return;
    }

    ISMRMRD::EncodingSpaceType e_space = h.encoding[0].encodedSpace;
    maxRO_ = e_space.matrixSize.x;
    GDEBUG_STREAM("max RO : " << maxRO_);
    return;
}

void PartialFourierAdjustROGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) {
    for (auto [header, acq, traj] : in) {
        bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
        size_t channels = m1->getObjectPtr()->active_channels;
        size_t samples = m1->getObjectPtr()->number_of_samples;
        size_t centre_column = m1->getObjectPtr()->center_sample;

        if (!is_noise) {
            // adjust the center echo
            int az = addPrePostZeros(centre_column, samples);

            if (az != 0 && samples < maxRO_) {
                std::vector<size_t> data_out_dims = *acq.get_dimensions();
                data_out_dims[0] = maxRO_;
                hoNDArray<std::complex<float>> temp;
                try {
                    temp = hoNDArray<std::complex<float>>(data_out_dims);
                } catch (...) {
                    GDEBUG("Unable to create new data array for downsampled data\n");
                }
                temp.fill(0);
                std::complex<float>* pM3 = temp.get_data_ptr();
                std::complex<float>* pM2 = acq.get_data_ptr();
                size_t c;
                if (az == 1) // pre zeros
                {
                    for (c = 0; c < channels; c++) {
                        memcpy(pM3 + c * maxRO_ + maxRO_ - samples, pM2 + c * samples,
                               sizeof(std::complex<float>) * samples);
                    }
                }

                if (az == 2) // post zeros
                {
                    for (c = 0; c < channels; c++) {
                        memcpy(pM3 + c * maxRO_, pM2 + c * samples, sizeof(std::complex<float>) * samples);
                    }
                }

                acq = temp;
                header.number_of_samples = (uint16_t)data_out_dims[0];
            }
        }
        out.push(Core::Acquisition{header, std::move(acq), std::move(traj)});
    }
}
GADGETRON_GADGET_EXPORT(PartialFourierAdjustROGadget)
} // namespace Gadgetron
