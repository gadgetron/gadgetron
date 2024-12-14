#include "AsymmetricEchoAdjustROGadget.h"

namespace {
    int addPrePostZeros(size_t centre_column, size_t samples) {
        // 1 : pre zeros
        // 2 : post zeros
        // 0 : no zeros
        if (2 * centre_column == samples || centre_column >= samples) {
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
    AsymmetricEchoAdjustROGadget::AsymmetricEchoAdjustROGadget(const Core::Context& context, const Core::GadgetProperties& props)
        : Core::ChannelGadget<mrd::Acquisition>(context, props)
    {
        auto current_mrd_header = (context.header);
        maxRO_.resize(current_mrd_header.encoding.size());
        for (size_t e = 0; e < current_mrd_header.encoding.size(); e++) {
            mrd::EncodingSpaceType e_space = current_mrd_header.encoding[e].encoded_space;
            maxRO_[e] = e_space.matrix_size.x;
            GDEBUG_STREAM("max RO for encoding space  " << e << " : " << maxRO_[e]);
        }
    }
    void AsymmetricEchoAdjustROGadget::process(Core::InputChannel<mrd::Acquisition>& in, Core::OutputChannel& out) {
        for (auto acq : in) {
            bool is_noise = acq.head.flags.HasFlags(mrd::AcquisitionFlags::kIsNoiseMeasurement);
            size_t channels = acq.Coils();
            size_t samples = acq.Samples();
            size_t centre_column = acq.head.center_sample.value_or(samples / 2);
            if (!is_noise) {
                unsigned int encoding_ref = acq.head.encoding_space_ref.value_or(0);
                // adjust the center echo
                int az = addPrePostZeros(centre_column, samples);
                if (az != 0 && samples < maxRO_[encoding_ref]) {
                    auto original_data = acq.data;

                    acq.data.create(maxRO_[encoding_ref], channels);
                    acq.data.fill(0);

                    auto pAdjusted = acq.data.data();
                    auto pOrig = original_data.data();

                    size_t numOfBytes = sizeof(std::complex<float>) * samples;
                    if (az == 1) // pre zeros
                    {
                        //#pragma omp parallel for default(none) private(c) shared(channels, pAdjusted, pOrig, samples, numOfBytes)
                        for (size_t c = 0; c < channels; c++) {
                            memcpy(pAdjusted + c * maxRO_[encoding_ref] + maxRO_[encoding_ref] - samples, pOrig + c * samples,
                                numOfBytes);
                        }
                        acq.head.discard_pre = maxRO_[encoding_ref] - samples;
                        acq.head.discard_post = 0;
                    }
                    if (az == 2) // post zeros
                    {
                        //#pragma omp parallel for default(none) private(c) shared(channels, pAdjusted, pOrig, samples, numOfBytes)
                        for (size_t c = 0; c < channels; c++) {
                            memcpy(pAdjusted + c * maxRO_[encoding_ref], pOrig + c * samples, numOfBytes);
                        }
                        acq.head.discard_pre = 0;
                        acq.head.discard_post = maxRO_[encoding_ref] - samples;
                    }
                    acq.head.center_sample = acq.Samples() / 2;
                }
            }
            out.push(std::move(acq));
        }
    }
    GADGETRON_GADGET_EXPORT(AsymmetricEchoAdjustROGadget)
} // namespace Gadgetron
