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
    AsymmetricEchoAdjustROGadget::AsymmetricEchoAdjustROGadget(const Core::Context& context, const Core::GadgetProperties& props): Core::ChannelGadget<Core::Acquisition>(context, props) {
        auto current_ismrmrd_header = (context.header);
        maxRO_.resize(current_ismrmrd_header.encoding.size());
        for (size_t e = 0; e < current_ismrmrd_header.encoding.size(); e++) {
            ISMRMRD::EncodingSpace e_space = current_ismrmrd_header.encoding[e].encodedSpace;
            maxRO_[e] = e_space.matrixSize.x;
            GDEBUG_STREAM("max RO for encoding space  " << e << " : " << maxRO_[e]);
        }
    }
    void AsymmetricEchoAdjustROGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) {
        for (auto [header, acq, traj] : in) {
            bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(header.flags);
            long long channels = (long long)header.active_channels;
            size_t samples = header.number_of_samples;
            size_t centre_column = header.center_sample;
            if (!is_noise) {
                unsigned int encoding_ref = header.encoding_space_ref;
                // adjust the center echo
                int az = addPrePostZeros(centre_column, samples);
                if (az != 0 && samples < maxRO_[encoding_ref]) {
                    std::vector<size_t> data_out_dims = *acq.get_dimensions();
                    data_out_dims[0] = maxRO_[encoding_ref];
                    hoNDArray<std::complex<float>> temp;
                    try {
                        temp = hoNDArray<std::complex<float>>(data_out_dims);
                    } catch (...) {
                        GDEBUG("Unable to create new data array for downsampled data\n");
                    }
                    temp.fill(0);
                    std::complex<float>* pM3 = temp.get_data_ptr();
                    std::complex<float>* pM2 = acq.get_data_ptr();
                    long long c;
                    size_t numOfBytes = sizeof(std::complex<float>) * samples;
                    if (az == 1) // pre zeros
                    {
                        //#pragma omp parallel for default(none) private(c) shared(channels, pM3, pM2, samples, numOfBytes)
                        for (c = 0; c < channels; c++) {
                            memcpy(pM3 + c * maxRO_[encoding_ref] + maxRO_[encoding_ref] - samples, pM2 + c * samples,
                                numOfBytes);
                        }
                        header.discard_pre = maxRO_[encoding_ref] - samples;
                        header.discard_post = 0;
                    }
                    if (az == 2) // post zeros
                    {
                        //#pragma omp parallel for default(none) private(c) shared(channels, pM3, pM2, samples, numOfBytes)
                        for (c = 0; c < channels; c++) {
                            memcpy(pM3 + c * maxRO_[encoding_ref], pM2 + c * samples, numOfBytes);
                        }
                        header.discard_pre = 0;
                        header.discard_post = maxRO_[encoding_ref] - samples;
                    }
                    acq = temp;
                    header.number_of_samples = (uint16_t)data_out_dims[0];
                    header.center_sample = header.number_of_samples / 2;
                }
            }
            out.push(Core::Acquisition{std::move(header), std::move(acq), std::move(traj)});
        }
    }
    GADGETRON_GADGET_EXPORT(AsymmetricEchoAdjustROGadget)
} // namespace Gadgetron
