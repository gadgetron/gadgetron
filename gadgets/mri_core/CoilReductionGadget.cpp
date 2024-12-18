#include "CoilReductionGadget.h"

namespace Gadgetron {

CoilReductionGadget::CoilReductionGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : Core::ChannelGadget<mrd::Acquisition>(context, props)
{
    auto h = context.header;
    coils_in_ = h.acquisition_system_information->receiver_channels.value_or(128);

    if (coil_mask.empty()) {
        if (coils_out <= 0) {
            GERROR("Invalid number of output coils %d\n", coils_out);
        }
        coil_mask_ = std::vector<unsigned short>(coils_out, 1);
    }
    else {
        std::vector<std::string> chm;
        boost::split(chm, coil_mask, boost::is_any_of(" "));
        for (size_t i = 0; i < chm.size(); i++) {
            std::string ch = boost::algorithm::trim_copy(chm[i]);
            if (ch.size() > 0) {
                size_t mv = static_cast<size_t>(std::stoi(ch));
                GDEBUG("Coil mask value: %d\n", mv);
                if (mv > 0) {
                    coil_mask_.push_back(1);
                } else {
                    coil_mask_.push_back(0);
                }
            }
        }
    }

    while (coil_mask_.size() < coils_in_)
        coil_mask_.push_back(0);
    while (coil_mask_.size() > coils_in_)
        coil_mask_.pop_back();

    if (coil_mask_.size() != coils_in_) {
        GERROR("Error configuring coils for coil reduction\n");
    }

    coils_out_ = 0;
    for (size_t i = 0; i < coil_mask_.size(); i++) {
        if (coil_mask_[i])
            coils_out_++;
    }

    GDEBUG("Coil reduction from %d to %d\n", coils_in_, coils_out_);
}

void CoilReductionGadget::process(Core::InputChannel<mrd::Acquisition>& in, Core::OutputChannel& out) {
    for (auto acq : in) {
        if (acq.Coils() == coils_out_) {
            // No need to do anything
            out.push(std::move(acq));
            continue;
        }

        if (acq.Coils() > coil_mask_.size()) {
            GERROR("Fatal error, too many coils for coil mask\n");
            continue;
        }

        auto nsamples = acq.Samples();
        auto nchannels = acq.Coils();

        std::vector<size_t> dims_out(2);
        dims_out[0] = nsamples;
        dims_out[1] = coils_out_;
        hoNDArray<std::complex<float>> reduced(dims_out);

        acq.head.channel_order.resize(coils_out_);
        auto s = acq.data.data();
        auto d = reduced.data();

        size_t coils_copied = 0;
        for (size_t c = 0; c < nchannels; c++) {
            if (coil_mask_[c]) {
                memcpy(d + coils_copied * nsamples, s + c * nsamples, sizeof(std::complex<float>) * nsamples);
                acq.head.channel_order[coils_copied] = coils_copied;
                coils_copied++;
            }
        }

        acq.data = std::move(reduced);

        out.push(std::move(acq));
    }
}
GADGETRON_GADGET_EXPORT(CoilReductionGadget)
}
