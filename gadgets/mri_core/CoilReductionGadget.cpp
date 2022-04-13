#include "CoilReductionGadget.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

namespace Gadgetron {

CoilReductionGadget::CoilReductionGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : Core::ChannelGadget<Core::Acquisition>(context, props) {
    auto h = (context.header);
    coils_in_ =
        h.acquisitionSystemInformation->receiverChannels ? *h.acquisitionSystemInformation->receiverChannels : 128;

    std::string coil_mask_int = coil_mask;

    if (coil_mask_int.compare(std::string("")) == 0) {
        if (coils_out <= 0) {
            GERROR("Invalid number of output coils %d\n", coils_out);
            return;
        }
        coil_mask_ = std::vector<unsigned short>(coils_out, 1);
    } else {
        std::vector<std::string> chm;
        boost::split(chm, coil_mask_int, boost::is_any_of(" "));
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
        return;
    }

    coils_out_ = 0;
    for (size_t i = 0; i < coil_mask_.size(); i++) {
        if (coil_mask_[i])
            coils_out_++;
    }

    GDEBUG("Coil reduction from %d to %d\n", coils_in_, coils_out_);
}

void CoilReductionGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) {
    for (auto [header, acq, traj] : in) {

        std::vector<size_t> dims_out(2);
        dims_out[0] = header.number_of_samples;
        dims_out[1] = coils_out_;

        hoNDArray<std::complex<float>> m3 = hoNDArray<std::complex<float>>();

        try {
            m3.create(dims_out);
        } catch (std::runtime_error& err) {
            GEXCEPTION(err, "Unable to create storage for reduced dataset size\n");
            return;
        }

        std::complex<float>* s = acq.get_data_ptr();
        std::complex<float>* d = m3.get_data_ptr();

        size_t samples = header.number_of_samples;
        size_t coils_copied = 0;
        for (int c = 0; c < header.active_channels; c++) {
            if (c > coil_mask_.size()) {
                GERROR("Fatal error, too many coils for coil mask\n");
                return;
            }
            if (coil_mask_[c]) {
                memcpy(d + coils_copied * samples, s + c * samples, sizeof(std::complex<float>) * samples);
                coils_copied++;
            }
        }
        header.active_channels = coils_out_;
        out.push(Core::Acquisition{header, std::move(m3), std::move(traj)});
    }
}
GADGETRON_GADGET_EXPORT(CoilReductionGadget)
} // namespace Gadgetron
