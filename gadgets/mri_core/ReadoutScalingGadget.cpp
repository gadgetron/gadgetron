#include "ReadoutScalingGadget.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron
{
    ReadoutScalingGadget::ReadoutScalingGadget(const Core::Context& context, const Core::GadgetProperties& props): Core::ChannelGadget<Core::Acquisition>(context, props)
    {
        auto current_ismrmrd_header = (context.header);
    }

    void ReadoutScalingGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out)
    {
        for (auto [header, acq, traj] : in)
        {
            bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(header.flags);
            long long channels = (long long)header.active_channels;
            size_t samples = header.number_of_samples;
            size_t centre_column = header.center_sample;
            if (!is_noise && noise_scaling_factor>0)
            {
                float scale_ratio = this->noise_scaling_factor * std::sqrt(512.0/(2.0*(samples-centre_column)));
                if (header.scan_counter<10)
                {
                    GDEBUG_STREAM("Applying the scaling " << scale_ratio << " to readout " << header.scan_counter);
                }
                Gadgetron::scal(scale_ratio, acq);
            }

            out.push(Core::Acquisition{std::move(header), std::move(acq), std::move(traj)});
        }
    }

    GADGETRON_GADGET_EXPORT(ReadoutScalingGadget)

} // namespace Gadgetron
