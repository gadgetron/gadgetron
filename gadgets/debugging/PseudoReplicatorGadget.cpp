#include "PseudoReplicatorGadget.h"
#include <random>

namespace Gadgetron {

PseudoReplicatorGadget::PseudoReplicatorGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : Core::ChannelGadget<mrd::ReconData>(context, props)
{
}

void PseudoReplicatorGadget::process(Core::InputChannel<mrd::ReconData>& input, Core::OutputChannel& output)
{
	std::mt19937 engine(5489UL);
	std::normal_distribution<float> distribution;

	for (auto reconData : input) {
		// First just send the normal data to obtain standard image
		output.push(reconData);

		for (int i = 0; i < repetitions; i++) {

			for (auto& rbit : reconData.buffers) {
				auto& data = rbit.data.data;
				auto dataptr = data.data();
				for (size_t k = 0; k < data.size(); k++) {
					dataptr[k] += std::complex<float>(distribution(engine), distribution(engine));
				}
			}

			GDEBUG("Sending out Pseudoreplica\n");
			output.push(reconData);
		}
	}
}

GADGETRON_GADGET_EXPORT(PseudoReplicatorGadget)
} // namespace Gadgetron
