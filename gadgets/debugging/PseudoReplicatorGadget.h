#pragma once

#include "Gadget.h"

namespace Gadgetron {

class PseudoReplicatorGadget : public Core::ChannelGadget<mrd::ReconData>
{
public:
	using Core::ChannelGadget<mrd::ReconData>::ChannelGadget;

	NODE_PROPERTY(repetitions, int, "Number of pseudoreplicas to produce", 10);
	PseudoReplicatorGadget(const Core::Context& context, const Core::GadgetProperties& props);

	void process(Core::InputChannel<mrd::ReconData>& input, Core::OutputChannel& output) override;
};

} /* namespace Gadgetron */
