/**
    \brief  Appends defined number of noisy pseudoreplicas to incoming IsmrmrdReconData
    \author Original: David Hansen
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Untested
*/

#pragma once
#include "Node.h"
#include "mri_core_data.h"
#include <random>

namespace Gadgetron{
    class PseudoReplicatorGadget : public Core::ChannelGadget<IsmrmrdReconData> {
    public:
      using Core::ChannelGadget<IsmrmrdReconData>::ChannelGadget;
      void process(Core::InputChannel<IsmrmrdReconData>& input, Core::OutputChannel& out) override;
      NODE_PROPERTY(repetitions,int,"Number of pseudoreplicas to produce",10);
    };
}


