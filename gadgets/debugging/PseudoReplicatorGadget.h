/**
    \brief  Appends defined number of noisy pseudoreplicas to incoming IsmrmrdReconData
    \test   pseudoreplica.cfg
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
      NODE_PROPERTY(seed,unsigned long,"Random number generator seed",5489UL);
    };
}