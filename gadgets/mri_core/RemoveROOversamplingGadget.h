/**
    \brief  Removes readout oversampling
    \author Original: Thomas Sangild Sorensen
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Tested by: simple_gre.cfg, generic_cartesian_cine_denoise.cfg, and others
*/

#pragma once

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "Node.h"
#include "Types.h"
#include "hoNDArray.h"

namespace Gadgetron {
class RemoveROOversamplingGadget : public Core::ChannelGadget<Core::Acquisition> {
  public:
    using Core::ChannelGadget<Core::Acquisition>::ChannelGadget;
    RemoveROOversamplingGadget(const Core::Context& context, const Core::GadgetProperties& props);
    ~RemoveROOversamplingGadget() override = default;
    void process(Core::InputChannel<Core::Acquisition>& input, Core::OutputChannel& output) override;

  protected:
    hoNDArray<std::complex<float>> fft_res_;
    hoNDArray<std::complex<float>> ifft_res_;

    hoNDArray<std::complex<float>> fft_buf_;
    hoNDArray<std::complex<float>> ifft_buf_;

    int encodeNx_;
    float encodeFOV_;
    int reconNx_;
    float reconFOV_;

    // if true the gadget performs the operation
    // otherwise, it just passes the data on
    bool dowork_;
};
} // namespace Gadgetron
