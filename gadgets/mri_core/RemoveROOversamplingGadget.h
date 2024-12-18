/**
    \brief  Removes readout oversampling
    \test   Tested by: simple_gre.cfg, generic_cartesian_cine_denoise.cfg, and others
*/

#pragma once

#include "Gadget.h"
#include "Node.h"
#include "hoNDArray.h"
#include "hoNDFFT.h"
#ifdef USE_OMP // TODO: Should this be removed? Its from the old version
#include "omp.h"
#endif // USE_OMP

namespace Gadgetron {
class RemoveROOversamplingGadget : public Core::ChannelGadget<mrd::Acquisition> {
  public:
    using Core::ChannelGadget<mrd::Acquisition>::ChannelGadget;
    RemoveROOversamplingGadget(const Core::Context& context, const Core::GadgetProperties& props);
    ~RemoveROOversamplingGadget() override = default;
    void process(Core::InputChannel<mrd::Acquisition>& input, Core::OutputChannel& output) override;

  protected:
    hoNDArray<std::complex<float>> fft_res_;
    hoNDArray<std::complex<float>> ifft_res_;

    hoNDArray<std::complex<float>> fft_buf_;
    hoNDArray<std::complex<float>> ifft_buf_;

    int encodeNx_;
    float encodeFOV_;
    int reconNx_;
    float reconFOV_;

    bool dowork_;    // if true the gadget performs the operation, otherwise, it just passes the data on
};
} // namespace Gadgetron
