#pragma once

#include "Node.h"
#include "hoNDArray.h"
#include "hoNDKLT.h"

#include <complex>
#include <map>

namespace Gadgetron {

  class PCACoilGadget : public Core::ChannelGadget<mrd::Acquisition>
  {
  public:
    PCACoilGadget(const Core::Context& context, const Core::GadgetProperties& props);
    ~PCACoilGadget() override;

    void process(Core::InputChannel<mrd::Acquisition>& input, Core::OutputChannel& output) override;

  protected:
    NODE_PROPERTY(uncombined_channels_by_name, std::string, "List of comma separated channels by name", "");
    // GADGET_PROPERTY(present_uncombined_channels, int, "Number of uncombined channels found", 0);

    void calculate_coefficients(int location);
    void do_pca(mrd::Acquisition& acq);

    std::vector<unsigned int> uncombined_channels_;

    //Map containing buffers, one for each location
    std::map< int, std::vector< mrd::Acquisition > > buffer_;

    //Keep track of whether we are buffering for a particular location
    std::map< int, bool> buffering_mode_;

    //Map for storing PCA coefficients for each location
    std::map<int, hoNDKLT<std::complex<float> >* > pca_coefficients_;

    const size_t max_buffered_profiles_ = 100;
    const size_t samples_to_use_ = 16;
  };
}