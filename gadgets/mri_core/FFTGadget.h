/**
    \brief  Breaks IsmrmrdReconData into separate images and performs FFT
    \author Original: Thomas Sangild Sorensen
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Tested by: epi_2d.cfg
*/

#pragma once
#include "Node.h"
#include "gadgetron_mricore_export.h"
#include "hoNDArray.h"

#include "mri_core_acquisition_bucket.h"
#include "mri_core_data.h"
#include <complex>
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

namespace Gadgetron{

    class FFTGadget : public Core::ChannelGadget<IsmrmrdReconData> {
    public:
        FFTGadget(const Core::Context& context, const Core::GadgetProperties& props);
        void process(Core::InputChannel<IsmrmrdReconData>& input, Core::OutputChannel& out) override;

    protected:
        ISMRMRD::IsmrmrdHeader header;
        long long image_counter_;      

    };
}


