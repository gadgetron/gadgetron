/**
    \brief  Performs basic FFT reconstruction on IsmrmrdReconData and passes along as IsmrmrdReconData
    \author Original: Souheil Inati
    \author ChannelGadget Conversion: Andrew Dupuis
    \test   Tested by: simple_gre.cfg, simple_gre_python_image_array_recon.cfg, and others
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

    class SimpleReconGadget : public Core::ChannelGadget<IsmrmrdReconData> {
    public:
        SimpleReconGadget(const Core::Context& context, const Core::GadgetProperties& props);
        void process(Core::InputChannel<IsmrmrdReconData>& input, Core::OutputChannel& out) override;

    protected:
        ISMRMRD::IsmrmrdHeader header;
        long long image_counter_;      

    };
}


