/**
    \brief  Checks standard deviation of repeated sets of IsmrmrdReconData
    \test   Tested by: 
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
#include "hoNDArray_math.h"
#include "hoNDFFT.h"

namespace Gadgetron{

    class StandardDeviationGadget : public Core::ChannelGadget<IsmrmrdReconData> {
    public:
        StandardDeviationGadget(const Core::Context& context, const Core::GadgetProperties& props);
        void process(Core::InputChannel<IsmrmrdReconData>& input, Core::OutputChannel& out) override;

    protected:
        ISMRMRD::IsmrmrdHeader header;
        long long image_counter_;     
        NODE_PROPERTY(errorThreshold,float,"Maximum allowable mean of standard deviations", 0.000001);
    };
}