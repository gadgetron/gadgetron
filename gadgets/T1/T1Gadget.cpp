//
// Created by dchansen on 3/9/20.
//
#include "demons_registration.h"

#include "PureGadget.h"
#include "mri_core_data.h"
#include "t1fit.h"

namespace Gadgetron {

    class T1Gadget : Core::PureGadget<Core::Image<float>, IsmrmrdImageArray> {

    public:
        T1Gadget(const Core::Context& context, const Core::GadgetProperties& properties)
            : Core::PureGadget<Core::Image<float>, IsmrmrdImageArray>(context, properties), TIs{*(context.header.sequenceParameters->TI)} {


        }

    private:
        Core::Image<float> process_function(IsmrmrdImageArray images) const final {

            return Core::Image<float>{};


        }


        std::vector<float> TIs;
    };
}