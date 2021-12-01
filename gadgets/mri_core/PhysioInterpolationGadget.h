#ifndef PhysioInterpolationGadget_H
#define PhysioInterpolationGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

enum class PhysioInterpolationMode {
    separate,
    complete
};

enum class PhysioInterpolationMethod {
    Spline,
    BSpline
};

class PhysioInterpolationGadget : public Core::ChannelGadget<Core::Image<std::complex<float>>>
    {
    public:

        using Core::ChannelGadget<Core::Image<std::complex<float>>>::ChannelGadget;

        ~PhysioInterpolationGadget() override = default;


    protected:
        NODE_PROPERTY(physiology_time_index, int, "Physiology time index", 0);
        //GADGET_PROPERTY_LIMITS(mode, int, "Mode, 0=seperate series for each RR, 1=First complete RR only", 0, GadgetPropertyLimitsEnumeration, 0, 1);
        NODE_PROPERTY(mode, PhysioInterpolationMode, "Mode, 0=seperate series for each RR, 1=First complete RR only", PhysioInterpolationMode::separate);
        NODE_PROPERTY(phases, unsigned short, "Number of cardiac phases", 30);
        NODE_PROPERTY(first_beat_on_trigger, bool, "Indicates that acquisition was started on trigger", false);
        NODE_PROPERTY(interp_method, PhysioInterpolationMethod, "Interpolation method", PhysioInterpolationMethod::Spline);
        NODE_PROPERTY(time_stamp_resolution_, double, "Time stamp resolution in ms", 2.5);

      public:
        void process(Core::InputChannel<Core::Image<std::complex<float>>>& in, Core::OutputChannel& out) override;
};
}

#endif //PhysioInterpolationGadget_H
