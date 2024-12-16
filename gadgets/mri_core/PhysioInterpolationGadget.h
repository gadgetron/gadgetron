#ifndef PhysioInterpolationGadget_H
#define PhysioInterpolationGadget_H

#include "Gadget.h"
#include "hoNDArray.h"

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


inline void from_string(const std::string& str, PhysioInterpolationMode& mode ){
    if (str == "separate" || str == "0" ) mode = PhysioInterpolationMode::separate;
    else if (str == "complete" || str == "1") mode = PhysioInterpolationMode::complete;
    else throw std::invalid_argument(str + " is not a valid PhysioInterpolationMode");
}

inline void from_string(const std::string& str, PhysioInterpolationMethod& method ){
    if (str == "Spline") method = PhysioInterpolationMethod::Spline;
    else if (str == "BSpline") method = PhysioInterpolationMethod::BSpline;
    else throw std::invalid_argument(str + " is not a valid input for PhysioInterpolationMethod");
}

class PhysioInterpolationGadget : public Core::ChannelGadget<mrd::Image<std::complex<float>>>
    {
    public:

        using Core::ChannelGadget<mrd::Image<std::complex<float>>>::ChannelGadget;

        ~PhysioInterpolationGadget() override = default;


    protected:
        NODE_PROPERTY(physiology_time_index, int, "Physiology time index", 0);
        NODE_PROPERTY(mode, PhysioInterpolationMode, "Mode, 0=separate series for each RR, 1=First complete RR only", PhysioInterpolationMode::separate);
        NODE_PROPERTY(phases, unsigned short, "Number of cardiac phases", 30);
        NODE_PROPERTY(first_beat_on_trigger, bool, "Indicates that acquisition was started on trigger", false);
        NODE_PROPERTY(interp_method, PhysioInterpolationMethod, "Interpolation method", PhysioInterpolationMethod::Spline);
        NODE_PROPERTY(time_stamp_resolution_, double, "Time stamp resolution in ms", 2.5);

      public:
        void process(Core::InputChannel<mrd::Image<std::complex<float>>>& in, Core::OutputChannel& out) override;

};
}

#endif //PhysioInterpolationGadget_H
