#pragma once

#include <ismrmrd/xml.h>
#include "log.h"
#include "Gadget.h"
#include "vds.h"

namespace Gadgetron {
namespace Spiral {

    class TrajectoryParameters {
    public:
        TrajectoryParameters() = default;
        TrajectoryParameters(const ISMRMRD::IsmrmrdHeader &h);

        std::pair<hoNDArray<floatd2>, hoNDArray<float>>
        calculate_trajectories_and_weight(const ISMRMRD::AcquisitionHeader &acq_header);

    private:
        Core::optional<hoNDArray<std::complex<float>>> girf_kernel;
        float girf_sampling_time_us;
        long Tsamp_ns_;
        long Nints_;
        double gmax_;
        double smax_;
        double krmax_;
        double fov_;
        float TE_;

        hoNDArray<floatd2> correct_gradients(const hoNDArray<floatd2> &gradients, float grad_samp_us,
                                             float girf_samp_us, const float *read_dir, const float *phase_dir,
                                             const float *slice_dir);


    };
}
}
