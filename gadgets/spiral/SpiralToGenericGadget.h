#ifndef SpiralToGenericGadget_H
#define SpiralToGenericGadget_H
#pragma once

#include "../../toolboxes/mri/spiral/TrajectoryParameters.h"
#include "Gadget.h"
#include "hoNDArray.h"

#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>

namespace Gadgetron {

    class SpiralToGenericGadget :
            public Gadget1<mrd::Acquisition> {

    public:
        SpiralToGenericGadget();

        virtual ~SpiralToGenericGadget();

    protected:

        virtual int process_config(const mrd::Header& header);

        virtual int process(GadgetContainerMessage<mrd::Acquisition> *m1);

    private:

        bool prepared_;

        int samples_to_skip_start_;
        int samples_to_skip_end_;
        hoNDArray<float> trajectory_and_weights_;
        Spiral::TrajectoryParameters trajectory_parameters_;
        void prepare_trajectory(const mrd::Acquisition &acq);

    };
}
#endif //SpiralToGenericGadget_H
