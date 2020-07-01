#ifndef SpiralToGenericGadget_H
#define SpiralToGenericGadget_H
#pragma once

#include "../../toolboxes/mri/spiral/TrajectoryParameters.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "gadgetron_spiral_export.h"
#include "hoNDArray.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <ismrmrd/xml.h>
#include <boost/optional.hpp>

namespace Gadgetron {

    class EXPORTGADGETS_SPIRAL SpiralToGenericGadget :
            public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float> > > {

    public:
        GADGET_DECLARE(SpiralToGenericGadget);

        SpiralToGenericGadget();

        virtual ~SpiralToGenericGadget();

    protected:

        virtual int process_config(ACE_Message_Block *mb);

        virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
                            GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2);

    private:

        bool prepared_;

        int samples_to_skip_start_;
        int samples_to_skip_end_;
        hoNDArray<float> trajectory_and_weights;
        Spiral::TrajectoryParameters trajectory_parameters;
        void prepare_trajectory(const ISMRMRD::AcquisitionHeader &acq_header);

    };
}
#endif //SpiralToGenericGadget_H
