#include "SpiralToGenericGadget.h"
#include "ismrmrd/xml.h"
#include "vds.h"

#include <algorithm>
#include <vector>
#include "mri_core_girf_correction.h"

namespace Gadgetron {

    SpiralToGenericGadget::SpiralToGenericGadget()
            : samples_to_skip_start_(0), samples_to_skip_end_(0), prepared_(false) {
    }

    SpiralToGenericGadget::~SpiralToGenericGadget() {}

    int SpiralToGenericGadget::process_config(ACE_Message_Block *mb) {
        // Start parsing the ISMRMRD XML header
        //
        ISMRMRD::IsmrmrdHeader h;
        ISMRMRD::deserialize(mb->rd_ptr(), h);


        if (h.encoding.size() != 1) {
            GDEBUG("This Gadget only supports one encoding space\n");
            return GADGET_FAIL;
        }


        trajectory_parameters = Spiral::TrajectoryParameters(h);


        return GADGET_OK;
    }

    int SpiralToGenericGadget::
    process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
            GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2) {
        // Noise should have been consumed by the noise adjust, but just in case...
        //

        bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
        if (is_noise) {
            m1->release();
            return GADGET_OK;
        }

        // Delete previously attached trajectories
        if (m2->cont()) {
            m2->cont()->release();
        }

        // Compute hoNDArray of trajectory and weights at first pass
        //

        if (!prepared_) {
            prepare_trajectory(*m1->getObjectPtr());
            prepared_ = true;
        }

        auto samples_per_interleave = trajectory_and_weights.get_size(1);
        // Adjustments based in the incoming data
        //

        if (samples_to_skip_end_ == -1) {
            samples_to_skip_end_ = m1->getObjectPtr()->number_of_samples - samples_per_interleave;
            GDEBUG("Adjusting samples_to_skip_end_ = %d\n", samples_to_skip_end_);
        }

        // Define some utility variables
        //

        unsigned int samples_to_copy = m1->getObjectPtr()->number_of_samples - samples_to_skip_end_;
        unsigned int interleave = m1->getObjectPtr()->idx.kspace_encode_step_1;

        // Prepare for a new array continuation for the trajectory/weights of the incoming profile
        std::vector<size_t> trajectory_dimensions = {trajectory_and_weights.get_size(0),trajectory_and_weights.get_size(1)};

        hoNDArray<float> *traj_source = new hoNDArray<float>
                (trajectory_dimensions, trajectory_and_weights.get_data_ptr() + 3 * samples_per_interleave * interleave);

        // Make a new array as continuation of m1, and pass along
        //

        GadgetContainerMessage<hoNDArray<float> > *cont = new GadgetContainerMessage<hoNDArray<float> >();
        *(cont->getObjectPtr()) = *traj_source;
        m2->cont(cont);

        //We need to make sure that the trajectory dimensions are attached.
        m1->getObjectPtr()->trajectory_dimensions = 3;

        if (this->next()->putq(m1) < 0) {
            GDEBUG("Failed to put job on queue.\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    void SpiralToGenericGadget::prepare_trajectory(const ISMRMRD::AcquisitionHeader &acq_header) {

        hoNDArray<floatd2> trajectory;
        hoNDArray<float> weights;
        std::tie(trajectory,weights) = trajectory_parameters.calculate_trajectories_and_weight(acq_header);

        auto traj_dims = *trajectory.get_dimensions();
        std::vector<size_t> dims = {3};
        dims.insert(dims.end(),traj_dims.begin(),traj_dims.end());

        trajectory_and_weights = hoNDArray<float>(dims);

        size_t elements = weights.get_number_of_elements();
        for (size_t i = 0; i < elements; i++){
            trajectory_and_weights[i*3] = trajectory[i][0];
            trajectory_and_weights[i*3+1] = trajectory[i][1];
            trajectory_and_weights[i*3+2] = weights[i];
        }

    }

    GADGET_FACTORY_DECLARE(SpiralToGenericGadget)
}
