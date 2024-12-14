#include "SpiralToGenericGadget.h"
#include "../../toolboxes/mri/spiral/vds.h"

#include <algorithm>
#include <vector>
#include "mri_core_girf_correction.h"

namespace Gadgetron {

    SpiralToGenericGadget::SpiralToGenericGadget()
            : samples_to_skip_start_(0), samples_to_skip_end_(0), prepared_(false) {
    }

    SpiralToGenericGadget::~SpiralToGenericGadget() {}

    int SpiralToGenericGadget::process_config(const mrd::Header& header) {
        if (header.encoding.size() != 1) {
            GDEBUG("This Gadget only supports one encoding space\n");
            return GADGET_FAIL;
        }

        trajectory_parameters_ = Spiral::TrajectoryParameters(header);

        return GADGET_OK;
    }

    int SpiralToGenericGadget::process(GadgetContainerMessage<mrd::Acquisition> *m1)
    {
        mrd::Acquisition& acq = *m1->getObjectPtr();
        // Noise should have been consumed by the noise adjust, but just in case...
        //

        bool is_noise = acq.head.flags.HasFlags(mrd::AcquisitionFlags::kIsNoiseMeasurement);
        if (is_noise) {
            m1->release();
            return GADGET_OK;
        }

        // Compute hoNDArray of trajectory and weights at first pass
        //
        if (!prepared_) {
            prepare_trajectory(acq);
            prepared_ = true;
        }

        auto samples_per_interleave = trajectory_and_weights_.get_size(1);
        // Adjustments based in the incoming data
        //
        if (samples_to_skip_end_ == -1) {
            samples_to_skip_end_ = m1->getObjectPtr()->Samples() - samples_per_interleave;
            GDEBUG("Adjusting samples_to_skip_end_ = %d\n", samples_to_skip_end_);
        }

        // Define some utility variables
        //
        unsigned int samples_to_copy = acq.Samples() - samples_to_skip_end_;
        unsigned int interleave = acq.head.idx.kspace_encode_step_1.value_or(0);

        // Prepare for a new array continuation for the trajectory/weights of the incoming profile
        std::vector<size_t> trajectory_dimensions = {trajectory_and_weights_.get_size(0),trajectory_and_weights_.get_size(1)};

        hoNDArray<float> *traj_source = new hoNDArray<float>
                (trajectory_dimensions, trajectory_and_weights_.get_data_ptr() + 3 * samples_per_interleave * interleave);

        // Attach new trajectory
        acq.trajectory = *traj_source;

        if (this->next()->putq(m1) < 0) {
            GDEBUG("Failed to put job on queue.\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    void SpiralToGenericGadget::prepare_trajectory(const mrd::Acquisition &acq) {

        hoNDArray<floatd2> trajectory;
        hoNDArray<float> weights;
        std::tie(trajectory,weights) = trajectory_parameters_.calculate_trajectories_and_weight(acq);

        auto traj_dims = trajectory.get_dimensions();
        std::vector<size_t> dims = {3};
        dims.insert(dims.end(),traj_dims.begin(),traj_dims.end());

        trajectory_and_weights_ = hoNDArray<float>(dims);

        size_t elements = weights.get_number_of_elements();
        for (size_t i = 0; i < elements; i++){
            trajectory_and_weights_[i*3] = trajectory[i][0];
            trajectory_and_weights_[i*3+1] = trajectory[i][1];
            trajectory_and_weights_[i*3+2] = weights[i];
        }

    }

    GADGET_FACTORY_DECLARE(SpiralToGenericGadget)
}
