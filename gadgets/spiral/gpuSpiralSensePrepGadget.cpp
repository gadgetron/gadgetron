#include "gpuSpiralSensePrepGadget.h"
#include "../../toolboxes/mri/spiral/vds.h"
#include "GPUTimer.h"
#include "GenericReconJob.h"
#include "b1_map.h"
#include "check_CUDA.h"
#include "cuNDArray_reductions.h"
#include "cuNDArray_utils.h"
#include "hoNDArray_fileio.h"
#include "ismrmrd/xml.h"
#include "mri_core_girf_correction.h"
#include "vector_td.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"
#include <armadillo>

#include <algorithm>
#include <vector>
#include <boost/range/algorithm/copy.hpp>

namespace Gadgetron {

    gpuSpiralSensePrepGadget::gpuSpiralSensePrepGadget()
            : samples_to_skip_start_(0), samples_to_skip_end_(-1), samples_per_interleave_(0), prepared_(false),
              use_multiframe_grouping_(false), acceleration_factor_(0) {
    }

    gpuSpiralSensePrepGadget::~gpuSpiralSensePrepGadget() {
        for (auto& buffer : this->buffer_) {
            for (auto& m1m2 : buffer){
                std::get<0>(m1m2)->release();
            }
        }

    }

    int gpuSpiralSensePrepGadget::process_config(ACE_Message_Block *mb) {

        int number_of_devices = 0;
        if (cudaGetDeviceCount(&number_of_devices) != cudaSuccess) {
            GDEBUG("Error: unable to query number of CUDA devices.\n");
            return GADGET_FAIL;
        }

        if (number_of_devices == 0) {
            GDEBUG("Error: No available CUDA devices.\n");
            return GADGET_FAIL;
        }

        device_number_ = deviceno.value();

        if (device_number_ >= number_of_devices) {
            GDEBUG("Adjusting device number from %d to %d\n", device_number_, (device_number_ % number_of_devices));
            device_number_ = (device_number_ % number_of_devices);
        }

        if (cudaSetDevice(device_number_) != cudaSuccess) {
            GDEBUG("Error: unable to set CUDA device.\n");
            return GADGET_FAIL;
        }

        cudaDeviceProp deviceProp;
        if (cudaGetDeviceProperties(&deviceProp, device_number_) != cudaSuccess) {
            GDEBUG("Error: unable to query device properties.\n");
            return GADGET_FAIL;
        }

        unsigned int warp_size = deviceProp.warpSize;

        propagate_csm_from_set_ = propagate_csm_from_set.value();

        if (propagate_csm_from_set_ > 0) {
            GDEBUG("Currently, only set 0 can propagate coil sensitivity maps. Set %d was specified.\n",
                   propagate_csm_from_set_);
            return GADGET_FAIL;
        }

        if (propagate_csm_from_set_ >= 0) {
            GDEBUG("Propagating csm from set %d to all sets\n", propagate_csm_from_set_);
        }

        buffer_using_solver_ = buffer_using_solver.value();
        use_multiframe_grouping_ = use_multiframe_grouping.value();

        if (buffer_using_solver_ && !use_multiframe_grouping_) {
            GDEBUG("Enabling 'buffer_using_solver' requires also enabling 'use_multiframe_grouping'.\n");
            return GADGET_FAIL;
        }

        // Start parsing the ISMRMRD XML header
        //

        ISMRMRD::IsmrmrdHeader h;
        ISMRMRD::deserialize(mb->rd_ptr(), h);


        if (h.encoding.size() != 1) {
            GDEBUG("This Gadget only supports one encoding space\n");
            return GADGET_FAIL;
        }

        // Get the encoding space and trajectory description

        ISMRMRD::TrajectoryDescription traj_desc;


        // Determine reconstruction matrix sizes
        //

        ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;

        kernel_width_ = buffer_convolution_kernel_width.value();
        oversampling_factor_ = buffer_convolution_oversampling_factor.value();

        image_dimensions_recon_.push_back(
                ((static_cast<unsigned int>(std::ceil(e_space.matrixSize.x * reconstruction_os_factor_x.value())) +
                  warp_size - 1) / warp_size) * warp_size);
        image_dimensions_recon_.push_back(
                ((static_cast<unsigned int>(std::ceil(e_space.matrixSize.y * reconstruction_os_factor_y.value())) +
                  warp_size - 1) / warp_size) * warp_size);

        image_dimensions_recon_os_ = uint64d2
                (((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[0] * oversampling_factor_)) + warp_size -
                   1) / warp_size) * warp_size,
                 ((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[1] * oversampling_factor_)) + warp_size -
                   1) / warp_size) * warp_size);

        // In case the warp_size constraint kicked in
        oversampling_factor_ = float(image_dimensions_recon_os_[0]) / float(image_dimensions_recon_[0]);



        ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
        ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

        fov_vec_.push_back(r_space.fieldOfView_mm.x);
        fov_vec_.push_back(r_space.fieldOfView_mm.y);
        fov_vec_.push_back(r_space.fieldOfView_mm.z);

        slices_ = e_limits.slice ? e_limits.slice->maximum + 1 : 1;
        sets_ = e_limits.set ? e_limits.set->maximum + 1 : 1;

        buffer_ = decltype(buffer_)(slices_ * sets_);

        image_headers_queue_ = decltype(image_headers_queue_)(slices_ * sets_);

        GDEBUG("recon matrix_size_x    : %d\n", image_dimensions_recon_[0]);
        GDEBUG("recon matrix_size_y    : %d\n", image_dimensions_recon_[1]);


        trajectoryParameters = Spiral::TrajectoryParameters(h);

        return GADGET_OK;
    }



    int gpuSpiralSensePrepGadget::
    process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
            GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2) {
        // Noise should have been consumed by the noise adjust, but just in case...
        //
        auto &header = *m1->getObjectPtr();
        bool is_noise = header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
        if (is_noise) {
            m1->release();
            return GADGET_OK;
        }


        if (!prepared_) {
            prepare_nfft(header);
            prepared_ = true;
        }

        // Allocate host data buffer if it is NULL
        setup_buffers(header);

        // Define some utility variables
        unsigned int samples_to_copy = header.number_of_samples - samples_to_skip_end_;
        unsigned int interleave = header.idx.kspace_encode_step_1;
        unsigned int slice = header.idx.slice;
        unsigned int set = header.idx.set;
        unsigned int samples_per_channel = host_data_buffer_[set * slices_ + slice].get_size(0);

        // Some book-keeping to keep track of the frame count
        interleaves_counter_singleframe_[set * slices_ + slice]++;
        interleaves_counter_multiframe_[set * slices_ + slice]++;

        // Duplicate the profile to avoid double deletion in case problems are encountered below.
        // Enqueue profile until all profiles for the reconstruction have been received.
        buffer_[set * slices_ + slice].push_back(std::make_pair(m1,m2));

        // Copy profile into the accumulation buffer for csm/regularization estimation

        if (samples_to_skip_end_ == -1) {
            samples_to_skip_end_ = header.number_of_samples - samples_per_interleave_;
            GDEBUG("Adjusting samples_to_skip_end_ = %d\n", samples_to_skip_end_);
        }

        {
            auto data_ptr = host_data_buffer_[set * slices_ + slice].get_data_ptr();

            std::complex<float> *profile_ptr = m2->getObjectPtr()->get_data_ptr();
            for (unsigned int c = 0; c < header.active_channels; c++) {
                memcpy(data_ptr + c * samples_per_channel + interleave * samples_to_copy,
                       profile_ptr + c * header.number_of_samples, samples_to_copy * sizeof(std::complex<float>));
            }
        }

        // Have we received sufficient data for a new frame?
        //

        bool is_last_scan_in_slice = header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

        if (is_last_scan_in_slice) {

            // This was the final profile of a frame
            //

            if (interleaves_ % interleaves_counter_singleframe_[set * slices_ + slice]) {
                GDEBUG("Unexpected number of interleaves encountered in frame\n");
                return GADGET_FAIL;
            }

            // Has the acceleration factor changed?
            //

            if (acceleration_factor_ != interleaves_ / interleaves_counter_singleframe_[set * slices_ + slice]) {
                change_acceleration_factor(header);
            }

            // Prepare an image header for this frame
            //

            image_headers_queue_[set * slices_ + slice].emplace_back(make_image_header(header));
            // Check if it is time to reconstruct.
            // I.e. prepare and pass a Sense job downstream...
            //
            if (!use_multiframe_grouping_ ||
                (use_multiframe_grouping_ && interleaves_counter_multiframe_[set * slices_ + slice] == interleaves_)) {

                unsigned int num_coils = header.active_channels;

                //

                // Compute coil images from the fully sampled data buffer
                cuNDArray<float_complext> reg_image = make_reg_image(host_data_buffer_[set * slices_ + slice], set,
                                                                     num_coils);
                boost::shared_ptr<hoNDArray<float_complext> > reg_host = reg_image.to_host();
                boost::shared_ptr<hoNDArray<float_complext> > csm_host = csm_->to_host();

                auto queue_data = get_data_from_queues(set, slice, num_coils);
                GadgetContainerMessage<GenericReconJob> *m4 = new GadgetContainerMessage<GenericReconJob>();

                m4->getObjectPtr()->dat_host_ = std::get<0>(queue_data);
                m4->getObjectPtr()->csm_host_ = csm_host;
                m4->getObjectPtr()->reg_host_ = reg_host;
                m4->getObjectPtr()->tra_host_ = std::get<1>(queue_data);
                m4->getObjectPtr()->dcw_host_ = std::get<2>(queue_data);

                // Pull the image headers out of the queue
                //

                long frames_per_reconstruction = (use_multiframe_grouping_) ? acceleration_factor_ : 1;

                if (image_headers_queue_[set * slices_ + slice].size() != frames_per_reconstruction) {
                    m4->release();
                    GDEBUG("Unexpected size of image header queue: %d, %d\n",
                           image_headers_queue_[set * slices_ + slice].size(), frames_per_reconstruction);
                    return GADGET_FAIL;
                }

                m4->getObjectPtr()->image_headers_ =
                        boost::shared_array<ISMRMRD::ImageHeader>(new ISMRMRD::ImageHeader[frames_per_reconstruction]);

                for (unsigned int i = 0; i < frames_per_reconstruction; i++) {
                    auto m = image_headers_queue_[set*slices_+slice][i];
                    m4->getObjectPtr()->image_headers_[i] = m;
                }
                image_headers_queue_[set* slices_ + slice].clear();

                // The Sense Job needs an image header as well.
                // Let us just copy the initial one...

                GadgetContainerMessage<ISMRMRD::ImageHeader> *m3 = new GadgetContainerMessage<ISMRMRD::ImageHeader>;
                *m3->getObjectPtr() = m4->getObjectPtr()->image_headers_[0];
                m3->cont(m4);

                if (this->next()->putq(m3) < 0) {
                    GDEBUG("Failed to put job on queue.\n");
                    m3->release();
                    return GADGET_FAIL;
                }
                interleaves_counter_multiframe_[set * slices_ + slice] = 0;
            }
            interleaves_counter_singleframe_[set * slices_ + slice] = 0;
        }
        return GADGET_OK;
    }

    void gpuSpiralSensePrepGadget::setup_buffers(const ISMRMRD::AcquisitionHeader &header) {
        if (host_data_buffer_.empty()) {

            std::vector<size_t> data_dimensions = {size_t(samples_per_interleave_ * interleaves_),
                                                   header.active_channels};
            host_data_buffer_ = std::vector<hoNDArray<float_complext>>(slices_ * sets_,
                                                                       hoNDArray<float_complext>(data_dimensions));

            for (auto &buffer : host_data_buffer_)
                std::fill(buffer.begin(), buffer.end(), 0);
        }

        // Allocate various counters if they are NULL
//

        if (image_counter_.empty()) {
            image_counter_ = std::vector<long>(slices_ * sets_, 0);
        }

        if (interleaves_counter_singleframe_.empty()) {
            interleaves_counter_singleframe_ = std::vector<long>(slices_ * sets_, 0);
        }

        if (interleaves_counter_multiframe_.empty()) {
            interleaves_counter_multiframe_ = std::vector<long>(slices_ * sets_, 0);
        }
    }

    cuNDArray<float_complext>
    gpuSpiralSensePrepGadget::make_reg_image(const hoNDArray<float_complext> &buffer, size_t set, size_t num_coils) {
        std::vector<size_t> image_dims{image_dimensions_recon_[0], image_dimensions_recon_[1], num_coils};
        cuNDArray<float_complext> image(image_dims);
        cuNDArray<float_complext> data(buffer);

        nfft_plan_->compute(data, image, dcw_buffer_.get(), NFFT_comp_mode::BACKWARDS_NC2C);

        // Check if we need to compute a new csm
        if (propagate_csm_from_set_ < 0 || propagate_csm_from_set_ == set || !csm_) {
            csm_ = boost::make_shared<cuNDArray<float_complext>>(estimate_b1_map<float, 2>(image)); // Estimates csm
        }
        E_->set_csm(csm_);

        // Compute regularization using basic coil combination
//

        image_dims.pop_back();
        cuNDArray<float_complext> reg_image(image_dims);
        E_->mult_csm_conj_sum(&image, &reg_image);

        if (buffer_using_solver_) {

            // Compute regularization using cg solver
            //

            // Define preconditioning weights
            boost::shared_ptr<cuNDArray<float> > _precon_weights = sum(abs_square(csm_.get()).get(), 2);
            reciprocal_sqrt_inplace(_precon_weights.get());
            boost::shared_ptr<cuNDArray<float_complext> > precon_weights = real_to_complex<float_complext>(
                    _precon_weights.get());
            _precon_weights.reset();
            D_->set_weights(precon_weights);

            // Solve from the plain coil combination
            reg_image = *cg_.solve_from_rhs(&reg_image);
        }
        return reg_image;
    }

    void gpuSpiralSensePrepGadget::change_acceleration_factor(const ISMRMRD::AcquisitionHeader &header) {
        GDEBUG("Change of acceleration factor detected\n");
        acceleration_factor_ = interleaves_/ interleaves_counter_singleframe_[header.idx.set * slices_ + header.idx.slice];

        // The encoding operator needs to have its domain/codomain dimensions set accordingly
        if (buffer_using_solver_) {

            std::vector<size_t> domain_dims = image_dimensions_recon_;

            std::vector<size_t> codomain_dims = host_traj_.get_dimensions();
            codomain_dims.push_back(header.active_channels);

            E_->set_domain_dimensions(domain_dims);
            E_->set_codomain_dimensions(codomain_dims);

            cuNDArray<floatd2> traj(host_traj_);
            E_->preprocess(&traj);
        }
    }

    ISMRMRD::ImageHeader
    gpuSpiralSensePrepGadget::make_image_header(const ISMRMRD::AcquisitionHeader &acq_header) {

        auto header = ISMRMRD::ImageHeader();

        auto set = acq_header.idx.set;
        auto slice = acq_header.idx.slice;

        header.version = acq_header.version;

        header.matrix_size[0] = image_dimensions_recon_[0];
        header.matrix_size[1] = image_dimensions_recon_[1];
        header.matrix_size[2] = acceleration_factor_;

        boost::copy(fov_vec_,header.field_of_view);

        header.channels = acq_header.active_channels;
        header.slice = acq_header.idx.slice;
        header.set = acq_header.idx.set;


        header.acquisition_time_stamp = acq_header.acquisition_time_stamp;

        boost::copy(acq_header.physiology_time_stamp,header.physiology_time_stamp);
        boost::copy(acq_header.position,header.position);
        boost::copy(acq_header.read_dir,header.read_dir);
        boost::copy(acq_header.phase_dir,header.phase_dir);
        boost::copy(acq_header.slice_dir,header.slice_dir);
        boost::copy(acq_header.patient_table_position,header.patient_table_position);


        header.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
        header.image_index = image_counter_[set * slices_ + slice]++;
        header.image_series_index = set * slices_ + slice;

        // Enqueue header until we are ready to assemble a Sense job
        //

        return header;
    }

    void gpuSpiralSensePrepGadget::prepare_nfft( const ISMRMRD::AcquisitionHeader& acq_header) {



        // Setup the NFFT plan
        //

        std::tie(host_traj_,host_weights_) = trajectoryParameters.calculate_trajectories_and_weight(acq_header);
        interleaves_ = host_traj_.get_size(1);
        samples_per_interleave_ = host_traj_.get_size(0);


        host_traj_.reshape({int64_t(samples_per_interleave_*interleaves_)});
        host_weights_.reshape({int64_t(samples_per_interleave_*interleaves_)});

        cuNDArray<floatd2> traj(host_traj_);
        dcw_buffer_ = boost::make_shared<cuNDArray<float>>(host_weights_);

        nfft_plan_ = NFFT<cuNDArray,float,2>::make_plan(from_std_vector<size_t, 2>(image_dimensions_recon_), image_dimensions_recon_os_,
                         kernel_width_);
        nfft_plan_->preprocess(traj, NFFT_prep_mode::NC2C);



        // Setup the non-Cartesian Sense encoding operator
        //



        // Setup cg solver if the csm/regularization image is to be based hereon
        //

        E_ = boost::shared_ptr<cuNonCartesianSenseOperator<float, 2> >(new cuNonCartesianSenseOperator<float, 2>);
        E_->setup(from_std_vector<size_t, 2>(image_dimensions_recon_), image_dimensions_recon_os_, kernel_width_);


        if (buffer_using_solver_) {

            E_->set_dcw(sqrt(dcw_buffer_.get()));

            D_ = boost::shared_ptr<cuCgPreconditioner<float_complext> >(new cuCgPreconditioner<float_complext>());
            cg_.set_encoding_operator(E_);
            cg_.set_preconditioner(D_);
            cg_.set_max_iterations(2);
            cg_.set_tc_tolerance(1e-6);
            cg_.set_output_mode(decltype(cg_)::OUTPUT_SILENT);
        }

    }

    GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *
    gpuSpiralSensePrepGadget::duplicate_profile(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *profile) {
        GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *copy =
                new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();

        GadgetContainerMessage<hoNDArray<std::complex<float> > > *cont_copy =
                new GadgetContainerMessage<hoNDArray<std::complex<float> > >();

        *copy->getObjectPtr() = *profile->getObjectPtr();
        *(cont_copy->getObjectPtr()) = *(AsContainerMessage<hoNDArray<std::complex<float> > >(
                profile->cont())->getObjectPtr());

        copy->cont(cont_copy);
        return copy;
    }

    std::tuple<boost::shared_ptr<hoNDArray<float_complext>>, boost::shared_ptr<hoNDArray<floatd2>>, boost::shared_ptr<hoNDArray<float>>>
    gpuSpiralSensePrepGadget::get_data_from_queues(size_t set, size_t slice, size_t num_coils) {
        unsigned int profiles_buffered = buffer_[set * slices_ + slice].size();


        auto data_host = boost::make_shared<hoNDArray<float_complext>>(
                samples_per_interleave_ * interleaves_counter_singleframe_[set * slices_ + slice] *
                ((use_multiframe_grouping_) ? acceleration_factor_ : 1),
                num_coils
        );

        std::vector<size_t> ddimensions = {(size_t) samples_per_interleave_ *
                                           interleaves_counter_singleframe_[set * slices_ + slice],
                                           (use_multiframe_grouping_) ? (size_t) acceleration_factor_ : 1};

        boost::shared_ptr<hoNDArray<floatd2> > traj_host(new hoNDArray<floatd2>(ddimensions));
        boost::shared_ptr<hoNDArray<float> > dcw_host(new hoNDArray<float>(ddimensions));

        for (unsigned int p = 0; p < profiles_buffered; p++) {

            GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *acq;
            GadgetContainerMessage<hoNDArray<std::complex<float> > > *daq;

            std::tie(acq,daq) = buffer_[set * slices_ + slice][p];

            for (unsigned int c = 0; c < num_coils; c++) {
                float_complext *data_ptr = data_host->get_data_ptr();
                data_ptr += c * samples_per_interleave_ * profiles_buffered + p * samples_per_interleave_;

                std::complex<float> *r_ptr = daq->getObjectPtr()->get_data_ptr();
                r_ptr += c * daq->getObjectPtr()->get_size(0);

                memcpy(data_ptr, r_ptr, samples_per_interleave_ * sizeof(float_complext));
            }

            floatd2 *traj_ptr = traj_host->get_data_ptr();
            traj_ptr += p * samples_per_interleave_;

            floatd2 *t_ptr = host_traj_.get_data_ptr();
            t_ptr += acq->getObjectPtr()->idx.kspace_encode_step_1 * samples_per_interleave_;

            memcpy(traj_ptr, t_ptr, samples_per_interleave_ * sizeof(floatd2));

            float *dcw_ptr = dcw_host->get_data_ptr();
            dcw_ptr += p * samples_per_interleave_;

            float *d_ptr = host_weights_.get_data_ptr();
            d_ptr += acq->getObjectPtr()->idx.kspace_encode_step_1 * samples_per_interleave_;

            memcpy(dcw_ptr, d_ptr, samples_per_interleave_ * sizeof(float));

            acq->release();
        }
        buffer_[set*slices_+slice].clear();
        return std::make_tuple(data_host, traj_host, dcw_host);
    }


    GADGET_FACTORY_DECLARE(gpuSpiralSensePrepGadget)
}
