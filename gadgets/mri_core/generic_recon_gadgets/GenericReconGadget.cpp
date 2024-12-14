
#include "GenericReconGadget.h"
#include "hoNDArray_reductions.h"
#include "mri_core_kspace_filter.h"

namespace Gadgetron {

    GenericReconGadget::GenericReconGadget() : BaseClass() {}

    GenericReconGadget::~GenericReconGadget() {}

    int GenericReconGadget::process_config(const mrd::Header& header) {
        GADGET_CHECK_RETURN(BaseClass::process_config(header) == GADGET_OK, GADGET_FAIL);

        auto& h = header;

        if (!h.acquisition_system_information) {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        this->initialize_encoding_space_limits(h);
        // -------------------------------------------------

        size_t NE            = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        acceFactorE1_.resize(NE, 1);
        acceFactorE2_.resize(NE, 1);
        calib_mode_.resize(NE, mrd::CalibrationMode::kNoacceleration);

        space_matrix_offset_E1_.resize(NE, 0);
        space_matrix_offset_E2_.resize(NE, 0);

        size_t e;
        for (e = 0; e < h.encoding.size(); e++) {

            if (!h.encoding[e].parallel_imaging) {
                GDEBUG_STREAM("Parallel Imaging section not found in header for encoding space " << e);
                calib_mode_[e]   = mrd::CalibrationMode::kNoacceleration;
                acceFactorE1_[e] = 1;
                acceFactorE2_[e] = 1;
            } else {
                mrd::ParallelImagingType p_imaging = *h.encoding[e].parallel_imaging;

                acceFactorE1_[e] = p_imaging.acceleration_factor.kspace_encoding_step_1;
                acceFactorE2_[e] = p_imaging.acceleration_factor.kspace_encoding_step_2;
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE1 is " << acceFactorE1_[e]);
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE2 is " << acceFactorE2_[e]);

                calib_mode_[e] = mrd::CalibrationMode::kNoacceleration;
                if (acceFactorE1_[e] > 1 || acceFactorE2_[e] > 1) {
                    calib_mode_[e] = p_imaging.calibration_mode.value_or(mrd::CalibrationMode::kNoacceleration);
                }
            }

            // -------------------------------------------------

            bool is_cartesian_sampling = (h.encoding[e].trajectory == mrd::Trajectory::kCartesian);
            bool is_epi_sampling       = (h.encoding[e].trajectory == mrd::Trajectory::kEpi);
            if (is_cartesian_sampling || is_epi_sampling) {
                if (h.encoding[e].encoding_limits.kspace_encoding_step_1.has_value()) {
                    space_matrix_offset_E1_[e] = (int)h.encoding[e].encoded_space.matrix_size.y / 2
                                                 - (int)h.encoding[e].encoding_limits.kspace_encoding_step_1->center;
                }

                if (h.encoding[e].encoding_limits.kspace_encoding_step_2.has_value()
                    && h.encoding[e].encoded_space.matrix_size.z > 1) {
                    space_matrix_offset_E2_[e] = (int)h.encoding[e].encoded_space.matrix_size.z / 2
                                                 - (int)h.encoding[e].encoding_limits.kspace_encoding_step_2->center;
                }
            }
        }

        system_field_strength_T_ = h.acquisition_system_information->system_field_strength_t.value_or(1.5);

        protocol_name_ = "Unknown";

        if (!h.measurement_information)
        {
            GDEBUG("measurementInformation not found in header");
        }
        else
        {
            if (!h.measurement_information->protocol_name)
            {
                GDEBUG("measurementInformation->protocolName not found in header");
            }
            else
            {
                protocol_name_ = h.measurement_information->protocol_name.value();
            }
        }

        std::string measurement_id = "";
        if (h.measurement_information)
        {
            if (h.measurement_information->measurement_id)
            {
                measurement_id = *h.measurement_information->measurement_id;
            }

            patient_position_ = h.measurement_information->patient_position;
        }

        this->measurement_id_ = measurement_id;

        // analyze measurement id
        if (measurement_id.size() > 0)
        {
            std::string mid = measurement_id;
            size_t ind = mid.find("_");
            if (ind != std::string::npos)
            {
                device_ = mid.substr(0, ind);
                mid = mid.substr(ind + 1, std::string::npos);

                ind = mid.find("_");
                if (ind != std::string::npos)
                {
                    patient_ = mid.substr(0, ind);
                    mid = mid.substr(ind + 1, std::string::npos);

                    ind = mid.find("_");
                    if (ind != std::string::npos)
                    {
                        study_ = mid.substr(0, ind);
                        measurement_ = mid.substr(ind + 1, std::string::npos);
                    }
                }
            }
        }

        if (h.acquisition_system_information->system_vendor)
        {
            vendor_ = h.acquisition_system_information->system_vendor.value();
        }

        return GADGET_OK;
    }

    int GenericReconGadget::process(Gadgetron::GadgetContainerMessage<mrd::ReconData>* m1) {
        process_called_times_++;

        mrd::ReconData* recon_data = m1->getObjectPtr();
        if (recon_data->buffers.size() > num_encoding_spaces_) {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : "
                         << recon_data->buffers.size() << " instead of " << num_encoding_spaces_);
        }

        // for every encoding space
        for (size_t e = 0; e < recon_data->buffers.size(); e++) {
            std::stringstream os;
            os << "_encoding_" << e;

            GDEBUG_CONDITION_STREAM(
                verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(
                verbose.value(), "======================================================================");

            // ---------------------------------------------------------------
            // export incoming data

            /**
                        if (!debug_folder_full_path_.empty())
                        {
                            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0)
                            {
                                gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_,
            debug_folder_full_path_ + "data" + os.str());
                            }
                        }

                        if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].data_.trajectory_)
                        {
                            if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0)
                            {
                                gt_exporter_.export_array(*(recon_bit_->rbit_[e].data_.trajectory_),
            debug_folder_full_path_ + "data_traj" + os.str());
                            }
                        }
            **/

            // ---------------------------------------------------------------
            // export incoming ref

            /**
                        if (recon_bit_->rbit_[e].ref_)
                        {
                            if (!debug_folder_full_path_.empty())
                            {
                                gt_exporter_.export_array_complex(recon_bit_->rbit_[e].ref_->data_,
            debug_folder_full_path_ + "ref" + os.str());
                            }

                            if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].ref_->trajectory_)
                            {
                                if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0)
                                {
                                    gt_exporter_.export_array(*(recon_bit_->rbit_[e].ref_->trajectory_),
            debug_folder_full_path_ + "ref_traj" + os.str());
                                }
                            }
                        }
            **/

            // ---------------------------------------------------------------
            // add recon code here ...
        }

        m1->release();
        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    void GenericReconGadget::make_ref_coil_map(mrd::ReconBuffer& ref, std::vector<size_t> recon_dims,
        hoNDArray<std::complex<float>>& ref_calib, hoNDArray<std::complex<float>>& ref_coil_map, size_t encoding)
    {
        hoNDArray<std::complex<float>>& ref_data = ref.data;

        // sampling limits
        size_t sRO = ref.sampling.sampling_limits.kspace_encoding_step_0.minimum;
        size_t eRO = ref.sampling.sampling_limits.kspace_encoding_step_0.maximum;
        size_t cRO = ref.sampling.sampling_limits.kspace_encoding_step_0.center;

        size_t sE1 = ref.sampling.sampling_limits.kspace_encoding_step_1.minimum;
        size_t eE1 = ref.sampling.sampling_limits.kspace_encoding_step_1.maximum;
        size_t cE1 = ref.sampling.sampling_limits.kspace_encoding_step_1.center;

        size_t sE2 = ref.sampling.sampling_limits.kspace_encoding_step_2.minimum;
        size_t eE2 = ref.sampling.sampling_limits.kspace_encoding_step_2.maximum;
        size_t cE2 = ref.sampling.sampling_limits.kspace_encoding_step_2.center;

        // recon size
        size_t recon_RO = recon_dims[0];
        size_t recon_E1 = recon_dims[1];
        size_t recon_E2 = recon_dims[2];

        // ref array size
        size_t ref_data_RO = ref_data.get_size(0);
        size_t ref_data_E1 = ref_data.get_size(1);
        size_t ref_data_E2 = ref_data.get_size(2);
        size_t CHA         = ref_data.get_size(3);
        size_t N           = ref_data.get_size(4);
        size_t S           = ref_data.get_size(5);
        size_t SLC         = ref_data.get_size(6);

        // determine the ref_coil_map size
        size_t RO = 2 * cRO;
        if (sRO > 0 || eRO < RO - 1) {
            RO = 2 * std::max(cRO - sRO, eRO - cRO + 1);
            if (RO > recon_RO)
                RO = recon_RO;
        }

        size_t E1 = eE1 - sE1 + 1;
        size_t E2 = eE2 - sE2 + 1;

        // cut the center region for ref coil map
        if ((calib_mode_[encoding] == mrd::CalibrationMode::kInterleaved)
            || (calib_mode_[encoding] == mrd::CalibrationMode::kNoacceleration)) {
            E1 = 2 * std::min(cE1 - sE1, eE1 - cE1 + 1) - 1;
            if (E1 > recon_E1)
                E1 = recon_E1;

            if (E2 > 1) {
                E2 = 2 * std::min(cE2 - sE2, eE2 - cE2 + 1) - 1;
                if (E2 > recon_E2)
                    E2 = recon_E2;
            }
        }

        ref_coil_map.create(RO, E1, E2, CHA, N, S, SLC);
        Gadgetron::clear(ref_coil_map);

        // make sure center aligned
        size_t ref_sE1 = cE1 - E1 / 2;
        size_t ref_eE1 = ref_sE1 + E1 - 1;

        size_t ref_sE2 = cE2 - E2 / 2;
        size_t ref_eE2 = ref_sE2 + E2 - 1;

        size_t slc, s, n, cha, e2, e1;
        for (slc = 0; slc < SLC; slc++) {
            for (s = 0; s < S; s++) {
                for (n = 0; n < N; n++) {
                    for (cha = 0; cha < CHA; cha++) {
                        for (e2 = ref_sE2; e2 <= ref_eE2; e2++) {
                            for (e1 = ref_sE1; e1 <= ref_eE1; e1++) {
                                std::complex<float>* pSrc = &(ref_data(0, e1, e2, cha, n, s, slc));
                                std::complex<float>* pDst
                                    = &(ref_coil_map(0, e1 - ref_sE1, e2 - ref_sE2, cha, n, s, slc));

                                memcpy(pDst + sRO, pSrc, sizeof(std::complex<float>) * (eRO - sRO + 1));
                            }
                        }
                    }
                }
            }
        }

        if (!debug_folder_full_path_.empty()) {
            std::stringstream os;
            os << "encoding_" << encoding;

            gt_exporter_.export_array_complex(
                ref_coil_map, debug_folder_full_path_ + "ref_coil_map_before_filtering_" + os.str());
        }

        // filter the ref_coil_map
        if (filter_RO_ref_coi_map_.get_size(0) != RO) {
            Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(0), ref.sampling.sampling_limits.kspace_encoding_step_0.minimum,
                ref.sampling.sampling_limits.kspace_encoding_step_0.maximum, filter_RO_ref_coi_map_);

            if (!debug_folder_full_path_.empty()) {
                std::stringstream os;
                os << "encoding_" << encoding;

                gt_exporter_.export_array_complex(
                    filter_RO_ref_coi_map_, debug_folder_full_path_ + "filter_RO_ref_coi_map_" + os.str());
            }
        }

        if (filter_E1_ref_coi_map_.get_size(0) != E1) {
            Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(1), 0, E1 - 1, filter_E1_ref_coi_map_);

            if (!debug_folder_full_path_.empty()) {
                std::stringstream os;
                os << "encoding_" << encoding;

                gt_exporter_.export_array_complex(
                    filter_E1_ref_coi_map_, debug_folder_full_path_ + "filter_E1_ref_coi_map_" + os.str());
            }
        }

        if ((E2 > 1) && (filter_E2_ref_coi_map_.get_size(0) != E2)) {
            Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(2), 0, E2 - 1, filter_E2_ref_coi_map_);

            if (!debug_folder_full_path_.empty()) {
                std::stringstream os;
                os << "encoding_" << encoding;

                gt_exporter_.export_array_complex(
                    filter_E2_ref_coi_map_, debug_folder_full_path_ + "filter_E2_ref_coi_map_" + os.str());
            }
        }

        hoNDArray<std::complex<float>> ref_recon_buf;

        if (E2 > 1) {
            Gadgetron::apply_kspace_filter_ROE1E2(
                ref_coil_map, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_, ref_recon_buf);
        } else {
            Gadgetron::apply_kspace_filter_ROE1(
                ref_coil_map, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, ref_recon_buf);
        }

        if (!debug_folder_full_path_.empty()) {
            std::stringstream os;
            os << "encoding_" << encoding;

            gt_exporter_.export_array_complex(
                ref_recon_buf, debug_folder_full_path_ + "ref_coil_map_after_filtering_" + os.str());
        }

        // pad the ref_coil_map into the data array
        Gadgetron::pad(recon_RO, recon_E1, recon_E2, ref_recon_buf, ref_coil_map);

        std::vector<size_t> dim = ref_data.get_dimensions();
        ref_calib.create(dim, ref_data.begin());

        if (!debug_folder_full_path_.empty()) {
            std::stringstream os;
            os << "encoding_" << encoding;

            gt_exporter_.export_array_complex(ref_coil_map, debug_folder_full_path_ + "ref_coil_map_" + os.str());
            gt_exporter_.export_array_complex(ref_calib, debug_folder_full_path_ + "ref_calib_" + os.str());
        }
    }

    void GenericReconGadget::perform_coil_map_estimation(
        const hoNDArray<std::complex<float>>& ref_coil_map, hoNDArray<std::complex<float>>& coil_map, size_t e)
    {
        try {
            coil_map = ref_coil_map;
            Gadgetron::clear(coil_map);

            size_t E2 = ref_coil_map.get_size(2);
            if (E2 > 1) {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(ref_coil_map, complex_im_recon_buf_);
            } else {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(ref_coil_map, complex_im_recon_buf_);
            }

            if (!debug_folder_full_path_.empty()) {
                std::stringstream os;
                os << "encoding_" << e;

                gt_exporter_.export_array_complex(
                    complex_im_recon_buf_, debug_folder_full_path_ + "complex_im_for_coil_map_" + os.str());
            }

            if (coil_map_algorithm.value() == "Inati") {
                size_t ks    = this->coil_map_kernel_size_readout.value();
                size_t kz    = this->coil_map_kernel_size_phase.value();
                size_t power = 3;

                Gadgetron::coil_map_Inati(complex_im_recon_buf_, coil_map, ks, kz, power);
            } else {
                size_t ks      = this->coil_map_kernel_size_readout.value();
                size_t kz      = this->coil_map_kernel_size_phase.value();
                size_t iterNum = this->coil_map_num_iter.value();
                float thres    = this->coil_map_thres_iter.value();

                Gadgetron::coil_map_Inati_Iter(complex_im_recon_buf_, coil_map, ks, kz, iterNum, thres);
            }

            if (!debug_folder_full_path_.empty()) {
                std::stringstream os;
                os << "encoding_" << e;

                gt_exporter_.export_array_complex(coil_map, debug_folder_full_path_ + "coil_map_" + os.str());
            }
        } catch (...) {
            GADGET_THROW("Errors happened in GenericReconGadget::perform_coil_map_estimation(...) ... ");
        }
    }

    void GenericReconGadget::compute_image_header(mrd::ReconAssembly& recon_bit, mrd::ImageArray& res, size_t e)
    {
        size_t RO  = res.data.get_size(0);
        size_t E1  = res.data.get_size(1);
        size_t E2  = res.data.get_size(2);
        size_t CHA = res.data.get_size(3);
        size_t N   = res.data.get_size(4);
        size_t S   = res.data.get_size(5);
        size_t SLC = res.data.get_size(6);

        GADGET_CHECK_THROW(N == recon_bit.data.headers.get_size(2));
        GADGET_CHECK_THROW(S == recon_bit.data.headers.get_size(3));
        GADGET_CHECK_THROW(SLC == recon_bit.data.headers.get_size(4));

        res.headers.create(N, S, SLC);
        res.meta.create(N, S, SLC);

        size_t n, s, slc;

        for (slc = 0; slc < SLC; slc++) {
            for (s = 0; s < S; s++) {
                for (n = 0; n < N; n++) {
                    size_t header_E1 = recon_bit.data.headers.get_size(0);
                    size_t header_E2 = recon_bit.data.headers.get_size(1);

                    // for every kspace, find the recorded header which is closest to the kspace center [E1/2 E2/2]
                    mrd::AcquisitionHeader acq_header;

                    // for every kspace, find the min and max of acquisition time, find the min and max of physio time
                    uint32_t min_acq_time(std::numeric_limits<uint32_t>::max()),
                        max_acq_time(0);

                    // Only the first three physio timestamps are saved in the ImageMeta
                    const size_t MAX_PHYSIO_TIMESTAMPS = 3;
                    std::vector<size_t> min_physio_time(MAX_PHYSIO_TIMESTAMPS, std::numeric_limits<uint32_t>::max());
                    std::vector<size_t> max_physio_time(MAX_PHYSIO_TIMESTAMPS, 0);

                    long long bestE1 = E1 + 1;
                    long long bestE2 = E2 + 1;

                    for (size_t e2 = 0; e2 < header_E2; e2++) {
                        for (size_t e1 = 0; e1 < header_E1; e1++) {
                            mrd::AcquisitionHeader& curr_header = recon_bit.data.headers(e1, e2, n, s, slc);

                            if (curr_header.acquisition_time_stamp>0) {
                                if (min_acq_time > curr_header.acquisition_time_stamp)
                                    min_acq_time = curr_header.acquisition_time_stamp.value_or(0);

                                if (max_acq_time < curr_header.acquisition_time_stamp)
                                    max_acq_time = curr_header.acquisition_time_stamp.value_or(0);

                                for (size_t ii = 0; ii < MAX_PHYSIO_TIMESTAMPS; ii++) {
                                    if (ii >= curr_header.physiology_time_stamp.size()) {
                                        break;
                                    }

                                    if (min_physio_time[ii] > curr_header.physiology_time_stamp[ii]) {
                                        min_physio_time[ii] = curr_header.physiology_time_stamp[ii];
                                    }

                                    if (max_physio_time[ii] < curr_header.physiology_time_stamp[ii]) {
                                        max_physio_time[ii] = curr_header.physiology_time_stamp[ii];
                                    }
                                }
                            }

                            long long e1_in_bucket = curr_header.idx.kspace_encode_step_1.value_or(0) + space_matrix_offset_E1_[e];

                            if (E2 > 1) {
                                long long e2_in_bucket
                                    = curr_header.idx.kspace_encode_step_2.value_or(0) + space_matrix_offset_E2_[e];

                                if (std::abs(e1_in_bucket - (long long)(E1 / 2)) < bestE1
                                    && std::abs(e2_in_bucket - (long long)(E2 / 2)) < bestE2) {
                                    bestE1 = std::abs(e1_in_bucket - (long long)E1 / 2);
                                    bestE2 = std::abs(e2_in_bucket - (long long)E2 / 2);

                                    acq_header = curr_header;
                                }
                            } else {
                                if (std::abs(e1_in_bucket - (long long)(E1 / 2)) < bestE1) {
                                    bestE1 = std::abs(e1_in_bucket - (long long)E1 / 2);

                                    acq_header = curr_header;
                                }
                            }
                        }
                    }

                    mrd::ImageHeader& im_header = res.headers(n, s, slc);
                    mrd::ImageMeta& meta = res.meta(n, s, slc);
                    mrd::SamplingDescription& sampling = recon_bit.data.sampling;

                    im_header.measurement_uid = acq_header.measurement_uid;

                    im_header.field_of_view[0] = sampling.recon_fov.x;
                    im_header.field_of_view[1] = sampling.recon_fov.y;
                    im_header.field_of_view[2] = sampling.recon_fov.z;

                    im_header.position = acq_header.position;

                    im_header.col_dir = acq_header.read_dir;

                    im_header.line_dir = acq_header.phase_dir;

                    im_header.slice_dir = acq_header.slice_dir;

                    im_header.patient_table_position = acq_header.patient_table_position;

                    im_header.average    = acq_header.idx.average;
                    im_header.slice      = acq_header.idx.slice;
                    im_header.contrast   = acq_header.idx.contrast;
                    im_header.phase      = acq_header.idx.phase;
                    im_header.repetition = acq_header.idx.repetition;
                    im_header.set        = acq_header.idx.set;

                    im_header.acquisition_time_stamp = acq_header.acquisition_time_stamp;

                    im_header.physiology_time_stamp = acq_header.physiology_time_stamp;

                    im_header.image_type         = mrd::ImageType::kMagnitude;
                    im_header.image_index        = n + s * N + slc * N * S;
                    im_header.image_series_index = 0;

                    im_header.user_int = acq_header.user_int;
                    im_header.user_float = acq_header.user_float;

                    meta["encoding"] = {(long)e};

                    meta["encoding_FOV"] = {sampling.encoded_fov.x, sampling.encoded_fov.y, sampling.encoded_fov.z};

                    meta["recon_FOV"] = {sampling.recon_fov.x, sampling.recon_fov.y, sampling.recon_fov.z};

                    meta["encoded_matrix"] = {(long)sampling.encoded_matrix.x, (long)sampling.encoded_matrix.y, (long)sampling.encoded_matrix.z};

                    meta["recon_matrix"] = {(long)sampling.recon_matrix.x, (long)sampling.recon_matrix.y, (long)sampling.recon_matrix.z};

                    meta["sampling_limits_RO"] = {(long)sampling.sampling_limits.kspace_encoding_step_0.minimum, (long)sampling.sampling_limits.kspace_encoding_step_0.center, (long)sampling.sampling_limits.kspace_encoding_step_0.maximum};

                    meta["sampling_limits_E1"] = {(long)sampling.sampling_limits.kspace_encoding_step_1.minimum, (long)sampling.sampling_limits.kspace_encoding_step_1.center, (long)sampling.sampling_limits.kspace_encoding_step_1.maximum};

                    meta["sampling_limits_E2"] = {(long)sampling.sampling_limits.kspace_encoding_step_2.minimum, (long)sampling.sampling_limits.kspace_encoding_step_2.center, (long)sampling.sampling_limits.kspace_encoding_step_2.maximum};

                    meta["PatientPosition"] = {im_header.position[0], im_header.position[1], im_header.position[2]};

                    meta["read_dir"] = {im_header.col_dir[0], im_header.col_dir[1], im_header.col_dir[2]};

                    meta["phase_dir"] = {im_header.line_dir[0], im_header.line_dir[1], im_header.line_dir[2]};

                    meta["slice_dir"] = {im_header.slice_dir[0], im_header.slice_dir[1], im_header.slice_dir[2]};

                    meta["patient_table_position"] = {im_header.patient_table_position[0], im_header.patient_table_position[1], im_header.patient_table_position[2]};

                    meta["acquisition_time_stamp"] = {(long)im_header.acquisition_time_stamp.value_or(0)};

                    std::transform(im_header.physiology_time_stamp.begin(), im_header.physiology_time_stamp.end(), std::back_inserter(meta["physiology_time_stamp"]), [](const auto& i) { return (long)i; });

                    meta["acquisition_time_range"] = {(long)min_acq_time, (long)max_acq_time};

                    meta["physiology_time_range"] = {(long)min_physio_time[0], (long)max_physio_time[0], (long)min_physio_time[1], (long)max_physio_time[1], (long)min_physio_time[2], (long)max_physio_time[2]};

                    for (size_t i = 0; i < im_header.user_int.size(); i++) {
                        std::stringstream str;
                        str << "user_int_" << i;
                        meta[str.str()] = {(long)im_header.user_int[i]};
                    }

                    for (size_t i = 0; i < im_header.user_float.size(); i++) {
                        std::stringstream str;
                        str << "user_float_" << i;
                        meta[str.str()] = {im_header.user_float[i]};
                    }

                    meta["gadgetron_sha1"] = {PINGVIN_SHA1};

                    meta["measurementID"] = {this->measurement_id_};
                    meta["protocolName"] = {this->protocol_name_};
                    meta["patientID"] = {this->patient_};
                    meta["studyID"] = {this->study_};
                    meta["measurementNumber"] = {this->measurement_};
                    meta["deviceID"] = {this->device_};
                    meta["patient_position"] = {(long)this->patient_position_};
                }
            }
        }
    }

    void GenericReconGadget::compute_snr_scaling_factor(mrd::ReconAssembly& recon_bit, float& effective_acce_factor, float& snr_scaling_ratio)
    {
        size_t RO     = recon_bit.data.data.get_size(0);
        size_t E1     = recon_bit.data.data.get_size(1);
        size_t E2     = recon_bit.data.data.get_size(2);
        size_t dstCHA = recon_bit.data.data.get_size(3);
        size_t N      = recon_bit.data.data.get_size(4);
        size_t S      = recon_bit.data.data.get_size(5);
        size_t SLC    = recon_bit.data.data.get_size(6);

        effective_acce_factor = 1;
        snr_scaling_ratio     = 1;

        size_t e1, e2, n, s;
        size_t num_readout_lines = 0;
        for (s = 0; s < S; s++) {
            for (n = 0; n < N; n++) {
                for (e2 = 0; e2 < E2; e2++) {
                    for (e1 = 0; e1 < E1; e1++) {
                        if (std::abs(recon_bit.data.data(RO / 2, e1, e2, 0, n, 0, 0)) > 0) {
                            num_readout_lines++;
                        }
                    }
                }
            }
        }

        if (num_readout_lines > 0) {
            float lenRO = RO;

            size_t start_RO = recon_bit.data.sampling.sampling_limits.kspace_encoding_step_0.minimum;
            size_t end_RO   = recon_bit.data.sampling.sampling_limits.kspace_encoding_step_0.maximum;

            if ((start_RO < RO) && (end_RO < RO) && (end_RO - start_RO + 1 < RO)) {
                lenRO = (end_RO - start_RO + 1);
            }
            if (this->verbose.value())
                GDEBUG_STREAM("length for RO : " << lenRO << " - " << lenRO / RO);

            effective_acce_factor = (float)(S * N * E1 * E2) / (num_readout_lines);
            if (this->verbose.value())
                GDEBUG_STREAM("effective_acce_factor : " << effective_acce_factor);

            float ROScalingFactor = (float)RO / (float)lenRO;

            snr_scaling_ratio = (float)(std::sqrt(ROScalingFactor * effective_acce_factor));

            if (this->verbose.value())
                GDEBUG_STREAM("snr_scaling_ratio : " << snr_scaling_ratio);
        } else {
            GWARN_STREAM("Cannot find any sampled lines ... ");
        }
    }

    void GenericReconGadget::send_out_image_array(mrd::ImageArray& res, size_t encoding, int series_num, const std::string& data_role)
    {
        this->prepare_image_array(res, encoding, series_num, data_role);
        this->next()->putq(new GadgetContainerMessage<mrd::ImageArray>(res));
    }

    GADGET_FACTORY_DECLARE(GenericReconGadget)
}
