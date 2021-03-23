
#include "GenericReconGadget.h"
#include "hoNDArray_reductions.h"
#include "mri_core_kspace_filter.h"

namespace Gadgetron {

    GenericReconGadget::GenericReconGadget() : BaseClass() {}

    GenericReconGadget::~GenericReconGadget() {}

    int GenericReconGadget::process_config(ACE_Message_Block* mb) {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try {
            deserialize(mb->rd_ptr(), h);
        } catch (...) {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!h.acquisitionSystemInformation) {
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
        calib_mode_.resize(NE, ISMRMRD_noacceleration);

        space_matrix_offset_E1_.resize(NE, 0);
        space_matrix_offset_E2_.resize(NE, 0);

        size_t e;
        for (e = 0; e < h.encoding.size(); e++) {

            if (!h.encoding[e].parallelImaging) {
                GDEBUG_STREAM("Parallel Imaging section not found in header for encoding space " << e);
                calib_mode_[e]   = ISMRMRD_noacceleration;
                acceFactorE1_[e] = 1;
                acceFactorE2_[e] = 1;
            } else {
                ISMRMRD::ParallelImaging p_imaging = *h.encoding[e].parallelImaging;

                acceFactorE1_[e] = p_imaging.accelerationFactor.kspace_encoding_step_1;
                acceFactorE2_[e] = p_imaging.accelerationFactor.kspace_encoding_step_2;
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE1 is " << acceFactorE1_[e]);
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE2 is " << acceFactorE2_[e]);

                std::string calib = *p_imaging.calibrationMode;

                bool separate    = (calib.compare("separate") == 0);
                bool embedded    = (calib.compare("embedded") == 0);
                bool external    = (calib.compare("external") == 0);
                bool interleaved = (calib.compare("interleaved") == 0);
                bool other       = (calib.compare("other") == 0);

                calib_mode_[e] = Gadgetron::ISMRMRD_noacceleration;
                if (acceFactorE1_[e] > 1 || acceFactorE2_[e] > 1) {
                    if (interleaved)
                        calib_mode_[e] = Gadgetron::ISMRMRD_interleaved;
                    else if (embedded)
                        calib_mode_[e] = Gadgetron::ISMRMRD_embedded;
                    else if (separate)
                        calib_mode_[e] = Gadgetron::ISMRMRD_separate;
                    else if (external)
                        calib_mode_[e] = Gadgetron::ISMRMRD_external;
                    else if (other)
                        calib_mode_[e] = Gadgetron::ISMRMRD_other;
                }
            }

            // -------------------------------------------------

            bool is_cartesian_sampling = (h.encoding[e].trajectory == ISMRMRD::TrajectoryType::CARTESIAN);
            bool is_epi_sampling       = (h.encoding[e].trajectory == ISMRMRD::TrajectoryType::EPI);
            if (is_cartesian_sampling || is_epi_sampling) {
                if (h.encoding[e].encodingLimits.kspace_encoding_step_1.is_present()) {
                    space_matrix_offset_E1_[e] = (int)h.encoding[e].encodedSpace.matrixSize.y / 2
                                                 - (int)h.encoding[e].encodingLimits.kspace_encoding_step_1->center;
                }

                if (h.encoding[e].encodingLimits.kspace_encoding_step_2.is_present()
                    && h.encoding[e].encodedSpace.matrixSize.z > 1) {
                    space_matrix_offset_E2_[e] = (int)h.encoding[e].encodedSpace.matrixSize.z / 2
                                                 - (int)h.encoding[e].encodingLimits.kspace_encoding_step_2->center;
                }
            }
        }

        if (!h.acquisitionSystemInformation->systemFieldStrength_T)
        {
            system_field_strength_T_ = 1.5;
        }
        else
        {
            system_field_strength_T_ = h.acquisitionSystemInformation.get().systemFieldStrength_T.get();
        }

        protocol_name_ = "Unknown";

        if (!h.measurementInformation)
        {
            GDEBUG("measurementInformation not found in header");
        }
        else
        {
            if (!h.measurementInformation->protocolName)
            {
                GDEBUG("measurementInformation->protocolName not found in header");
            }
            else
            {
                protocol_name_ = h.measurementInformation.get().protocolName.get();
            }
        }

        std::string measurement_id = "";
        if (h.measurementInformation)
        {
            if (h.measurementInformation->measurementID)
            {
                measurement_id = *h.measurementInformation->measurementID;
            }

            patient_position_ = h.measurementInformation->patientPosition;
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

        if (h.acquisitionSystemInformation)
        {
            if (h.acquisitionSystemInformation->systemVendor)
            {
                vendor_ = *h.acquisitionSystemInformation->systemVendor;
            }
        }

        return GADGET_OK;
    }

    int GenericReconGadget::process(Gadgetron::GadgetContainerMessage<IsmrmrdReconData>* m1) {
        process_called_times_++;

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_) {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : "
                         << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        // for every encoding space
        for (size_t e = 0; e < recon_bit_->rbit_.size(); e++) {
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

    void GenericReconGadget::make_ref_coil_map(IsmrmrdDataBuffered& ref_, std::vector<size_t> recon_dims,
        hoNDArray<std::complex<float>>& ref_calib, hoNDArray<std::complex<float>>& ref_coil_map, size_t encoding) {

        hoNDArray<std::complex<float>>& ref_data = ref_.data_;

        // sampling limits
        size_t sRO = ref_.sampling_.sampling_limits_[0].min_;
        size_t eRO = ref_.sampling_.sampling_limits_[0].max_;
        size_t cRO = ref_.sampling_.sampling_limits_[0].center_;

        size_t sE1 = ref_.sampling_.sampling_limits_[1].min_;
        size_t eE1 = ref_.sampling_.sampling_limits_[1].max_;
        size_t cE1 = ref_.sampling_.sampling_limits_[1].center_;

        size_t sE2 = ref_.sampling_.sampling_limits_[2].min_;
        size_t eE2 = ref_.sampling_.sampling_limits_[2].max_;
        size_t cE2 = ref_.sampling_.sampling_limits_[2].center_;

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
        if ((calib_mode_[encoding] == Gadgetron::ISMRMRD_interleaved)
            || (calib_mode_[encoding] == Gadgetron::ISMRMRD_noacceleration)) {
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
            Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(0), ref_.sampling_.sampling_limits_[0].min_,
                ref_.sampling_.sampling_limits_[0].max_, filter_RO_ref_coi_map_);

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

        std::vector<size_t> dim = *ref_data.get_dimensions();
        ref_calib.create(dim, ref_data.begin());

        if (!debug_folder_full_path_.empty()) {
            std::stringstream os;
            os << "encoding_" << encoding;

            gt_exporter_.export_array_complex(ref_coil_map, debug_folder_full_path_ + "ref_coil_map_" + os.str());
            gt_exporter_.export_array_complex(ref_calib, debug_folder_full_path_ + "ref_calib_" + os.str());
        }
    }

    void GenericReconGadget::perform_coil_map_estimation(
        const hoNDArray<std::complex<float>>& ref_coil_map, hoNDArray<std::complex<float>>& coil_map, size_t e) {
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

    void GenericReconGadget::compute_image_header(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t e) {

        size_t RO  = res.data_.get_size(0);
        size_t E1  = res.data_.get_size(1);
        size_t E2  = res.data_.get_size(2);
        size_t CHA = res.data_.get_size(3);
        size_t N   = res.data_.get_size(4);
        size_t S   = res.data_.get_size(5);
        size_t SLC = res.data_.get_size(6);

        GADGET_CHECK_THROW(N == recon_bit.data_.headers_.get_size(2));
        GADGET_CHECK_THROW(S == recon_bit.data_.headers_.get_size(3));

        res.headers_.create(N, S, SLC);
        res.meta_.resize(N * S * SLC);

        size_t n, s, slc;

        for (slc = 0; slc < SLC; slc++) {
            for (s = 0; s < S; s++) {
                for (n = 0; n < N; n++) {
                    size_t header_E1 = recon_bit.data_.headers_.get_size(0);
                    size_t header_E2 = recon_bit.data_.headers_.get_size(1);

                    // for every kspace, find the recorded header which is closest to the kspace center [E1/2 E2/2]
                    ISMRMRD::AcquisitionHeader acq_header;

                    long long bestE1 = E1 + 1;
                    long long bestE2 = E2 + 1;

                    size_t e1, e2;
                    for (e2 = 0; e2 < header_E2; e2++) {
                        for (e1 = 0; e1 < header_E1; e1++) {
                            ISMRMRD::AcquisitionHeader& curr_header = recon_bit.data_.headers_(e1, e2, n, s, slc);

                            long long e1_in_bucket = curr_header.idx.kspace_encode_step_1 + space_matrix_offset_E1_[e];

                            if (E2 > 1) {
                                long long e2_in_bucket
                                    = curr_header.idx.kspace_encode_step_2 + space_matrix_offset_E2_[e];

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

                    ISMRMRD::ImageHeader& im_header = res.headers_(n, s, slc);
                    ISMRMRD::MetaContainer& meta    = res.meta_[n + s * N + slc * N * S];

                    im_header.version         = acq_header.version;
                    im_header.data_type       = ISMRMRD::ISMRMRD_CXFLOAT;
                    im_header.flags           = acq_header.flags;
                    im_header.measurement_uid = acq_header.measurement_uid;

                    im_header.matrix_size[0] = (uint16_t)RO;
                    im_header.matrix_size[1] = (uint16_t)E1;
                    im_header.matrix_size[2] = (uint16_t)E2;

                    im_header.field_of_view[0] = recon_bit.data_.sampling_.recon_FOV_[0];
                    im_header.field_of_view[1] = recon_bit.data_.sampling_.recon_FOV_[1];
                    im_header.field_of_view[2] = recon_bit.data_.sampling_.recon_FOV_[2];

                    im_header.channels = (uint16_t)CHA;

                    im_header.position[0] = acq_header.position[0];
                    im_header.position[1] = acq_header.position[1];
                    im_header.position[2] = acq_header.position[2];

                    im_header.read_dir[0] = acq_header.read_dir[0];
                    im_header.read_dir[1] = acq_header.read_dir[1];
                    im_header.read_dir[2] = acq_header.read_dir[2];

                    im_header.phase_dir[0] = acq_header.phase_dir[0];
                    im_header.phase_dir[1] = acq_header.phase_dir[1];
                    im_header.phase_dir[2] = acq_header.phase_dir[2];

                    im_header.slice_dir[0] = acq_header.slice_dir[0];
                    im_header.slice_dir[1] = acq_header.slice_dir[1];
                    im_header.slice_dir[2] = acq_header.slice_dir[2];

                    im_header.patient_table_position[0] = acq_header.patient_table_position[0];
                    im_header.patient_table_position[1] = acq_header.patient_table_position[1];
                    im_header.patient_table_position[2] = acq_header.patient_table_position[2];

                    im_header.average    = acq_header.idx.average;
                    im_header.slice      = acq_header.idx.slice;
                    im_header.contrast   = acq_header.idx.contrast;
                    im_header.phase      = acq_header.idx.phase;
                    im_header.repetition = acq_header.idx.repetition;
                    im_header.set        = acq_header.idx.set;

                    im_header.acquisition_time_stamp = acq_header.acquisition_time_stamp;

                    im_header.physiology_time_stamp[0] = acq_header.physiology_time_stamp[0];
                    im_header.physiology_time_stamp[1] = acq_header.physiology_time_stamp[1];
                    im_header.physiology_time_stamp[2] = acq_header.physiology_time_stamp[2];

                    im_header.image_type         = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
                    im_header.image_index        = (uint16_t)(n + s * N + slc * N * S);
                    im_header.image_series_index = 0;

                    memcpy(im_header.user_int, acq_header.user_int, sizeof(int32_t) * ISMRMRD::ISMRMRD_USER_INTS);
                    memcpy(im_header.user_float, acq_header.user_float, sizeof(float) * ISMRMRD::ISMRMRD_USER_FLOATS);

                    im_header.attribute_string_len = 0;

                    meta.set("encoding", (long)e);

                    meta.set("encoding_FOV", recon_bit.data_.sampling_.encoded_FOV_[0]);
                    meta.append("encoding_FOV", recon_bit.data_.sampling_.encoded_FOV_[1]);
                    meta.append("encoding_FOV", recon_bit.data_.sampling_.encoded_FOV_[2]);

                    meta.set("recon_FOV", recon_bit.data_.sampling_.recon_FOV_[0]);
                    meta.append("recon_FOV", recon_bit.data_.sampling_.recon_FOV_[1]);
                    meta.append("recon_FOV", recon_bit.data_.sampling_.recon_FOV_[2]);

                    meta.set("encoded_matrix", (long)recon_bit.data_.sampling_.encoded_matrix_[0]);
                    meta.append("encoded_matrix", (long)recon_bit.data_.sampling_.encoded_matrix_[1]);
                    meta.append("encoded_matrix", (long)recon_bit.data_.sampling_.encoded_matrix_[2]);

                    meta.set("recon_matrix", (long)recon_bit.data_.sampling_.recon_matrix_[0]);
                    meta.append("recon_matrix", (long)recon_bit.data_.sampling_.recon_matrix_[1]);
                    meta.append("recon_matrix", (long)recon_bit.data_.sampling_.recon_matrix_[2]);

                    meta.set("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].min_);
                    meta.append("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].center_);
                    meta.append("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].max_);

                    meta.set("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].min_);
                    meta.append("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].center_);
                    meta.append("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].max_);

                    meta.set("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].min_);
                    meta.append("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].center_);
                    meta.append("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].max_);

                    meta.set("PatientPosition", (double)res.headers_(n, s, slc).position[0]);
                    meta.append("PatientPosition", (double)res.headers_(n, s, slc).position[1]);
                    meta.append("PatientPosition", (double)res.headers_(n, s, slc).position[2]);

                    meta.set("read_dir", (double)res.headers_(n, s, slc).read_dir[0]);
                    meta.append("read_dir", (double)res.headers_(n, s, slc).read_dir[1]);
                    meta.append("read_dir", (double)res.headers_(n, s, slc).read_dir[2]);

                    meta.set("phase_dir", (double)res.headers_(n, s, slc).phase_dir[0]);
                    meta.append("phase_dir", (double)res.headers_(n, s, slc).phase_dir[1]);
                    meta.append("phase_dir", (double)res.headers_(n, s, slc).phase_dir[2]);

                    meta.set("slice_dir", (double)res.headers_(n, s, slc).slice_dir[0]);
                    meta.append("slice_dir", (double)res.headers_(n, s, slc).slice_dir[1]);
                    meta.append("slice_dir", (double)res.headers_(n, s, slc).slice_dir[2]);

                    meta.set("patient_table_position", (double)res.headers_(n, s, slc).patient_table_position[0]);
                    meta.append("patient_table_position", (double)res.headers_(n, s, slc).patient_table_position[1]);
                    meta.append("patient_table_position", (double)res.headers_(n, s, slc).patient_table_position[2]);

                    meta.set("acquisition_time_stamp", (long)res.headers_(n, s, slc).acquisition_time_stamp);

                    meta.set("physiology_time_stamp", (long)res.headers_(n, s, slc).physiology_time_stamp[0]);
                    meta.append("physiology_time_stamp", (long)res.headers_(n, s, slc).physiology_time_stamp[1]);
                    meta.append("physiology_time_stamp", (long)res.headers_(n, s, slc).physiology_time_stamp[2]);

                    size_t ui;
                    for (ui = 0; ui < ISMRMRD::ISMRMRD_USER_INTS; ui++)
                    {
                        std::ostringstream str;
                        str << "user_int_" << ui;
                        meta.append(str.str().c_str(), (long)res.headers_(n, s, slc).user_int[ui]);
                    }

                    for (ui = 0; ui < ISMRMRD::ISMRMRD_USER_FLOATS; ui++)
                    {
                        std::ostringstream str;
                        str << "user_float_" << ui;
                        meta.append(str.str().c_str(), (long)res.headers_(n, s, slc).user_float[ui]);
                    }

                    meta.set("gadgetron_sha1", GADGETRON_SHA1);

                    meta.set("measurementID", this->measurement_id_.c_str());
                    meta.set("protocolName", this->protocol_name_.c_str());
                    meta.set("patientID", this->patient_.c_str());
                    meta.set("studyID", this->study_.c_str());
                    meta.set("measurementNumber", this->measurement_.c_str());
                    meta.set("deviceID", this->device_.c_str());
                    meta.set("patient_position", this->patient_position_.c_str());
                }
            }
        }
    }

    void GenericReconGadget::compute_snr_scaling_factor(
        IsmrmrdReconBit& recon_bit, float& effective_acce_factor, float& snr_scaling_ratio) {

        size_t RO     = recon_bit.data_.data_.get_size(0);
        size_t E1     = recon_bit.data_.data_.get_size(1);
        size_t E2     = recon_bit.data_.data_.get_size(2);
        size_t dstCHA = recon_bit.data_.data_.get_size(3);
        size_t N      = recon_bit.data_.data_.get_size(4);
        size_t S      = recon_bit.data_.data_.get_size(5);
        size_t SLC    = recon_bit.data_.data_.get_size(6);

        effective_acce_factor = 1;
        snr_scaling_ratio     = 1;

        size_t e1, e2, n, s;
        size_t num_readout_lines = 0;
        for (s = 0; s < S; s++) {
            for (n = 0; n < N; n++) {
                for (e2 = 0; e2 < E2; e2++) {
                    for (e1 = 0; e1 < E1; e1++) {
                        if (std::abs(recon_bit.data_.data_(RO / 2, e1, e2, 0, n, 0, 0)) > 0) {
                            num_readout_lines++;
                        }
                    }
                }
            }
        }

        if (num_readout_lines > 0) {
            float lenRO = RO;

            size_t start_RO = recon_bit.data_.sampling_.sampling_limits_[0].min_;
            size_t end_RO   = recon_bit.data_.sampling_.sampling_limits_[0].max_;

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
    void GenericReconGadget::send_out_image_array(
        IsmrmrdImageArray& res, size_t encoding, int series_num, const std::string& data_role) {
        this->prepare_image_array(res, encoding, series_num, data_role);
        this->next()->putq(new GadgetContainerMessage<IsmrmrdImageArray>(res));
    }

    GADGET_FACTORY_DECLARE(GenericReconGadget)
}
