//
// Created by dchansen on 10/23/18.
//

#include "NonCartesianTools.h"

void Gadgetron::NonCartesian::append_image_header( IsmrmrdImageArray &res, const IsmrmrdReconBit &recon_bit,size_t encoding) {

    size_t RO = res.data_.get_size(0);
    size_t E1 = res.data_.get_size(1);
    size_t E2 = res.data_.get_size(2);
    size_t CHA = res.data_.get_size(3);
    size_t N = res.data_.get_size(4);
    size_t S = res.data_.get_size(5);
    size_t SLC = res.data_.get_size(6);


    res.headers_.create(N, S, SLC);
    res.meta_.resize(N * S * SLC);

    size_t n, s, slc;

    for (slc = 0; slc < SLC; slc++) {
        for (s = 0; s < S; s++) {
            for (n = 0; n < N; n++) {

                const ISMRMRD::AcquisitionHeader &acq_header = recon_bit.data_.headers_(0, 0, n, s, slc);
                ISMRMRD::ImageHeader &im_header = res.headers_(n, s, slc);
                ISMRMRD::MetaContainer &meta = res.meta_[n + s * N + slc * N * S];

                im_header.version = acq_header.version;
                im_header.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                im_header.flags = acq_header.flags;
                im_header.measurement_uid = acq_header.measurement_uid;

                uint16_t matrix_size[3] = {(uint16_t) RO, (uint16_t) E1, (uint16_t) E2};
                std::copy(matrix_size, std::end(matrix_size), im_header.matrix_size);

                std::copy(std::begin(recon_bit.data_.sampling_.recon_FOV_), std::end(recon_bit.data_.sampling_.recon_FOV_),
                          std::begin(im_header.field_of_view));

                im_header.channels = (uint16_t) CHA;

                std::copy(acq_header.position, std::end(acq_header.position), im_header.position);

                std::copy(acq_header.read_dir, std::end(acq_header.read_dir), im_header.read_dir);

                std::copy(acq_header.phase_dir, std::end(acq_header.phase_dir), im_header.phase_dir);
                std::copy(acq_header.slice_dir, std::end(acq_header.slice_dir), im_header.slice_dir);
                std::copy(acq_header.patient_table_position, std::end(acq_header.patient_table_position),
                          im_header.patient_table_position);


                im_header.average = acq_header.idx.average;
                im_header.slice = acq_header.idx.slice;
                im_header.contrast = acq_header.idx.contrast;
                im_header.phase = acq_header.idx.phase;
                im_header.repetition = acq_header.idx.repetition;
                im_header.set = acq_header.idx.set;

                im_header.acquisition_time_stamp = acq_header.acquisition_time_stamp;

                std::copy(acq_header.physiology_time_stamp, std::end(acq_header.physiology_time_stamp),
                          im_header.physiology_time_stamp);

                im_header.image_type = ISMRMRD::ISMRMRD_IMTYPE_COMPLEX;
                im_header.image_index = (uint16_t) (n + s * N + slc * N * S);
                im_header.image_series_index = 0;

                std::copy(acq_header.user_float, std::end(acq_header.user_float), im_header.user_float);
                std::copy(acq_header.user_int, std::end(acq_header.user_int), im_header.user_int);

                im_header.attribute_string_len = 0;

                meta.set("encoding", (long) encoding);

                meta.set("encoding_FOV", recon_bit.data_.sampling_.encoded_FOV_[0]);
                meta.append("encoding_FOV", recon_bit.data_.sampling_.encoded_FOV_[1]);
                meta.append("encoding_FOV", recon_bit.data_.sampling_.encoded_FOV_[2]);

                meta.set("recon_FOV", recon_bit.data_.sampling_.recon_FOV_[0]);
                meta.append("recon_FOV", recon_bit.data_.sampling_.recon_FOV_[1]);
                meta.append("recon_FOV", recon_bit.data_.sampling_.recon_FOV_[2]);

                meta.set("encoded_matrix", (long) recon_bit.data_.sampling_.encoded_matrix_[0]);
                meta.append("encoded_matrix", (long) recon_bit.data_.sampling_.encoded_matrix_[1]);
                meta.append("encoded_matrix", (long) recon_bit.data_.sampling_.encoded_matrix_[2]);

                meta.set("recon_matrix", (long) recon_bit.data_.sampling_.recon_matrix_[0]);
                meta.append("recon_matrix", (long) recon_bit.data_.sampling_.recon_matrix_[1]);
                meta.append("recon_matrix", (long) recon_bit.data_.sampling_.recon_matrix_[2]);

                meta.set("sampling_limits_RO", (long) recon_bit.data_.sampling_.sampling_limits_[0].min_);
                meta.append("sampling_limits_RO", (long) recon_bit.data_.sampling_.sampling_limits_[0].center_);
                meta.append("sampling_limits_RO", (long) recon_bit.data_.sampling_.sampling_limits_[0].max_);

                meta.set("sampling_limits_E1", (long) recon_bit.data_.sampling_.sampling_limits_[1].min_);
                meta.append("sampling_limits_E1", (long) recon_bit.data_.sampling_.sampling_limits_[1].center_);
                meta.append("sampling_limits_E1", (long) recon_bit.data_.sampling_.sampling_limits_[1].max_);

                meta.set("sampling_limits_E2", (long) recon_bit.data_.sampling_.sampling_limits_[2].min_);
                meta.append("sampling_limits_E2", (long) recon_bit.data_.sampling_.sampling_limits_[2].center_);
                meta.append("sampling_limits_E2", (long) recon_bit.data_.sampling_.sampling_limits_[2].max_);

                meta.set("PatientPosition", (double) res.headers_(n, s, slc).position[0]);
                meta.append("PatientPosition", (double) res.headers_(n, s, slc).position[1]);
                meta.append("PatientPosition", (double) res.headers_(n, s, slc).position[2]);

                meta.set("read_dir", (double) res.headers_(n, s, slc).read_dir[0]);
                meta.append("read_dir", (double) res.headers_(n, s, slc).read_dir[1]);
                meta.append("read_dir", (double) res.headers_(n, s, slc).read_dir[2]);

                meta.set("phase_dir", (double) res.headers_(n, s, slc).phase_dir[0]);
                meta.append("phase_dir", (double) res.headers_(n, s, slc).phase_dir[1]);
                meta.append("phase_dir", (double) res.headers_(n, s, slc).phase_dir[2]);

                meta.set("slice_dir", (double) res.headers_(n, s, slc).slice_dir[0]);
                meta.append("slice_dir", (double) res.headers_(n, s, slc).slice_dir[1]);
                meta.append("slice_dir", (double) res.headers_(n, s, slc).slice_dir[2]);

                meta.set("patient_table_position", (double) res.headers_(n, s, slc).patient_table_position[0]);
                meta.append("patient_table_position", (double) res.headers_(n, s, slc).patient_table_position[1]);
                meta.append("patient_table_position", (double) res.headers_(n, s, slc).patient_table_position[2]);

                meta.set("acquisition_time_stamp", (long) res.headers_(n, s, slc).acquisition_time_stamp);

                meta.set("physiology_time_stamp", (long) res.headers_(n, s, slc).physiology_time_stamp[0]);
                meta.append("physiology_time_stamp", (long) res.headers_(n, s, slc).physiology_time_stamp[1]);
                meta.append("physiology_time_stamp", (long) res.headers_(n, s, slc).physiology_time_stamp[2]);
            }
        }
    }
}

