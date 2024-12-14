
/** \file   mri_core_utility.h
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "hoNDArray.h"
#include "hoNDKLT.h"

#include "mrd/types.h"


namespace Gadgetron
{
    // --------------------------------------------------------------------------
    /// detect whether a readout has been sampled or not
    // --------------------------------------------------------------------------
    /// data: [RO E1 E2 CHA N S SLC]

    /// sampled: [E1 E2 N S SLC], if a readout is sampled, corresponding sampled location is 1; otherwise, 0
    template <typename T> hoNDArray<bool> detect_readout_sampling_status(const hoNDArray<T>& data);

    /// detect sampled region along E1
    template <typename T> std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray<T>& data);

    /// detect sampled region along E2
    template <typename T> std::tuple<size_t, size_t> detect_sampled_region_E2(const hoNDArray<T>& data);

    /// zero padding resize for kspace and complex images
    /// if sizeE2<=1, 2D zero padding resize is performed
    template <typename T> void zero_pad_resize(const hoNDArray<T>& complexIm, size_t sizeRO, size_t sizeE1, size_t sizeE2, hoNDArray<T>& complexImResized);

    /// compute averaged kspace along N and S dimensions
    /// if count_sampling_freq == true, number of times where a line is sampled along N is counted, instread of an uniformed averaging
    template <typename T> void compute_averaged_data_N_S(const hoNDArray<T>& data, bool average_N, bool average_S, bool count_sampling_freq, hoNDArray<T>& res);

    /// select specific N and/or S
    template <typename T> void select_data_N_S(const hoNDArray<T>& data, bool select_N, size_t n, bool select_S, size_t s, hoNDArray<T>& res);

    /// compute KL coefficients from the data
    /// the KLT has the size [N S SLC]
    /// if average_N==true or average_S==true, data will first be averaged along N or S
    /// if coil_compression_thres>0 or compression_num_modesKept>0, the number of kept channels is determine; compression_num_modesKept has the priority if it is set
    /// for all N, S and SLC, the same number of channels is kept. This number is either set by compression_num_modesKept or automatically determined in the first KLT prepare call
    template <typename T> void compute_eigen_channel_coefficients(const hoNDArray<T>& data, bool average_N, bool average_S, bool count_sampling_freq, size_t N, size_t S, double coil_compression_thres, size_t compression_num_modesKept, std::vector< std::vector< std::vector< hoNDKLT<T> > > >& KLT);

    /// apply eigen channel coefficients
    /// apply KLT coefficients to data for every N, S, and SLC
    template <typename T> void apply_eigen_channel_coefficients(const std::vector< std::vector< std::vector< hoNDKLT<T> > > >& KLT, hoNDArray<T>& data);

    /// get the path of debug folder
    // environmental variable GADGETRON_DEBUG_FOLDER is used
    void get_debug_folder_path(const std::string& debugFolder, std::string& debugFolderPath);

    // find the encoding limits from protocol
    void find_encoding_limits(mrd::Header& h, mrd::EncodingCounters& meas_max_idx, bool verbose);

    // find encoding matrix size and FOV
    void find_matrix_size_encoding(mrd::Header& h, size_t matrix_size_encoding[3]);
    void find_FOV_encoding(mrd::Header& h, float field_of_view_encoding[3]);

    // find recon matrix size and FOV
    void find_matrix_size_recon(mrd::Header& h, size_t matrix_size_recon[3]);
    void find_FOV_recon(mrd::Header& h, float field_of_view_recon[3]);

    // return the current moment as a string
    void get_current_moment(std::string& procTime);

    // set common meta attributes from MRD header
    void set_meta_from_mrd_header(const mrd::ImageHeader& header, mrd::ImageMeta& attrib);

    // get a vector of values from MRD ImageMeta
    template <typename T> void get_mrd_meta_values(const mrd::ImageMeta& attrib, const std::string& name, std::vector<T>& v);

    template <typename T> void set_mrd_meta_values(mrd::ImageMeta& attrib, const std::string& name, const std::vector<T>& v);

    template <typename T> void append_mrd_meta_values(mrd::ImageMeta& attrib, const std::string& name, const std::vector<T>& v);

    // perform the patient to device coordinate transformation
    void PatientCoordinateSystem_to_DeviceCoordinateSystem(double& x, double& y, double& z, const std::string& position);
    void DeviceCoordinateSystem_to_PatientCoordinateSystem(double& x, double& y, double& z, const std::string& position);

    // whether to images have identical slice prescription
    // if so, return true; otherwise, return false;
    bool check_idential_slice_prescription(mrd::ImageHeader a, mrd::ImageHeader b);

    std::map<std::string,std::int64_t> to_map(const std::vector<mrd::UserParameterLongType>&);
    std::map<std::string,double> to_map(const std::vector<mrd::UserParameterDoubleType>&);
    std::map<std::string,std::string> to_map(const std::vector<mrd::UserParameterStringType>&);

    mrd::ImageHeader image_header_from_acquisition(const mrd::AcquisitionHeader& acq_header, const mrd::Header& header);

    void add_acquisition_to_bucket(mrd::AcquisitionBucket& bucket, mrd::Acquisition acq);
}
