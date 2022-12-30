
/** \file   mri_core_utility.h
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"
#include "mri_core_data.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"
#include "hoNDKLT.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------------
    /// detect whether a readout has been sampled or not
    // --------------------------------------------------------------------------
    /// data: [RO E1 E2 CHA N S SLC]

    /// sampled: [E1 E2 N S SLC], if a readout is sampled, corresponding sampled location is 1; otherwise, 0
    template <typename T> EXPORTMRICORE hoNDArray<bool> detect_readout_sampling_status(const hoNDArray<T>& data);

    /// detect sampled region along E1
    template <typename T> EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E1(const hoNDArray<T>& data);

    /// detect sampled region along E2
    template <typename T> EXPORTMRICORE std::tuple<size_t, size_t> detect_sampled_region_E2(const hoNDArray<T>& data);

    /// zero padding resize for kspace and complex images
    /// if sizeE2<=1, 2D zero padding resize is performed
    template <typename T> EXPORTMRICORE void zero_pad_resize(const hoNDArray<T>& complexIm, size_t sizeRO, size_t sizeE1, size_t sizeE2, hoNDArray<T>& complexImResized);

    /// compute averaged kspace along N and S dimensions
    /// if count_sampling_freq == true, number of times where a line is sampled along N is counted, instread of an uniformed averaging
    template <typename T> EXPORTMRICORE void compute_averaged_data_N_S(const hoNDArray<T>& data, bool average_N, bool average_S, bool count_sampling_freq, hoNDArray<T>& res);

    /// select specific N and/or S
    template <typename T> EXPORTMRICORE void select_data_N_S(const hoNDArray<T>& data, bool select_N, size_t n, bool select_S, size_t s, hoNDArray<T>& res);

    /// compute KL coefficients from the data
    /// the KLT has the size [N S SLC]
    /// if average_N==true or average_S==true, data will first be averaged along N or S
    /// if coil_compression_thres>0 or compression_num_modesKept>0, the number of kept channels is determine; compression_num_modesKept has the priority if it is set
    /// for all N, S and SLC, the same number of channels is kept. This number is either set by compression_num_modesKept or automatically determined in the first KLT prepare call
    template <typename T> EXPORTMRICORE void compute_eigen_channel_coefficients(const hoNDArray<T>& data, bool average_N, bool average_S, bool count_sampling_freq, size_t N, size_t S, double coil_compression_thres, size_t compression_num_modesKept, std::vector< std::vector< std::vector< hoNDKLT<T> > > >& KLT);

    /// apply eigen channel coefficients
    /// apply KLT coefficients to data for every N, S, and SLC
    template <typename T> EXPORTMRICORE void apply_eigen_channel_coefficients(const std::vector< std::vector< std::vector< hoNDKLT<T> > > >& KLT, hoNDArray<T>& data);

    /// get the path of debug folder
    // environmental variable GADGETRON_DEBUG_FOLDER is used 
    EXPORTMRICORE void get_debug_folder_path(const std::string& debugFolder, std::string& debugFolderPath);

    // find the calibration mode from protocol
    void EXPORTMRICORE find_calib_mode(ISMRMRD::IsmrmrdHeader& h, Gadgetron::ismrmrdCALIBMODE& CalibMode, Gadgetron::IsmrmrdDIM& InterleaveDim, double& acceFactorE1, double& acceFactorE2, bool verbose = false);

    // find the encoding limits from protocol
    void EXPORTMRICORE find_encoding_limits(ISMRMRD::IsmrmrdHeader& h, ISMRMRD::EncodingCounters& meas_max_idx, bool verbose = false);

    // find encoding matrix size and FOV
    void EXPORTMRICORE find_matrix_size_encoding(ISMRMRD::IsmrmrdHeader& h, size_t matrix_size_encoding[3]);
    void EXPORTMRICORE find_FOV_encoding(ISMRMRD::IsmrmrdHeader& h, float field_of_view_encoding[3]);

    // find recon matrix size and FOV
    void EXPORTMRICORE find_matrix_size_recon(ISMRMRD::IsmrmrdHeader& h, size_t matrix_size_recon[3]);
    void EXPORTMRICORE find_FOV_recon(ISMRMRD::IsmrmrdHeader& h, float field_of_view_recon[3]);

    // return the current moment as a string
    void EXPORTMRICORE get_current_moment(std::string& procTime);

    // get a vector of values from ismrmrd meta
    void EXPORTMRICORE get_ismrmrd_meta_values(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<long>& v);
    void EXPORTMRICORE get_ismrmrd_meta_values(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<double>& v);
    void EXPORTMRICORE get_ismrmrd_meta_values(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<std::string>& v);

    template <typename T> EXPORTMRICORE void set_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v);
    void EXPORTMRICORE set_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v);

    template <typename T> EXPORTMRICORE void append_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v);
    void EXPORTMRICORE append_ismrmrd_meta_values(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v);

    // perform the patient to device coordinate transformation
    void EXPORTMRICORE PatientCoordinateSystem_to_DeviceCoordinateSystem(double& x, double& y, double& z, const std::string& position);
    void EXPORTMRICORE DeviceCoordinateSystem_to_PatientCoordinateSystem(double& x, double& y, double& z, const std::string& position);

    // whether to images have identical slice prescription
    // if so, return true; otherwise, return false;
    bool EXPORTMRICORE check_idential_slice_prescription(ISMRMRD::ISMRMRD_ImageHeader a, ISMRMRD::ISMRMRD_ImageHeader b);

    // get ismrmd dim name
    std::string EXPORTMRICORE get_ismrmrd_dim_name(const IsmrmrdDIM& dim);
    // given the name, get the ismrmrd dim
    IsmrmrdDIM EXPORTMRICORE get_ismrmrd_dim_from_name(const std::string& name);

    EXPORTMRICORE std::map<std::string,std::int64_t> to_map(const std::vector<ISMRMRD::UserParameterLong>&);
    EXPORTMRICORE std::map<std::string,double> to_map(const std::vector<ISMRMRD::UserParameterDouble>&);
    EXPORTMRICORE std::map<std::string,std::string> to_map(const std::vector<ISMRMRD::UserParameterString>&);

    ISMRMRD::ImageHeader image_header_from_acquisition(const ISMRMRD::AcquisitionHeader& acq_header,const ISMRMRD::IsmrmrdHeader& header, const hoNDArray<std::complex<float>>& data );
}
