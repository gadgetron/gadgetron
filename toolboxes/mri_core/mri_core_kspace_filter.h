
/** \file   mri_core_kspace_filter.h
    \brief  Implementation kspace filter functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron
{
    /// ------------------------------------------------------------------------
    /// filter types
    /// ------------------------------------------------------------------------
    // define the kspace filter type
    enum MRDKSPACEFILTER
    {
        MRD_FILTER_GAUSSIAN = 0,
        MRD_FILTER_HANNING,
        MRD_FILTER_TAPERED_HANNING,
        MRD_FILTER_NONE
    };

    // get the kspace filter algorithm from name
    MRDKSPACEFILTER get_kspace_filter_type(const std::string& name);
    std::string get_kspace_filter_name(MRDKSPACEFILTER filterType);

    /// ------------------------------------------------------------------------
    /// filter generation
    /// ------------------------------------------------------------------------
    /// generate symmetric filter
    /// sigma: for Gaussian, in the unit of pixel
    /// width: for TaperedHanning filter etc., the length of transition band
    /// filterType: "Gaussian" or "Hanning" or "TaperedHanning" or "None"
    template <typename T> void generate_symmetric_filter(size_t len, hoNDArray<T>& filter, MRDKSPACEFILTER filterType, double sigma = 1.5, size_t width = 15);

    /// generate asymmetric filter, used for partial fourier/asymmetric echo filtering
    /// start, end: the data range within len
    /// filterType: "TaperedHanning" or "None"
    /// if filterType==None and densityComp==true, the 0-1-2 filter will be generated
    /// if filterType==TaperedHanning and the densityComp is true, the density compensation version of tapered filter will be generated
    /// where unacquired region has filter values 0 and symmetric region 1 and nonsymmetric region 2
    /// if densityComp==false, the one side tapered filter will be generated
    template <typename T> void generate_asymmetric_filter(size_t len, size_t start, size_t end, hoNDArray<T>& filter, MRDKSPACEFILTER filterType, size_t width, bool densityComp = false);

    /// generate kspace filter for reference data
    /// a hanning filter is generated for the ref data, to make sure the ref image is free of ringing after zero-padding
    /// start, end: the data range within len
    template <typename T> void generate_symmetric_filter_ref(size_t len, size_t start, size_t end, hoNDArray<T>& filter);

    /// ------------------------------------------------------------------------
    /// applying filter
    /// ------------------------------------------------------------------------

    /// compute 2D filter from two 1D filters
    template <typename T> void compute_2d_filter(const hoNDArray<T>& fx, const hoNDArray<T>& fy, hoNDArray<T>& fxy);
    void compute_2d_filter(const hoNDArray<float>& fx, const hoNDArray<float>& fy, hoNDArray< std::complex<float> >& fxy);
    void compute_2d_filter(const hoNDArray<double>& fx, const hoNDArray<double>& fy, hoNDArray< std::complex<double> >& fxy);

    /// compute 3D filter from three 1D filters
    template <typename T> void compute_3d_filter(const hoNDArray<T>& fx, const hoNDArray<T>& fy, const hoNDArray<T>& fz, hoNDArray<T>& fxyz);
    void compute_3d_filter(const hoNDArray<float>& fx, const hoNDArray<float>& fy, const hoNDArray<float>& fz, hoNDArray< std::complex<float> >& fxyz);
    void compute_3d_filter(const hoNDArray<double>& fx, const hoNDArray<double>& fy, const hoNDArray<double>& fz, hoNDArray< std::complex<double> >& fxyz);

    /// apply kspace filter
    /// data: in kspace, [RO E1 E2 ...]
    /// all functions support in-place operation
    template <typename T> void apply_kspace_filter_RO(hoNDArray<T>& data, const hoNDArray<T>& fRO);
    template <typename T> void apply_kspace_filter_RO(const hoNDArray<T>& data, const hoNDArray<T>& fRO, hoNDArray<T>& dataFiltered);

    template <typename T> void apply_kspace_filter_E1(const hoNDArray<T>& data, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered);

    template <typename T> void apply_kspace_filter_ROE1(const hoNDArray<T>& data, const hoNDArray<T>& fROE1, hoNDArray<T>& dataFiltered);
    template <typename T> void apply_kspace_filter_ROE1(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered);

    template <typename T> void apply_kspace_filter_E2(const hoNDArray<T>& data, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    template <typename T> void apply_kspace_filter_ROE2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    template <typename T> void apply_kspace_filter_E1E2(const hoNDArray<T>& data, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    template <typename T> void apply_kspace_filter_ROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fROE1E2, hoNDArray<T>& dataFiltered);
    template <typename T> void apply_kspace_filter_ROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);

    /// ------------------------------------------------------------------------
    /// filter utility functions
    /// ------------------------------------------------------------------------

    /// find the symmetric sampled region
    /// given the start, end and center
    void find_symmetric_sampled_region(size_t start, size_t end, size_t center, size_t& startSym, size_t& endSym);

    /// compute the filter SNR unit scale factor
    template <typename T> void compute_filter_SNR_unit_scale_factor(const hoNDArray<T>& filter, T& scalFactor);
}
