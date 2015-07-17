
/** \file   mri_core_KarhunenLoeve_Transform.h
    \brief  Implementation Karhunen-Loeve transform (KLT) functionalities
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"

namespace Gadgetron
{
    /// ------------------------------------------------------------------------
    /// KL transform
    /// ------------------------------------------------------------------------
    /// data: M*N data array
    /// the KLT direction is along the N
    /// eigenVectors: N*N eigen vectors, every column is a eigen vector
    /// eigenValues: N*1 eigen values, descending order
    template <typename T> EXPORTMRICORE void eigen_analysis(const hoNDArray<T>& data, hoNDArray<T>& eigenVectors, hoNDArray<T>& eigenValues);

    /// apply the eigen transform
    /// data: M*N data matrix
    /// eigenVectors: N*K eigen vector matrix, every column is a eigen vector
    /// dataEigen: M*K eigen data matrix
    template <typename T> EXPORTMRICORE void apply_eigen_vector(const hoNDArray<T>& data, const hoNDArray<T>& eigenVectors, hoNDArray<T>& dataEigen);

    /// number of kept eigen modes
    /// all modes with eigen values greater than thres*max(eigenValues) are kept
    template <typename T> EXPORTMRICORE void compute_number_of_eigen_modes_kept(const hoNDArray<T>& eigenValues, double thres, size_t& numOfModesKept);

    /// prune the eigen vector matrixes to keep the last numOfModesKept columns
    template <typename T> EXPORTMRICORE void prune_eigen_vector_matrix(const hoNDArray<T>& eigenVectors, size_t numOfModesKept, hoNDArray<T>& eigenVectorsPruned);

    /// ------------------------------------------------------------------------
    /// KL transform coil compression
    /// ------------------------------------------------------------------------
    /// data: the first dimsion is dimension 0
    /// the KL transform is applied along the KLTDim dimension
    /// coeff: square eigen vector matrix
    /// KLTDim : the dimension to perform KL transform
    template <typename T> EXPORTMRICORE void compute_KLT_coeff(const hoNDArray<T>& data, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues, size_t KLTDim);

    /// coeff: CHA*numOfModesKept eigen vector matrix
    /// eigenValues: CHA*1 eigen values
    /// thres <0 or numOfModesKept==-1, keep all modes
    /// if isChaLastDim==true, the CHA is the last dimension
    /// data: [RO E1 CHA ...]
    template <typename T> EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray<T>& data, double thres, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues);
    template <typename T> EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray<T>& data, int numOfModesKept, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues);

    /// data: [RO E1 E2 CHA ...]
    template <typename T> EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray<T>& data, double thres, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues);
    template <typename T> EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray<T>& data, int numOfModesKept, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues);

    /// apply coil compression coefficients
    /// data : [RO E1 srcCHA ...] for 2D or [RO E1 E2 srcCHA ...] for 3D
    /// coeff: [srcCHA dstCHA]
    /// dataEigen: [RO E1 dstCHA ...] or [RO E1 E2 dstCHA ...]
    template <typename T> EXPORTMRICORE void appy_KLT_coil_compression_coeff_2D(const hoNDArray<T>& data, const hoNDArray<T>& coeff, hoNDArray<T>& dataEigen);
    template <typename T> EXPORTMRICORE void appy_KLT_coil_compression_coeff_3D(const hoNDArray<T>& data, const hoNDArray<T>& coeff, hoNDArray<T>& dataEigen);

    /// compute KL transform and perform filtering
    /// the KLT dimension is the last dimension
    template <typename T> EXPORTMRICORE void compute_KLT_filter(const hoNDArray<T>& data, size_t numOfModesKept, hoNDArray<T>& dataKLF);
}
