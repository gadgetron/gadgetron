
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
    template <typename T> EXPORTMRICORE void apply_eigen_vector(const hoNDArray<T>& data, hoNDArray<T>& dataEigen, const hoNDArray<T>& eigenVectors);

    /// number of kept eigen modes
    /// all modes with eigen values greater than thres*max(eigenValues) are kept
    template <typename T> EXPORTMRICORE void compute_number_of_eigen_modes_kept(const hoNDArray<T>& eigenValues, double thres, long long& numOfModesKept);

    /// prune the eigen vector matrixes to keep the last numOfModesKept columns
    template <typename T> EXPORTMRICORE void prune_eigen_vector_matrix(const hoNDArray<T>& eigenVectors, long long numOfModesKept, hoNDArray<T>& eigenVectorsPruned);

    /// ------------------------------------------------------------------------
    /// KL transform coil compression
    /// ------------------------------------------------------------------------
    /// data: at least 3D [RO E1 CHA ...]
    /// the KL transform is applied along CHA
    /// coeff: CHA*numOfModesKept eigen vector matrix
    /// eigenValues: CHA*1 eigen values
    /// thres <0 or numOfModesKept==-1, keep all modes
    /// if isChaLastDim==true, the CHA is the last dimension
    template <typename T> EXPORTMRICORE void computeKLCoilCompressionCoeff(const hoNDArray<T>& data, double thres, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues, bool isChaLastDim = false);
    template <typename T> EXPORTMRICORE void computeKLCoilCompressionCoeff(const hoNDArray<T>& data, int numOfModesKept, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues, bool isChaLastDim = false);
    /// coeff: CHA*CHA eigen vector matrix
    template <typename T> EXPORTMRICORE void computeKLTCoeff(const hoNDArray<T>& data, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues, bool isChaLastDim = false);

    /// dataEigen: [RO E1 numOfModesKept ...] 
    template <typename T> EXPORTMRICORE void computeKLCoilCompression(const hoNDArray<T>& data, double thres, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues, hoNDArray<T>& dataEigen, bool isChaLastDim = false);
    template <typename T> EXPORTMRICORE void computeKLCoilCompression(const hoNDArray<T>& data, int numOfModesKept, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues, hoNDArray<T>& dataEigen, bool isChaLastDim = false);

    /// apply coil compression coefficients
    template <typename T> EXPORTMRICORE void appyKLCoilCompressionCoeff(const hoNDArray<T>& data, const hoNDArray<T>& coeff, hoNDArray<T>& dataEigen, bool isChaLastDim = false);

    /// apply coil compression coefficients on array [RO E1 srcCHA ...]
    /// dataEigen: [RO E1 dstCHA ...]
    /// coeff: [srcCHA dstCHA] matrixes for every last dimension
    /// every last dimension has different compression coefficients
    template <typename T> EXPORTMRICORE void applyKLCoilCompressionCoeff(const hoNDArray<T>& data, const std::vector<hoNDArray<T> >& coeff, hoNDArray<T>& dataEigen, bool isChaLastDim = false);

    /// compute KL transform and perform filtering
    /// the KL dimension is the last dimension
    template <typename T> EXPORTMRICORE void computeKLFilter(const hoNDArray<T>& data, size_t numOfModesKept, hoNDArray<T>& dataKLF);
}
