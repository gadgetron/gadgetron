
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
    /// the KLT direction is along the second dimension
    /// eigenVectors: N*N eigen vectors, every column is a eigen vector
    /// eigenValues: N*1 eigen values, ascending order
    template <typename T> EXPORTMRICORE void eigen_analysis(const hoNDArray<T>& data, hoNDArray<T>& eigenVectors, hoNDArray<T>& eigenValues);

    /// apply the eigen transform
    /// data: M*N data matrix
    /// eigenVectors: N*K eigen vector matrix, every column is a eigen vector
    /// dataEigen: M*K eigen data matrix
    template <typename T> EXPORTMRICORE void apply_eigen_vector(const hoNDArray<T>& data, const hoNDArray<T>& eigenVectors, hoNDArray<T>& dataEigen);

    /// data: the first dimsion is dimension 0
    /// the KL transform is applied along the KLTDim dimension
    /// coeff: square eigen vector matrix
    /// KLTDim : the dimension to perform KL transform
    template <typename T> EXPORTMRICORE void compute_KLT_coeff(const hoNDArray<T>& data, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues, size_t KLTDim);
}
