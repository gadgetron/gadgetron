
/** \file   mri_core_KarhunenLoeve_Transform.cpp
    \brief  Implementation Karhunen-Loeve Transform (KLT) functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_KarhunenLoeve_Transform.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"
#include "hoMatrix.h"

namespace Gadgetron
{

template <typename T> 
void eigen_analysis(const hoNDArray<T>& data, hoNDArray<T>& eigenVectors, hoNDArray<T>& eigenValues)
{
    try
    {
        typedef typename realType<T>::Type ValueType;

        size_t M = data.get_size(0);
        size_t N = data.get_size(1);

        eigenVectors.create(N, N);
        eigenValues.create(N, 1);
        Gadgetron::clear(eigenVectors);
        Gadgetron::clear(eigenValues);

        char uplo = 'L';
        bool isAHA = true;
        Gadgetron::herk(eigenVectors, data, uplo, isAHA);

        hoMatrix<T> m;
        m.createMatrix(N, N, eigenVectors.begin());
        m.copyLowerTriToUpper();

        hoMatrix<T> mean(N, 1);

        hoNDArray<T> dataMean(1, N);
        Gadgetron::sum_over_dimension(data, dataMean, 0);
        memcpy(mean.begin(), dataMean.begin(), sizeof(T)*N);

        Gadgetron::scal((ValueType)1.0 / M, mean);

        hoNDArray<T> MMH(N, N);
        Gadgetron::clear(MMH);

        hoNDArray<T> meanCopy(mean);
        Gadgetron::gemm(MMH, meanCopy, false, mean, true);
        Gadgetron::scal((ValueType)M, MMH);
        Gadgetron::subtract(eigenVectors, MMH, eigenVectors);
        Gadgetron::scal((ValueType)1.0 / (M - 1), eigenVectors);

        hoMatrix<T> EH(N, N);
        conjugatetrans(m, EH);
        Gadgetron::add(eigenVectors, EH, eigenVectors);
        Gadgetron::scal((ValueType)(0.5), eigenVectors);

        Gadgetron::heev(eigenVectors, eigenValues);
    }
    catch (...)
    {
        GADGET_THROW("Errors in eigen_analysis(...) ... ");
    }
}

template EXPORTMRICORE void eigen_analysis(const hoNDArray<float>& data, hoNDArray<float>& eigenVectors, hoNDArray<float>& eigenValues);
template EXPORTMRICORE void eigen_analysis(const hoNDArray<double>& data, hoNDArray<double>& eigenVectors, hoNDArray<double>& eigenValues);
template EXPORTMRICORE void eigen_analysis(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& eigenVectors, hoNDArray< std::complex<float> >& eigenValues);
template EXPORTMRICORE void eigen_analysis(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& eigenVectors, hoNDArray< std::complex<double> >& eigenValues);

// ------------------------------------------------------------------------

template <typename T> 
void apply_eigen_vector(const hoNDArray<T>& data, const hoNDArray<T>& eigenVectors, hoNDArray<T>& dataEigen)
{
    try
    {
        size_t M = data.get_size(0);
        size_t N = data.get_size(1);

        GADGET_CHECK_THROW(eigenVectors.get_size(0) == N);

        size_t K = eigenVectors.get_size(1);

        if ((dataEigen.get_size(0) != M) || (dataEigen.get_size(1) != K))
        {
            dataEigen.create(M, K);
        }

        // M*N multiplies N*K
        Gadgetron::gemm(dataEigen, data, false, eigenVectors, false);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_eigen_vector(...) ... ");
    }
}

template EXPORTMRICORE void apply_eigen_vector(const hoNDArray<float>& data, const hoNDArray<float>& eigenVectors, hoNDArray<float>& dataEigen);
template EXPORTMRICORE void apply_eigen_vector(const hoNDArray<double>& data, const hoNDArray<double>& eigenVectors, hoNDArray<double>& dataEigen);
template EXPORTMRICORE void apply_eigen_vector(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& eigenVectors, hoNDArray< std::complex<float> >& dataEigen);
template EXPORTMRICORE void apply_eigen_vector(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& eigenVectors, hoNDArray< std::complex<double> >& dataEigen);

// ------------------------------------------------------------------------

template <typename T> 
void compute_KLT_coeff(const hoNDArray<T>& data, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues, size_t KLTDim)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(NDim >= 2);

        hoNDArray<T> eigenVectors;
        hoNDArray<T> A;

        if (KLTDim==NDim-1) // if the KLT dimension is the last dimension
        {
            size_t CHA = data.get_size(NDim - 1);
            size_t N = data.get_number_of_elements() / CHA;

            A.create(N, CHA, const_cast<T*>(data.begin()));
            Gadgetron::eigen_analysis(A, coeff, eigenValues);
        }
        else
        {
            std::vector<size_t> dimOrder(NDim);

            size_t l;
            for (l = 0; l<KLTDim; l++)
            {
                dimOrder[l] = l;
            }

            for (l = KLTDim; l<NDim-1; l++)
            {
                dimOrder[l] = l+1;
            }

            dimOrder[NDim - 1] = KLTDim;

            hoNDArray<T> dataP(data);
            Gadgetron::permute(&dataP, &dimOrder);

            size_t CHA = data.get_size(KLTDim);
            size_t N = data.get_number_of_elements() / CHA;
            A.create(N, CHA, dataP.begin());

            Gadgetron::eigen_analysis(A, coeff, eigenValues);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_KLT_coeff(...) ... ");
    }
}

template EXPORTMRICORE void compute_KLT_coeff(const hoNDArray<float>& data, hoNDArray<float>& coeff, hoNDArray<float>& eigenValues, size_t KLTDim);
template EXPORTMRICORE void compute_KLT_coeff(const hoNDArray<double>& data, hoNDArray<double>& coeff, hoNDArray<double>& eigenValues, size_t KLTDim);
template EXPORTMRICORE void compute_KLT_coeff(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coeff, hoNDArray< std::complex<float> >& eigenValues, size_t KLTDim);
template EXPORTMRICORE void compute_KLT_coeff(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coeff, hoNDArray< std::complex<double> >& eigenValues, size_t KLTDim);

}
