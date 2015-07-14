
/** \file   mri_core_KarhunenLoeve_Transform.cpp
    \brief  Implementation Karhunen-Loeve Transform (KLT) functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#include "mri_core_KarhunenLoeve_Transform.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_linalg.h"
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
    void apply_eigen_vector(const hoNDArray<T>& data, hoNDArray<T>& dataEigen, const hoNDArray<T>& eigenVectors)
    {
        try
        {
            size_t M = data.get_size(0);
            size_t N = data.get_size(1);

            GADGET_CHECK_THROW(eigenVectors.get_size(0) == N);

            size_t K = eigenVectors.get_size(1);

            dataEigen.create(M, K);
            Gadgetron::clear(dataEigen);

            // M*N multiplies N*K
            Gadgetron::gemm(dataEigen, data, false, eigenVectors, false);
        }
        catch (...)
        {
            GADGET_THROW("Errors in apply_eigen_vector(...) ... ");
        }
    }

    template EXPORTMRICORE void apply_eigen_vector(const hoNDArray<float>& data, hoNDArray<float>& dataEigen, hoNDArray<float>& eigenVectors);
    template EXPORTMRICORE void apply_eigen_vector(const hoNDArray<double>& data, hoNDArray<double>& dataEigen, hoNDArray<double>& eigenVectors);
    template EXPORTMRICORE void apply_eigen_vector(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& dataEigen, hoNDArray< std::complex<float> >& eigenVectors);
    template EXPORTMRICORE void apply_eigen_vector(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& dataEigen, hoNDArray< std::complex<double> >& eigenVectors);

    // ------------------------------------------------------------------------

    template <typename T> 
    void compute_number_of_eigen_modes_kept(const hoNDArray<T>& eigenValues, double thres, long long& numOfModesKept)
    {
    }

    // ------------------------------------------------------------------------

    template <typename T> 
    void prune_eigen_vector_matrix(const hoNDArray<T>& eigenVectors, long long numOfModesKept, hoNDArray<T>& eigenVectorsPruned)
    {
    }

    // ------------------------------------------------------------------------

}
