
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
void compute_number_of_eigen_modes_kept(const hoNDArray<T>& eigenValues, double thres, size_t& numOfModesKept)
{
    try
    {
        size_t M = eigenValues.get_size(0);

        if (thres <= 0)
        {
            numOfModesKept = M;
            return;
        }

        long long m;
        for (m = M - 2; m >= 0; m--)
        {
            if (std::abs(eigenValues(m, 0)) < thres*std::abs(eigenValues(M - 1, 0)))
            {
                break;
            }
        }

        numOfModesKept = M - m - 1;
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_number_of_eigen_modes_kept(...) ... ");
    }
}

template EXPORTMRICORE void compute_number_of_eigen_modes_kept(const hoNDArray<float>& eigenValues, double thres, size_t& numOfModesKept);
template EXPORTMRICORE void compute_number_of_eigen_modes_kept(const hoNDArray<double>& eigenValues, double thres, size_t& numOfModesKept);
template EXPORTMRICORE void compute_number_of_eigen_modes_kept(const hoNDArray< std::complex<float> >& eigenValues, double thres, size_t& numOfModesKept);
template EXPORTMRICORE void compute_number_of_eigen_modes_kept(const hoNDArray< std::complex<double> >& eigenValues, double thres, size_t& numOfModesKept);

// ------------------------------------------------------------------------

template <typename T> 
void prune_eigen_vector_matrix(const hoNDArray<T>& eigenVectors, size_t numOfModesKept, hoNDArray<T>& eigenVectorsPruned)
{
    try
    {
        size_t M = eigenVectors.get_size(0);
        size_t N = eigenVectors.get_size(1);

        if (numOfModesKept <= 0 || numOfModesKept>N)
        {
            GWARN_STREAM("prune_eigen_vector_matrix(...) - numOfModesKept<=0 || numOfModesKept>N : " << numOfModesKept);
            eigenVectorsPruned = eigenVectors;
            return;
        }

        eigenVectorsPruned.create(M, numOfModesKept);

        size_t n;
        for (n = N - numOfModesKept; n < N; n++)
        {
            memcpy(eigenVectorsPruned.begin() + (n - N + numOfModesKept)*M, eigenVectors.begin() + n*M, sizeof(T)*M);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in prune_eigen_vector_matrix(...) ... ");
    }
}

template EXPORTMRICORE void prune_eigen_vector_matrix(const hoNDArray<float>& eigenVectors, size_t numOfModesKept, hoNDArray<float>& eigenVectorsPruned);
template EXPORTMRICORE void prune_eigen_vector_matrix(const hoNDArray<double>& eigenVectors, size_t numOfModesKept, hoNDArray<double>& eigenVectorsPruned);
template EXPORTMRICORE void prune_eigen_vector_matrix(const hoNDArray< std::complex<float> >& eigenVectors, size_t numOfModesKept, hoNDArray< std::complex<float> >& eigenVectorsPruned);
template EXPORTMRICORE void prune_eigen_vector_matrix(const hoNDArray< std::complex<double> >& eigenVectors, size_t numOfModesKept, hoNDArray< std::complex<double> >& eigenVectorsPruned);

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

// ------------------------------------------------------------------------

template <typename T> 
void compute_KLT_coil_compression_coeff_2D(const hoNDArray<T>& data, double thres, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(NDim >= 3);

        hoNDArray<T> eigenVectors;
        compute_KLT_coeff(data, eigenVectors, eigenValues, 2); // RO, E1, CHA

        size_t numOfModesKept;
        compute_number_of_eigen_modes_kept(eigenValues, thres, numOfModesKept);

        prune_eigen_vector_matrix(eigenVectors, numOfModesKept, coeff);
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_KLT_coil_compression_coeff_2D(thres) ... ");
    }
}

template EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray<float>& data, double thres, hoNDArray<float>& coeff, hoNDArray<float>& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray<double>& data, double thres, hoNDArray<double>& coeff, hoNDArray<double>& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray< std::complex<float> >& data, double thres, hoNDArray< std::complex<float> >& coeff, hoNDArray< std::complex<float> >& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray< std::complex<double> >& data, double thres, hoNDArray< std::complex<double> >& coeff, hoNDArray< std::complex<double> >& eigenValues);

// ------------------------------------------------------------------------

template <typename T> 
void compute_KLT_coil_compression_coeff_2D(const hoNDArray<T>& data, size_t numOfModesKept, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(NDim >= 3);

        hoNDArray<T> eigenVectors;
        compute_KLT_coeff(data, eigenVectors, eigenValues, 2); // RO, E1, CHA

        size_t CHA = data.get_size(2);

        if (numOfModesKept >= CHA - 1)
        {
            coeff = eigenVectors;
            return;
        }

        if (numOfModesKept == 0) numOfModesKept = 1;

        prune_eigen_vector_matrix(eigenVectors, numOfModesKept, coeff);
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_KLT_coil_compression_coeff_2D(numOfModesKept) ... ");
    }
}

template EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray<float>& data, size_t numOfModesKept, hoNDArray<float>& coeff, hoNDArray<float>& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray<double>& data, size_t numOfModesKept, hoNDArray<double>& coeff, hoNDArray<double>& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray< std::complex<float> >& data, size_t numOfModesKept, hoNDArray< std::complex<float> >& coeff, hoNDArray< std::complex<float> >& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_2D(const hoNDArray< std::complex<double> >& data, size_t numOfModesKept, hoNDArray< std::complex<double> >& coeff, hoNDArray< std::complex<double> >& eigenValues);

// ------------------------------------------------------------------------

/// data: [RO E1 E2 CHA ...]
template <typename T> 
void compute_KLT_coil_compression_coeff_3D(const hoNDArray<T>& data, double thres, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(NDim >= 4);

        hoNDArray<T> eigenVectors;
        compute_KLT_coeff(data, eigenVectors, eigenValues, 3); // RO, E1, E2, CHA

        size_t numOfModesKept;
        compute_number_of_eigen_modes_kept(eigenValues, thres, numOfModesKept);

        prune_eigen_vector_matrix(eigenVectors, numOfModesKept, coeff);
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_KLT_coil_compression_coeff_3D(thres) ... ");
    }
}

template EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray<float>& data, double thres, hoNDArray<float>& coeff, hoNDArray<float>& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray<double>& data, double thres, hoNDArray<double>& coeff, hoNDArray<double>& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray< std::complex<float> >& data, double thres, hoNDArray< std::complex<float> >& coeff, hoNDArray< std::complex<float> >& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray< std::complex<double> >& data, double thres, hoNDArray< std::complex<double> >& coeff, hoNDArray< std::complex<double> >& eigenValues);

// ------------------------------------------------------------------------

template <typename T> 
void compute_KLT_coil_compression_coeff_3D(const hoNDArray<T>& data, size_t numOfModesKept, hoNDArray<T>& coeff, hoNDArray<T>& eigenValues)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(NDim >= 4);

        hoNDArray<T> eigenVectors;
        compute_KLT_coeff(data, eigenVectors, eigenValues, 3); // RO, E1, E2, CHA

        size_t CHA = data.get_size(3);

        if (numOfModesKept >= CHA - 1)
        {
            coeff = eigenVectors;
            return;
        }

        if (numOfModesKept == 0) numOfModesKept = 1;

        prune_eigen_vector_matrix(eigenVectors, numOfModesKept, coeff);
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_KLT_coil_compression_coeff_3D(numOfModesKept) ... ");
    }
}

template EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray<float>& data, size_t numOfModesKept, hoNDArray<float>& coeff, hoNDArray<float>& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray<double>& data, size_t numOfModesKept, hoNDArray<double>& coeff, hoNDArray<double>& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray< std::complex<float> >& data, size_t numOfModesKept, hoNDArray< std::complex<float> >& coeff, hoNDArray< std::complex<float> >& eigenValues);
template EXPORTMRICORE void compute_KLT_coil_compression_coeff_3D(const hoNDArray< std::complex<double> >& data, size_t numOfModesKept, hoNDArray< std::complex<double> >& coeff, hoNDArray< std::complex<double> >& eigenValues);

// ------------------------------------------------------------------------

template <typename T> 
void appy_KLT_coil_compression_coeff_2D(const hoNDArray<T>& data, const hoNDArray<T>& coeff, hoNDArray<T>& dataEigen)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(NDim >= 3);

        std::vector<size_t> dim;
        data.get_dimensions(dim);

        size_t RO = dim[0];
        size_t E1 = dim[1];

        size_t CHA = data.get_size(2);
        size_t N = data.get_number_of_elements() / (RO*E1*CHA);

        size_t dstCHA = coeff.get_size(1);

        std::vector<size_t> dimEigen(dim);
        dimEigen[2] = dstCHA;

        dataEigen.create(dimEigen);

        long long n;

#pragma omp parallel for default(none) private(n) shared(N, RO, E1, CHA, dstCHA, data, dataEigen, coeff) if(N>1)
        for (n = 0; n < (long long)N; n++)
        {
            hoNDArray<T> data2D(RO*E1, CHA, const_cast<T*>(data.begin()) + n*RO*E1*CHA);
            hoNDArray<T> dataEigen2D(RO*E1, dstCHA, const_cast<T*>(dataEigen.begin()) + n*RO*E1*dstCHA);

            Gadgetron::apply_eigen_vector(data2D, coeff, dataEigen2D);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in appy_KLT_coil_compression_coeff_2D(...) ... ");
    }
}

template EXPORTMRICORE void appy_KLT_coil_compression_coeff_2D(const hoNDArray<float>& data, const hoNDArray<float>& coeff, hoNDArray<float>& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_2D(const hoNDArray<double>& data, const hoNDArray<double>& coeff, hoNDArray<double>& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_2D(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& coeff, hoNDArray< std::complex<float> >& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_2D(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& coeff, hoNDArray< std::complex<double> >& dataEigen);

// ------------------------------------------------------------------------

template <typename T> 
void appy_KLT_coil_compression_coeff_3D(const hoNDArray<T>& data, const hoNDArray<T>& coeff, hoNDArray<T>& dataEigen)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(NDim >= 4);

        std::vector<size_t> dim;
        data.get_dimensions(dim);

        size_t RO = dim[0];
        size_t E1 = dim[1];
        size_t E2 = dim[2];

        size_t CHA = data.get_size(3);
        size_t N = data.get_number_of_elements() / (RO*E1*E2*CHA);

        size_t dstCHA = coeff.get_size(1);

        std::vector<size_t> dimEigen(dim);
        dimEigen[3] = dstCHA;

        dataEigen.create(dimEigen);

        long long n;

#pragma omp parallel for default(none) private(n) shared(N, RO, E1, E2, CHA, dstCHA, data, dataEigen, coeff) if(N>1)
        for (n = 0; n < (long long)N; n++)
        {
            hoNDArray<T> data3D(RO*E1*E2, CHA, const_cast<T*>(data.begin()) + n*RO*E1*E2*CHA);
            hoNDArray<T> dataEigen3D(RO*E1*E2, dstCHA, const_cast<T*>(dataEigen.begin()) + n*RO*E1*E2*dstCHA);

            Gadgetron::apply_eigen_vector(data3D, coeff, dataEigen3D);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in appy_KLT_coil_compression_coeff_3D(...) ... ");
    }
}

template EXPORTMRICORE void appy_KLT_coil_compression_coeff_3D(const hoNDArray<float>& data, const hoNDArray<float>& coeff, hoNDArray<float>& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_3D(const hoNDArray<double>& data, const hoNDArray<double>& coeff, hoNDArray<double>& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_3D(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& coeff, hoNDArray< std::complex<float> >& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_3D(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& coeff, hoNDArray< std::complex<double> >& dataEigen);

// ------------------------------------------------------------------------

template <typename T> 
void appy_KLT_coil_compression_coeff_2D(const hoNDArray<T>& data, const std::vector< hoNDArray<T> >& coeff, hoNDArray<T>& dataEigen)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(NDim >= 3);

        GADGET_CHECK_THROW(coeff.size() >= data.get_size(NDim - 1));

        size_t LastDim = coeff.size();
        size_t dstCHA = coeff[0].get_size(1);

        size_t n;
        for (n = 1; n<LastDim; n++)
        {
            GADGET_CHECK_THROW(coeff[n].get_size(1) == dstCHA);
        }

        size_t LastDimData = data.get_size(NDim - 1);
        std::vector<size_t> dim;
        data.get_dimensions(dim);
        long long N = data.get_number_of_elements() / LastDimData;

        std::vector<size_t> dimEigen(dim);
        dimEigen[2] = dstCHA;

        dataEigen.create(&dimEigen);
        long long eigenN = dataEigen.get_number_of_elements() / LastDimData;

        std::vector<size_t> dimLastDim(NDim - 1);
        for (n = 0; n<NDim - 1; n++)
        {
            dimLastDim[n] = dim[n];
        }

        std::vector<size_t> dimEigenLastDim(dimLastDim);
        dimEigenLastDim[2] = dstCHA;

        if (LastDimData>1)
        {
            hoNDArray<T> dataEigenLastDim;
            for (n = 0; n < LastDimData; n++)
            {
                hoNDArray<T> dataLastDim(&dimLastDim, const_cast<T*>(data.begin() + n*N));
                Gadgetron::appy_KLT_coil_compression_coeff_2D(dataLastDim, coeff[n], dataEigenLastDim);
                memcpy(dataEigen.begin() + n*eigenN, dataEigenLastDim.begin(), dataEigenLastDim.get_number_of_bytes());
            }
        }
        else
        {
            hoNDArray<T> dataLastDim(&dimLastDim, const_cast<T*>(data.begin()));
            hoNDArray<T> dataEigenLastDim(&dimEigenLastDim, dataEigen.begin());
            Gadgetron::appy_KLT_coil_compression_coeff_2D(dataLastDim, coeff[n], dataEigenLastDim);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in appy_KLT_coil_compression_coeff_2D(std::vector<hoNDArray<T> >& coeff) ... ");
    }
}

template EXPORTMRICORE void appy_KLT_coil_compression_coeff_2D(const hoNDArray<float>& data, const std::vector< hoNDArray<float> >& coeff, hoNDArray<float>& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_2D(const hoNDArray<double>& data, const std::vector< hoNDArray<double> >& coeff, hoNDArray<double>& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_2D(const hoNDArray< std::complex<float> >& data, const std::vector< hoNDArray< std::complex<float> > >& coeff, hoNDArray< std::complex<float> >& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_2D(const hoNDArray< std::complex<double> >& data, const std::vector< hoNDArray< std::complex<double> > >& coeff, hoNDArray< std::complex<double> >& dataEigen);

// ------------------------------------------------------------------------

template <typename T> 
void appy_KLT_coil_compression_coeff_3D(const hoNDArray<T>& data, const std::vector< hoNDArray<T> >& coeff, hoNDArray<T>& dataEigen)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_THROW(NDim >= 4);

        GADGET_CHECK_THROW(coeff.size() >= data.get_size(NDim - 1));

        size_t LastDim = coeff.size();
        size_t dstCHA = coeff[0].get_size(1);

        size_t n;
        for (n = 1; n<LastDim; n++)
        {
            GADGET_CHECK_THROW(coeff[n].get_size(1) == dstCHA);
        }

        size_t LastDimData = data.get_size(NDim - 1);
        std::vector<size_t> dim;
        data.get_dimensions(dim);
        long long N = data.get_number_of_elements() / LastDimData;

        std::vector<size_t> dimEigen(dim);
        dimEigen[3] = dstCHA;

        dataEigen.create(&dimEigen);
        long long eigenN = dataEigen.get_number_of_elements() / LastDimData;

        std::vector<size_t> dimLastDim(NDim - 1);
        for (n = 0; n<NDim - 1; n++)
        {
            dimLastDim[n] = dim[n];
        }

        std::vector<size_t> dimEigenLastDim(dimLastDim);
        dimEigenLastDim[3] = dstCHA;

        if (LastDimData>1)
        {
            hoNDArray<T> dataEigenLastDim;
            for (n = 0; n<LastDimData; n++)
            {
                hoNDArray<T> dataLastDim(&dimLastDim, const_cast<T*>(data.begin() + n*N));
                Gadgetron::appy_KLT_coil_compression_coeff_3D(dataLastDim, coeff[n], dataEigenLastDim);
                memcpy(dataEigen.begin() + n*eigenN, dataEigenLastDim.begin(), dataEigenLastDim.get_number_of_bytes());
            }
        }
        else
        {
            hoNDArray<T> dataLastDim(&dimLastDim, const_cast<T*>(data.begin()));
            hoNDArray<T> dataEigenLastDim(&dimEigenLastDim, dataEigen.begin());
            Gadgetron::appy_KLT_coil_compression_coeff_2D(dataLastDim, coeff[n], dataEigenLastDim);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in appy_KLT_coil_compression_coeff_3D(std::vector<hoNDArray<T> >& coeff) ... ");
    }
}

template EXPORTMRICORE void appy_KLT_coil_compression_coeff_3D(const hoNDArray<float>& data, const std::vector< hoNDArray<float> >& coeff, hoNDArray<float>& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_3D(const hoNDArray<double>& data, const std::vector< hoNDArray<double> >& coeff, hoNDArray<double>& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_3D(const hoNDArray< std::complex<float> >& data, const std::vector< hoNDArray< std::complex<float> > >& coeff, hoNDArray< std::complex<float> >& dataEigen);
template EXPORTMRICORE void appy_KLT_coil_compression_coeff_3D(const hoNDArray< std::complex<double> >& data, const std::vector< hoNDArray< std::complex<double> > >& coeff, hoNDArray< std::complex<double> >& dataEigen);

// ------------------------------------------------------------------------

template <typename T> 
void compute_KLT_filter(const hoNDArray<T>& data, size_t numOfModesKept, hoNDArray<T>& dataKLF)
{
    try
    {
        if (!dataKLF.dimensions_equal(&data))
        {
            dataKLF = data;
        }

        size_t NDim = data.get_number_of_dimensions();
        size_t M = data.get_size(NDim - 1);
        size_t N = data.get_number_of_elements() / M;

        if (numOfModesKept > M) numOfModesKept = M;

        hoNDArray<T> A(N, M, const_cast<T*>(data.begin()));

        hoMatrix<T> eigenVectors;
        hoNDArray<T> eigenValues;
        eigen_analysis(A, eigenVectors, eigenValues);

        hoNDArray<T> E(eigenVectors);
        size_t r, c;
        for (c = 0; c<M - numOfModesKept + 1; c++)
        {
            for (r = 0; r<M; r++)
            {
                E(r, c) = T(0);
            }
        }

        hoMatrix<T> ET;
        ET.createMatrix(M, M);
        memcpy(ET.begin(), eigenVectors.begin(), sizeof(T)*M*M);

        Gadgetron::conjugatetrans(eigenVectors, ET);

        hoMatrix<T> EET(M, M);
        Gadgetron::clear(EET);
        Gadgetron::gemm(EET, E, false, ET, false);

        hoMatrix<T> R(N, M, dataKLF.begin());
        Gadgetron::gemm(R, A, false, EET, false);
    }
    catch (...)
    {
        GADGET_THROW("Errors in compute_KLT_filter(...) ... ");
    }
}

template EXPORTMRICORE void compute_KLT_filter(const hoNDArray<float>& data, size_t numOfModesKept, hoNDArray<float>& dataKLF);
template EXPORTMRICORE void compute_KLT_filter(const hoNDArray<double>& data, size_t numOfModesKept, hoNDArray<double>& dataKLF);
template EXPORTMRICORE void compute_KLT_filter(const hoNDArray< std::complex<float> >& data, size_t numOfModesKept, hoNDArray< std::complex<float> >& dataKLF);
template EXPORTMRICORE void compute_KLT_filter(const hoNDArray< std::complex<double> >& data, size_t numOfModesKept, hoNDArray< std::complex<double> >& dataKLF);

// ------------------------------------------------------------------------
}
