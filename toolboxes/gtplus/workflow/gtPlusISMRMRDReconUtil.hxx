
#include "gtPlusISMRMRDReconUtil.h"
#include "hoNDArray_elemwise.h"
#include <algorithm>

#define GT_IMAGING_GEOMETRY_DELTA 0.001

namespace Gadgetron { namespace gtPlus {

template <typename T> 
gtPlusISMRMRDReconUtil<T>::gtPlusISMRMRDReconUtil() {}

template <typename T> 
gtPlusISMRMRDReconUtil<T>::~gtPlusISMRMRDReconUtil() {}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::printInfo(std::ostream& os)
{
    using namespace std;
    os << "-------------- GTPlus ISMRMRD Recon Util -------------" << endl;
    os << "Implementation of recon utilities for ISMRMRD format" << endl;
    os << "------------------------------------------------------" << endl;
}

// ------------------------------------------------------------------------
// coil compression and KarhunenLoeverTransform
// ------------------------------------------------------------------------
template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
KLT_eigenAnalysis(const hoMatrix<T>& data, hoMatrix<T>& eigenVectors, hoMatrix<T>& eigenValues)
{
    try
    {
        typedef typename realType<T>::Type ValueType;

        size_t M = data.rows();
        size_t N = data.cols();

        GADGET_CHECK_RETURN_FALSE(eigenVectors.createMatrix(N, N));
        GADGET_CHECK_RETURN_FALSE(eigenValues.createMatrix(N, 1));
        Gadgetron::clear(eigenVectors);
        Gadgetron::clear(eigenValues);

        //hoMatrix<T> dataCopy(data);
        //GADGET_CHECK_RETURN_FALSE(Gadgetron::gemm(eigenVectors, data, true, dataCopy, false));

        char uplo = 'L';
        bool isAHA = true;
        Gadgetron::herk(eigenVectors, data, uplo, isAHA);
        eigenVectors.copyLowerTriToUpper();

        //eigenVectors.print(std::cout);

        hoMatrix<T> mean(N, 1);
        GADGET_CHECK_RETURN_FALSE(data.sumOverCol(mean));
        Gadgetron::scal((ValueType)1.0/M, mean);

        //mean.print(std::cout);

        hoMatrix<T> MMH(N, N);
        Gadgetron::clear(MMH);

        hoMatrix<T> meanCopy(mean);
        Gadgetron::gemm(MMH, meanCopy, false, mean, true);
        Gadgetron::scal((ValueType)M, MMH);
        Gadgetron::subtract(eigenVectors, MMH, eigenVectors);
        Gadgetron::scal((ValueType)1.0/(M-1), eigenVectors);

        //MMH.print(std::cout);
        //eigenVectors.print(std::cout);

        hoMatrix<T> EH(eigenVectors);
        conjugatetrans(eigenVectors, EH);
        Gadgetron::add(eigenVectors, EH, eigenVectors);
        Gadgetron::scal( (ValueType)(0.5), eigenVectors);

        //eigenVectors.print(std::cout);

        Gadgetron::heev(eigenVectors, eigenValues);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::KLT_eigenAnalysis(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
KLT_applyEigen(const hoMatrix<T>& data, hoMatrix<T>& dataEigen, const hoMatrix<T>& eigenVectors)
{
    try
    {
        size_t M = data.rows();
        size_t N = data.cols();

        GADGET_CHECK_RETURN_FALSE(eigenVectors.rows()==N);

        size_t K = eigenVectors.cols();

        GADGET_CHECK_RETURN_FALSE(dataEigen.createMatrix(M, K));
        Gadgetron::clear(dataEigen);

        // M*N multiplies N*K
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::gemm(dataEigen, data, false, eigenVectors, false));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::KLT_applyEigen(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
KLT_applyEigen(const hoNDArray<T>& data, hoNDArray<T>& dataEigen, const hoMatrix<T>& eigenVectors)
{
    try
    {
        size_t M = data.get_size(0);
        size_t N = data.get_size(1);

        GADGET_CHECK_RETURN_FALSE(eigenVectors.rows()==N);

        size_t K = eigenVectors.cols();

        dataEigen.create(M, K);

        hoNDArray<T> eigenVec(eigenVectors.get_dimensions(), const_cast<T*>(eigenVectors.begin()));

        // M*N multiplies N*K
        Gadgetron::clear(dataEigen);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::gemm(dataEigen, data, false, eigenVec, false));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::KLT_applyEigen(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
KLT_numberOfKeptModes(const hoMatrix<T>& eigenValues, double thres, long long& numOfModesKept)
{
    try
    {
        size_t M = eigenValues.rows();

        if ( thres <= 0 )
        {
            numOfModesKept = (long long)M;
            return true;
        }

        long long m;
        for ( m=M-2; m>=0; m-- )
        {
            if ( std::abs(eigenValues(m,0)) < thres*std::abs(eigenValues(M-1,0)) )
            {
                break;
            }
        }

        numOfModesKept = M - m -1;

        if ( numOfModesKept <= 0 )
        {
            GWARN_STREAM("KLT_numberOfKeptModes(...) - numOfModesKept <= 0 : " << thres);
            GWARN_STREAM("KLT_numberOfKeptModes(...) - keep all modes : " << M);
            numOfModesKept = (long long)M;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::KLT_numberOfKeptModes(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
pruneEigenVectorMatrix(const hoMatrix<T>& eigenVectors, long long numOfModesKept, hoMatrix<T>& eigenVectorsPruned)
{
    try
    {
        size_t M = eigenVectors.rows();
        size_t N = eigenVectors.cols();

        if ( numOfModesKept<=0 || numOfModesKept>(long long)N )
        {
            GWARN_STREAM("gtPlusISMRMRDReconUtil<T>::pruneEigenVectorMatrix(...) - numOfModesKept<=0 || numOfModesKept>N : " << numOfModesKept);
            eigenVectorsPruned = eigenVectors;
            return true;
        }

        GADGET_CHECK_RETURN_FALSE(eigenVectorsPruned.createMatrix(M, numOfModesKept));
        GADGET_CHECK_RETURN_FALSE(eigenVectors.subMatrix(eigenVectorsPruned, 0, M-1, N-numOfModesKept, N-1));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::pruneEigenVectorMatrix(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
computeKLTCoeff(const hoNDArray<T>& data, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, bool isChaLastDim)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        hoMatrix<T> eigenVectors;
        hoMatrix<T> A;

        if ( isChaLastDim )
        {
            size_t CHA = data.get_size(NDim-1);
            size_t N = data.get_number_of_elements()/CHA;

            GADGET_CHECK_RETURN_FALSE(A.createMatrix(N, CHA, const_cast<T*>(data.begin())));
            GADGET_CHECK_RETURN_FALSE(KLT_eigenAnalysis(A, coeff, eigenValues));
        }
        else
        {
            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t CHA = data.get_size(2);

            if ( NDim == 3 )
            {
                GADGET_CHECK_RETURN_FALSE(A.createMatrix(RO*E1, CHA, const_cast<T*>(data.begin())));
                GADGET_CHECK_RETURN_FALSE(KLT_eigenAnalysis(A, coeff, eigenValues));
            }
            else if ( NDim == 4 )
            {
                size_t N = data.get_size(3);
                hoNDArray<T> dataP(RO, E1, N, CHA);

                std::vector<size_t> dimOrder(4);
                dimOrder[0] = 0;
                dimOrder[1] = 1;
                dimOrder[2] = 3;
                dimOrder[3] = 2;

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(const_cast< hoNDArray<T>* >(&data), &dataP, &dimOrder));

                GADGET_CHECK_RETURN_FALSE(A.createMatrix(RO*E1*N, CHA, dataP.begin()));

                GADGET_CHECK_RETURN_FALSE(KLT_eigenAnalysis(A, coeff, eigenValues));
            }
            else if ( NDim >= 5 )
            {
                std::vector<size_t> dimOrder(NDim);
                size_t l;
                for ( l=0; l<NDim; l++ )
                {
                    dimOrder[l] = l;
                }
                dimOrder[2] = NDim-1;
                dimOrder[NDim-1] = 2;

                hoNDArray<T> dataP(data);
                permute(&dataP, &dimOrder);

                size_t num = data.get_number_of_elements()/CHA;
                GADGET_CHECK_RETURN_FALSE(A.createMatrix(num, CHA, dataP.begin()));

                GADGET_CHECK_RETURN_FALSE(KLT_eigenAnalysis(A, coeff, eigenValues));
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::computeKLTCoeff(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
computeKLCoilCompressionCoeff(const hoNDArray<T>& data, double thres, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, bool isChaLastDim)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        hoMatrix<T> eigenVectors;
        GADGET_CHECK_RETURN_FALSE(computeKLTCoeff(data, eigenVectors, eigenValues, isChaLastDim));

        long long numOfModesKept;
        GADGET_CHECK_RETURN_FALSE(KLT_numberOfKeptModes(eigenValues, thres, numOfModesKept));
        GADGET_CHECK_RETURN_FALSE(pruneEigenVectorMatrix(eigenVectors, numOfModesKept, coeff));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::computeKLCoilCompressionCoeff(thres) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
computeKLCoilCompressionCoeff(const hoNDArray<T>& data, int numOfModesKept, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, bool isChaLastDim)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        hoMatrix<T> eigenVectors;
        GADGET_CHECK_RETURN_FALSE(computeKLTCoeff(data, eigenVectors, eigenValues, isChaLastDim));
        GADGET_CHECK_RETURN_FALSE(pruneEigenVectorMatrix(eigenVectors, numOfModesKept, coeff));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::computeKLCoilCompressionCoeff(numOfModesKept) ... ");
        return false;
    }

    return true;
}

template <typename T> 
inline bool gtPlusISMRMRDReconUtil<T>::
computeKLCoilCompression(const hoNDArray<T>& data, double thres, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, hoNDArray<T>& dataEigen, bool isChaLastDim)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(computeKLCoilCompressionCoeff(data, thres, coeff, eigenValues, isChaLastDim));
        GADGET_CHECK_RETURN_FALSE(appyKLCoilCompressionCoeff(data, coeff, dataEigen, isChaLastDim));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::computeKLCoilCompression(thres) ... ");
        return false;
    }

    return true;
}

template <typename T> 
inline bool gtPlusISMRMRDReconUtil<T>::
computeKLCoilCompression(const hoNDArray<T>& data, int numOfModesKept, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, hoNDArray<T>& dataEigen, bool isChaLastDim)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(computeKLCoilCompressionCoeff(data, numOfModesKept, coeff, eigenValues, isChaLastDim));
        GADGET_CHECK_RETURN_FALSE(appyKLCoilCompressionCoeff(data, coeff, dataEigen, isChaLastDim));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::computeKLCoilCompression(numOfModesKept) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
appyKLCoilCompressionCoeff(const hoNDArray<T>& data, const hoMatrix<T>& coeff, hoNDArray<T>& dataEigen, bool isChaLastDim)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        boost::shared_ptr< std::vector<size_t> > dim = data.get_dimensions();

        size_t dstCHA = coeff.cols();

        // D = A * V
        hoMatrix<T> A;
        hoMatrix<T> D;

        if ( isChaLastDim )
        {
            size_t CHA = data.get_size(NDim-1);
            size_t N = data.get_number_of_elements()/CHA;

            hoNDArray<T> A_tmp(N, CHA, const_cast<T*>(data.begin()));
            // GADGET_CHECK_RETURN_FALSE(A.createMatrix(CHA, N, const_cast<T*>(data.begin())));

            std::vector<size_t> dimEigen(*dim);
            dimEigen[NDim-1] = dstCHA;
            dataEigen.create(&dimEigen);

            hoNDArray<T> D_tmp(N, dstCHA, dataEigen.begin());
            // GADGET_CHECK_RETURN_FALSE(D.createMatrix(dstCHA, N, dataEigen.begin()));

            GADGET_CHECK_RETURN_FALSE(KLT_applyEigen(A_tmp, D_tmp, coeff));
        }
        else
        {
            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t CHA = data.get_size(2);

            if ( NDim == 3 )
            {
                GADGET_CHECK_RETURN_FALSE(A.createMatrix(RO*E1, CHA, const_cast<T*>(data.begin())));

                dataEigen.create(RO, E1, dstCHA);
                GADGET_CHECK_RETURN_FALSE(D.createMatrix(RO*E1, dstCHA, dataEigen.begin()));

                GADGET_CHECK_RETURN_FALSE(KLT_applyEigen(A, D, coeff));
            }
            else if ( NDim == 4 )
            {
                size_t N = data.get_size(3);
                hoNDArray<T> dataP(RO, E1, N, CHA);

                std::vector<size_t> dimOrder(4);
                dimOrder[0] = 0;
                dimOrder[1] = 1;
                dimOrder[2] = 3;
                dimOrder[3] = 2;

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(const_cast< hoNDArray<T>* >(&data), &dataP, &dimOrder));

                GADGET_CHECK_RETURN_FALSE(A.createMatrix(RO*E1*N, CHA, dataP.begin()));

                hoNDArray<T> dataEigenP(RO, E1, N, dstCHA);
                GADGET_CHECK_RETURN_FALSE(D.createMatrix(RO*E1*N, dstCHA, dataEigenP.begin()));

                GADGET_CHECK_RETURN_FALSE(KLT_applyEigen(A, D, coeff));

                dataEigen.create(RO, E1, dstCHA, N);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&dataEigenP, &dataEigen, &dimOrder));
            }
            else if ( NDim >= 5 )
            {
                std::vector<size_t> dimOrder(NDim);
                size_t l;
                for ( l=0; l<NDim; l++ )
                {
                    dimOrder[l] = l;
                }
                dimOrder[2] = NDim-1;
                dimOrder[NDim-1] = 2;

                boost::shared_ptr< hoNDArray<T> > dataP = permute(const_cast< hoNDArray<T>* >(&data), &dimOrder);

                size_t num = data.get_number_of_elements()/CHA;
                GADGET_CHECK_RETURN_FALSE(A.createMatrix(num, CHA, dataP->begin()));

                boost::shared_ptr< std::vector<size_t> > dimP = dataP->get_dimensions();
                (*dimP)[NDim-1] = dstCHA;

                dataEigen.create(dimP);
                GADGET_CHECK_RETURN_FALSE(D.createMatrix(num, dstCHA, dataEigen.begin()));

                GADGET_CHECK_RETURN_FALSE(KLT_applyEigen(A, D, coeff));

                dataP = permute(&dataEigen, &dimOrder);
                dataEigen =  *dataP;
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::appyKLCoilCompressionCoeff(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
applyKLCoilCompressionCoeff(const hoNDArray<T>& data, const std::vector<hoMatrix<T> >& coeff, hoNDArray<T>& dataEigen, bool isChaLastDim)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        GADGET_CHECK_RETURN_FALSE(coeff.size()>=data.get_size(NDim-1));

        size_t LastDim = coeff.size();
        size_t dstCHA = coeff[0].cols();

        size_t n;
        for ( n=1; n<LastDim; n++ )
        {
            GADGET_CHECK_RETURN_FALSE(coeff[n].cols()==dstCHA);
        }

        size_t LastDimData = data.get_size(NDim-1);
        boost::shared_ptr< std::vector<size_t> > dim = data.get_dimensions();
        long long N = data.get_number_of_elements()/LastDimData;

        std::vector<size_t> dimEigen(*dim);

        if ( isChaLastDim )
        {
            dimEigen[NDim-2] = dstCHA;
        }
        else
        {
            dimEigen[2] = dstCHA;
        }

        dataEigen.create(&dimEigen);
        long long eigenN = dataEigen.get_number_of_elements()/LastDimData;

        std::vector<size_t> dimLastDim(NDim-1);
        for ( n=0; n<NDim-1; n++ )
        {
            dimLastDim[n] = (*dim)[n];
        }

        hoNDArray<T> dataEigenLastDim;
        for ( n=0; n<LastDimData; n++ )
        {
            hoNDArray<T> dataLastDim(&dimLastDim, const_cast<T*>(data.begin()+n*N));
            GADGET_CHECK_RETURN_FALSE(appyKLCoilCompressionCoeff(dataLastDim, coeff[n], dataEigenLastDim, isChaLastDim));
            memcpy(dataEigen.begin()+n*eigenN, dataEigenLastDim.begin(), dataEigenLastDim.get_number_of_bytes());
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::applyKLCoilCompressionCoeff(std::vector<hoMatrix<T> >& coeff) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::computeKLFilter(const hoNDArray<T>& data, size_t numOfModesKept, hoNDArray<T>& dataKLF)
{
    try
    {
        if ( !dataKLF.dimensions_equal(&data) )
        {
            dataKLF = data;
        }

        size_t NDim = data.get_number_of_dimensions();
        size_t M = data.get_size(NDim-1);
        size_t N = data.get_number_of_elements()/M;

        if ( numOfModesKept > M ) numOfModesKept = M;

        hoMatrix<T> A(N, M, const_cast<T*>(data.begin()));

        hoMatrix<T> eigenVectors, eigenValues;
        GADGET_CHECK_RETURN_FALSE(KLT_eigenAnalysis(A, eigenVectors, eigenValues));

        hoMatrix<T> E(eigenVectors);
        size_t r, c;
        for ( c=0; c<M-numOfModesKept+1; c++ )
        {
            for ( r=0; r<M; r++ )
            {
                E(r, c) = T(0);
            }
        }

        hoMatrix<T> ET(eigenVectors);
        GADGET_CHECK_RETURN_FALSE(Gadgetron::conjugatetrans(eigenVectors, ET));

        hoMatrix<T> EET(M, M);
        Gadgetron::clear(EET);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::gemm(EET, E, false, ET, false));

        hoMatrix<T> R(N, M, dataKLF.begin());
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::gemm(R, A, false, EET, false));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::computeKLFilter(...) ... ");
        return false;
    }

    return true;
}

// ------------------------------------------------------------------------
// zero-padding resize
// ------------------------------------------------------------------------
template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
zpadRange(size_t srcSize, size_t dstSize, size_t& start, size_t& end)
{
    try
    {
        if ( srcSize >= dstSize )
        {
            start = 0;
            end = srcSize-1;
            return true;
        }

        //unsigned srcCenterInd = srcSize/2;
        //unsigned dstCenterInd = dstSize/2;

        start = (dstSize/2) - (srcSize/2);
        end = srcSize + start -1;

        //start = std::floor((double)dstSize/2.0)+1+std::ceil(-1.0 * (double)srcSize/2.0)-1;
        //end = std::floor((double)dstSize/2.0)+std::ceil(srcSize/2.0)-1;
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::zpadRange(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
zeropad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataPadded)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);

        if ( RO==sizeX && E1==sizeY )
        {
            dataPadded = data;
            return true;
        }

        size_t sRO, eRO, sE1, eE1;
        GADGET_CHECK_RETURN_FALSE(zpadRange(RO, sizeX, sRO, eRO));
        GADGET_CHECK_RETURN_FALSE(zpadRange(E1, sizeY, sE1, eE1));

        boost::shared_ptr< std::vector<size_t> > dimPadded = data.get_dimensions();
        (*dimPadded)[0] = sizeX;
        (*dimPadded)[1] = sizeY;
        dataPadded.create(dimPadded);
        Gadgetron::clear(&dataPadded);

        size_t num = data.get_number_of_elements()/(RO*E1);

        long long n;

        const T* pData = data.begin();
        T* pDataPadded = dataPadded.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE1, eE1, sRO, RO, E1, pData, pDataPadded, sizeX, sizeY)
        for ( n=0; n<(long long)num; n++ )
        {
            for ( size_t y=sE1; y<=eE1; y++ )
            {
                memcpy(pDataPadded+n*sizeX*sizeY+y*sizeX+sRO, pData+n*RO*E1+(y-sE1)*RO, sizeof(T)*RO);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::zeropad2D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
zeropad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);
        GADGET_CHECK_RETURN_FALSE(sizeZ>=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataPadded = data;
            return true;
        }

        size_t sRO, eRO, sE1, eE1, sE2, eE2;
        GADGET_CHECK_RETURN_FALSE(zpadRange(RO, sizeX, sRO, eRO));
        GADGET_CHECK_RETURN_FALSE(zpadRange(E1, sizeY, sE1, eE1));
        GADGET_CHECK_RETURN_FALSE(zpadRange(E2, sizeZ, sE2, eE2));

        boost::shared_ptr< std::vector<size_t> > dimPadded = data.get_dimensions();
        (*dimPadded)[0] = sizeX;
        (*dimPadded)[1] = sizeY;
        (*dimPadded)[2] = sizeZ;
        dataPadded.create(dimPadded);
        Gadgetron::clear(&dataPadded);

        size_t num = data.get_number_of_elements()/(RO*E1*E2);

        long long n;

        const T* pData = data.begin();
        T* pDataPadded = dataPadded.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, pData, pDataPadded, sizeX, sizeY, sizeZ)
        for ( n=0; n<(long long)num; n++ )
        {
            T* pDst = pDataPadded+n*sizeX*sizeY*sizeZ;
            T* pSrc = const_cast<T*>(pData)+n*RO*E1*E2;

            long long z;
            // #pragma omp parallel for default(none) private(z) shared(pDst, pSrc, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, sizeX, sizeY, sizeZ) num_threads(2)
            for ( z=(long long)sE2; z<=(long long)eE2; z++ )
            {
                long long o1 = z*sizeX*sizeY + sRO;
                long long o2 = (z-sE2)*RO*E1;
                for ( size_t y=sE1; y<=eE1; y++ )
                {
                    memcpy(pDst+o1+y*sizeX, pSrc+o2+(y-sE1)*RO, sizeof(T)*RO);
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::zeropad3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
zeropad3DNoPresetZeros(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t srcCHA = data.get_size(3);
        size_t dstCHA = data.get_size(4);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);
        GADGET_CHECK_RETURN_FALSE(sizeZ>=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataPadded = data;
            return true;
        }

        size_t sRO, eRO, sE1, eE1, sE2, eE2;
        GADGET_CHECK_RETURN_FALSE(zpadRange(RO, sizeX, sRO, eRO));
        GADGET_CHECK_RETURN_FALSE(zpadRange(E1, sizeY, sE1, eE1));
        GADGET_CHECK_RETURN_FALSE(zpadRange(E2, sizeZ, sE2, eE2));

        GADGET_CHECK_RETURN_FALSE(dataPadded.get_size(0)==sizeX);
        GADGET_CHECK_RETURN_FALSE(dataPadded.get_size(1)==sizeY);
        GADGET_CHECK_RETURN_FALSE(dataPadded.get_size(2)==sizeZ);
        GADGET_CHECK_RETURN_FALSE(dataPadded.get_size(3)==srcCHA);
        GADGET_CHECK_RETURN_FALSE(dataPadded.get_size(4)==dstCHA);

        size_t num = data.get_number_of_elements()/(RO*E1*E2);

        long long n;

        const T* pData = data.begin();
        T* pDataPadded = dataPadded.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, pData, pDataPadded, sizeX, sizeY, sizeZ)
        for ( n=0; n<(long long)num; n++ )
        {
            T* pDst = pDataPadded+n*sizeX*sizeY*sizeZ;
            T* pSrc = const_cast<T*>(pData)+n*RO*E1*E2;

            long long z;
            //#pragma omp parallel for default(none) private(z) shared(pDst, pSrc, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, sizeX, sizeY, sizeZ) num_threads(2)
            for ( z=(long long)sE2; z<=(long long)eE2; z++ )
            {
                long long o1 = z*sizeX*sizeY + sRO;
                long long o2 = (z-sE2)*RO*E1;
                for ( size_t y=sE1; y<=eE1; y++ )
                {
                    memcpy(pDst+o1+y*sizeX, pSrc+o2+(y-sE1)*RO, sizeof(T)*RO);
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::zeropad3DNoPresetZeros(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
cutpad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataCut)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);

        GADGET_CHECK_RETURN_FALSE(sizeX<=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY<=E1);

        if ( RO==sizeX && E1==sizeY )
        {
            dataCut = data;
            return true;
        }

        size_t sRO, eRO, sE1, eE1;
        GADGET_CHECK_RETURN_FALSE(zpadRange(sizeX, RO, sRO, eRO));
        GADGET_CHECK_RETURN_FALSE(zpadRange(sizeY, E1, sE1, eE1));

        boost::shared_ptr< std::vector<size_t> > dim = data.get_dimensions();
        (*dim)[0] = sizeX;
        (*dim)[1] = sizeY;
        dataCut.create(dim);

        size_t num = data.get_number_of_elements()/(RO*E1);

        long long n;

        const T* pData = data.begin();
        T* pDataCut = dataCut.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE1, eE1, sRO, RO, E1, pData, pDataCut, sizeX, sizeY)
        for ( n=0; n<(long long)num; n++ )
        {
            for ( size_t y=sE1; y<=eE1; y++ )
            {
                memcpy(pDataCut+n*sizeX*sizeY+(y-sE1)*sizeX, pData+n*RO*E1+y*RO+sRO, sizeof(T)*sizeX);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::cutpad2D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
cutpad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataCut)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);

        GADGET_CHECK_RETURN_FALSE(sizeX<=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY<=E1);
        GADGET_CHECK_RETURN_FALSE(sizeZ<=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataCut = data;
            return true;
        }

        size_t sRO, eRO, sE1, eE1, sE2, eE2;
        GADGET_CHECK_RETURN_FALSE(zpadRange(sizeX, RO, sRO, eRO));
        GADGET_CHECK_RETURN_FALSE(zpadRange(sizeY, E1, sE1, eE1));
        GADGET_CHECK_RETURN_FALSE(zpadRange(sizeZ, E2, sE2, eE2));

        boost::shared_ptr< std::vector<size_t> > dim = data.get_dimensions();
        (*dim)[0] = sizeX;
        (*dim)[1] = sizeY;
        (*dim)[2] = sizeZ;
        dataCut.create(dim);

        size_t num = data.get_number_of_elements()/(RO*E1*E2);

        long long n;

        const T* pData = data.begin();
        T* pDataCut = dataCut.begin();

        #pragma omp parallel for default(none) private(n) shared(num, sE2, eE2, sE1, eE1, sRO, RO, E1, E2, pData, pDataCut, sizeX, sizeY, sizeZ)
        for ( n=0; n<(long long)num; n++ )
        {
            for ( size_t z=sE2; z<=eE2; z++ )
            {
                for ( size_t y=sE1; y<=eE1; y++ )
                {
                    memcpy(pDataCut+n*sizeX*sizeY*sizeZ+(z-sE2)*sizeX*sizeY+(y-sE1)*sizeX, pData+n*RO*E1*E2+z*RO*E1+y*RO+sRO, sizeof(T)*sizeX);
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::cutpad3D(...) ... ");
        return false;
    }

    return true;
}

// ------------------------------------------------------------------------
// kspace filter
// ------------------------------------------------------------------------
template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
compute2DFilterFromTwo1D(const hoNDArray<T>& fx, const hoNDArray<T>& fy, hoNDArray<T>& fxy)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);

        fxy.create(RO, E1);
        T* pFxy = fxy.begin();

        size_t x, y;

        for ( y=0; y<E1; y++ )
        {
            for ( x=0; x<RO; x++ )
            {
                pFxy[y*RO+x] = fx(x) * fy(y);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::compute2DFilterFromTwo1D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
compute2DFilterFromTwo1D(const hoNDArray<float>& fx, const hoNDArray<float>& fy, hoNDArray< std::complex<float> >& fxy)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);

        fxy.create(RO, E1);
         std::complex<float> * pFxy = fxy.begin();

        size_t x, y;

        for ( y=0; y<E1; y++ )
        {
            for ( x=0; x<RO; x++ )
            {
                pFxy[y*RO+x] =  std::complex<float> (fx(x) * fy(y));
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::compute2DFilterFromTwo1D(float) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
compute2DFilterFromTwo1D(const hoNDArray<double>& fx, const hoNDArray<double>& fy, hoNDArray< std::complex<double> >& fxy)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);

        fxy.create(RO, E1);
         std::complex<double> * pFxy = fxy.begin();

        size_t x, y;

        for ( y=0; y<E1; y++ )
        {
            for ( x=0; x<RO; x++ )
            {
                pFxy[y*RO+x] =  std::complex<double> (fx(x) * fy(y));
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::compute2DFilterFromTwo1D(double) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
compute3DFilterFromThree1D(const hoNDArray<T>& fx, const hoNDArray<T>& fy, const hoNDArray<T>& fz, hoNDArray<T>& fxyz)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);
        size_t E2 = fz.get_size(0);

        fxyz.create(RO, E1, E2);
        T* pFxyz = fxyz.begin();

        const T* px = fx.begin();
        const T* py = fy.begin();
        const T* pz = fz.begin();

        size_t x, y, z;

        T vz, vy, vx;

        size_t ind = 0;
        for ( z=0; z<E2; z++ )
        {
            vz = pz[z];
            for ( y=0; y<E1; y++ )
            {
                vy = py[y];
                for ( x=0; x<RO; x++ )
                {
                    vx = px[x];
                    pFxyz[ind] = (vx*vz*vy);
                    ind++;
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::compute3DFilterFromThree1D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
compute3DFilterFromThree1D(const hoNDArray<float>& fx, const hoNDArray<float>& fy, const hoNDArray<float>& fz, hoNDArray< std::complex<float> >& fxyz)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);
        size_t E2 = fz.get_size(0);

        fxyz.create(RO, E1, E2);
         std::complex<float> * pFxyz = fxyz.begin();

        size_t x, y, z;

        for ( z=0; z<E2; z++ )
        {
            for ( y=0; y<E1; y++ )
            {
                for ( x=0; x<RO; x++ )
                {
                    pFxyz[z+RO*E1+y*RO+x] =  std::complex<float> (fx(x)*fy(y)*fz(z));
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::compute3DFilterFromThree1D(float) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
compute3DFilterFromThree1D(const hoNDArray<double>& fx, const hoNDArray<double>& fy, const hoNDArray<double>& fz, hoNDArray< std::complex<double> >& fxyz)
{
    try
    {
        size_t RO = fx.get_size(0);
        size_t E1 = fy.get_size(0);
        size_t E2 = fz.get_size(0);

        fxyz.create(RO, E1, E2);
         std::complex<double> * pFxyz = fxyz.begin();

        size_t x, y, z;

        for ( z=0; z<E2; z++ )
        {
            for ( y=0; y<E1; y++ )
            {
                for ( x=0; x<RO; x++ )
                {
                    pFxyz[z+RO*E1+y*RO+x] =  std::complex<double> (fx(x)*fy(y)*fz(z));
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::compute3DFilterFromThree1D(double) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterRO(hoNDArray<T>& data, const hoNDArray<T>& fRO)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, fRO, data));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterRO(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterRO(const hoNDArray<T>& data, const hoNDArray<T>& fRO, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, fRO, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterRO(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterROE1(const hoNDArray<T>& data, const hoNDArray<T>& fROE1, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)*data.get_size(1)==fROE1.get_number_of_elements());
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, fROE1, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterROE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterROE1(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());

        hoNDArray<T> fxy;
        compute2DFilterFromTwo1D(fRO, fE1, fxy);

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, fxy, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterROE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterE1(const hoNDArray<T>& data, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());

        hoNDArray<T> fRO(data.get_size(0));
        fRO.fill(T(1.0));

        hoNDArray<T> fxy;
        compute2DFilterFromTwo1D(fRO, fE1, fxy);

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, fxy, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterE2(const hoNDArray<T>& data, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(4)==fE2.get_number_of_elements());

        hoNDArray<T> fRO(data.get_size(0));
        fRO.fill(T(1.0));

        hoNDArray<T> fE1(data.get_size(1));
        fE1.fill(T(1.0));

        hoNDArray<T> fxyz;
        compute3DFilterFromThree1D(fRO, fE1, fE2, fxyz);

        GADGET_CHECK_RETURN_FALSE(kspacefilterROE1E2(data, fxyz, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterE2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterROE2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(4)==fE2.get_number_of_elements());

        hoNDArray<T> fE1(data.get_size(1));
        fE1.fill(T(1.0));

        hoNDArray<T> fxyz;
        compute3DFilterFromThree1D(fRO, fE1, fE2, fxyz);

        GADGET_CHECK_RETURN_FALSE(kspacefilterROE1E2(data, fxyz, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterROE2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(4)==fE2.get_number_of_elements());

        hoNDArray<T> fRO(data.get_size(0));
        fRO.fill(T(1.0));

        hoNDArray<T> fxyz;
        compute3DFilterFromThree1D(fRO, fE1, fE2, fxyz);

        GADGET_CHECK_RETURN_FALSE(kspacefilterROE1E2(data, fxyz, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fROE1E2, hoNDArray<T>& dataFiltered)
{
    try
    {
        if ( data.get_size(2)==1 && data.get_size(3)==1 )
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, fROE1E2, dataFiltered));
        }
        else
        {
            size_t NDim = data.get_number_of_dimensions();
            std::vector<size_t> order(data.get_number_of_dimensions(), 1);

            size_t ii;
            for ( ii=0; ii<NDim; ii++ )
            {
                order[ii] = ii;
            }

            order[0] = 0;
            order[1] = 1;
            order[2] = 4;
            order[3] = 2;
            order[4] = 3;

            boost::shared_ptr< hoNDArray<T> > data_permuted = Gadgetron::permute(const_cast<hoNDArray<T>*>(&data), &order);

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, fROE1E2, dataFiltered));

            order[0] = 0;
            order[1] = 1;
            order[2] = 3;
            order[3] = 4;
            order[4] = 2;

            data_permuted = Gadgetron::permute(&dataFiltered, &order);
            dataFiltered = *data_permuted;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterROE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspacefilterROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(4)==fE2.get_number_of_elements());

        hoNDArray<T> fxyz;
        compute3DFilterFromThree1D(fRO, fE1, fE2, fxyz);

        GADGET_CHECK_RETURN_FALSE(kspacefilterROE1E2(data, fxyz, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspacefilterROE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspace3DfilterE2(const hoNDArray<T>& data, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        hoNDArray<T> fRO(data.get_size(0));
        fRO.fill(T(1.0));

        hoNDArray<T> fE1(data.get_size(1));
        fE1.fill(T(1.0));

        hoNDArray<T> fxyz;
        compute3DFilterFromThree1D(fRO, fE1, fE2, fxyz);

        GADGET_CHECK_RETURN_FALSE(kspace3DfilterROE1E2(data, fxyz, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspace3DfilterE2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspace3DfilterROE2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        hoNDArray<T> fE1(data.get_size(1));
        fE1.fill(T(1.0));

        hoNDArray<T> fxyz;
        compute3DFilterFromThree1D(fRO, fE1, fE2, fxyz);

        GADGET_CHECK_RETURN_FALSE(kspace3DfilterROE1E2(data, fxyz, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspace3DfilterROE2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspace3DfilterE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        hoNDArray<T> fRO(data.get_size(0));
        fRO.fill(T(1.0));

        hoNDArray<T> fxyz;
        compute3DFilterFromThree1D(fRO, fE1, fE2, fxyz);

        GADGET_CHECK_RETURN_FALSE(kspace3DfilterROE1E2(data, fxyz, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspace3DfilterE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspace3DfilterROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fROE1E2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fROE1E2.get_size(0));
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fROE1E2.get_size(1));
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fROE1E2.get_size(2));

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, fROE1E2, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspace3DfilterROE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
kspace3DfilterROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        hoNDArray<T> fxyz(fRO.get_number_of_elements(), fE1.get_number_of_elements(), fE2.get_number_of_elements());
        compute3DFilterFromThree1D(fRO, fE1, fE2, fxyz);

        GADGET_CHECK_RETURN_FALSE(kspace3DfilterROE1E2(data, fxyz, dataFiltered));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::kspace3DfilterROE1E2(...) ... ");
        return false;
    }

    return true;
}

// ------------------------------------------------------------------------
// compute kspace filters
// ------------------------------------------------------------------------

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
generateSymmetricFilter(size_t len, size_t start, size_t end, hoNDArray<T>& filter, ISMRMRDKSPACEFILTER filterType, double sigma, size_t width)
{
    try
    {
        if ( len == 0 ) return true;

        if ( start > len-1 ) start = 0;
        if ( end > len-1 ) end = len-1;

        if ( start > end )
        {
            start = 0;
            end = len-1;
        }

        filter.create(len);
        Gadgetron::fill(filter, T(1.0));

        if ( width==0 || width>=len ) width = 1;

        size_t ii;
        switch (filterType)
        {
            case ISMRMRD_FILTER_GAUSSIAN:
                {
                    double r = -1.0*sigma*sigma/2;

                    if ( len%2 == 0 )
                    {
                        // to make sure the zero points match and boundary of filters are symmetric
                        double stepSize = 2.0/(len-2);
                        std::vector<double> x(len-1);

                        for ( ii=0; ii<len-1; ii++ )
                        {
                            x[ii] = -1 + ii*stepSize;
                        }

                        for ( ii=0; ii<len-1; ii++ )
                        {
                            filter(ii+1) = T( (value_type)(std::exp(r*(x[ii]*x[ii]))) );
                        }

                        filter(0) = T(0);
                    }
                    else
                    {
                        double stepSize = 2.0/(len-1);
                        std::vector<double> x(len);

                        for ( ii=0; ii<len; ii++ )
                        {
                            x[ii] = -1 + ii*stepSize;
                        }

                        for ( ii=0; ii<len; ii++ )
                        {
                            filter(ii) = T( (value_type)(std::exp(r*(x[ii]*x[ii]))) );
                        }
                    }
                }
            break;

            case ISMRMRD_FILTER_TAPERED_HANNING:
                 {
                    hoNDArray<T> w(width);

                    for ( ii=1; ii<=width; ii++ )
                    {
                        w(ii-1) = T( (value_type)(0.5 * ( 1 - std::cos( 2.0*M_PI*ii/(2*width+1) ) )) );
                    }

                    if ( len%2 == 0 )
                    {
                        for ( ii=1; ii<=width; ii++ )
                        {
                            filter(ii) = w(ii-1);
                            filter(len-ii) = filter(ii);
                        }

                        filter(0) = T(0);
                    }
                    else
                    {
                        for ( ii=1; ii<=width; ii++ )
                        {
                            filter(ii-1) = w(ii-1);
                            filter(len-ii) = filter(ii-1);
                        }
                    }
                }
            break;

            // symmetric hanning
            //  does not include the first and last zero sample
            case ISMRMRD_FILTER_HANNING:
                 {
                    if ( len%2 == 0 )
                    {
                        size_t N = len-1;
                        double halfLen = (double)( (N+1)/2 );
                        for ( ii=1; ii<=halfLen; ii++ )
                        {
                            filter(ii) = T( (value_type)(0.5 * ( 1 - std::cos( 2.0*M_PI*ii/(N+1) ) )) );
                        }

                        for ( ii=(size_t)halfLen; ii<N; ii++ )
                        {
                            filter(ii+1) = filter(N-ii);
                        }

                        filter(0) = T(0);
                    }
                    else
                    {
                        double halfLen = (double)( (len+1)/2 );
                        for ( ii=1; ii<=(size_t)halfLen; ii++ )
                        {
                            filter(ii-1) = T( (value_type)(0.5 * ( 1 - std::cos( 2.0*M_PI*ii/(len+1) ) )) );
                        }

                        for ( ii=(size_t)halfLen; ii<len; ii++ )
                        {
                            filter(ii) = filter(len-1-ii);
                        }
                    }
                }
            break;

            default:
            break;
        }

        T sos = 0.0f;
        for ( ii=0; ii<len; ii++ )
        {
            sos += filter(ii)*filter(ii);
        }

        T r = (value_type)( 1.0/std::sqrt( std::abs(sos)/(len) ) );
        for ( ii=0; ii<len; ii++ )
        {
            filter(ii) *= r;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::generateSymmetricFilter(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
generateAsymmetricFilter(size_t len, size_t start, size_t end, hoNDArray<T>& filter, ISMRMRDKSPACEFILTER filterType, size_t width, bool densityComp)
{
    try
    {
        if ( len == 0 ) return true;

        if ( start > len-1 ) start = 0;
        if ( end > len-1 ) end = len-1;

        if ( start > end )
        {
            start = 0;
            end = len-1;
        }

        filter.create(len);
        Gadgetron::clear(filter);

        size_t ii;
        for ( ii=start; ii<=end; ii++ )
        {
            filter(ii) = T(1.0);
        }

        if ( width==0 || width>=len ) width = 1;

        hoNDArray<T> w(width);

        switch (filterType)
        {
            case ISMRMRD_FILTER_TAPERED_HANNING:
                 {
                    for ( ii=1; ii<=width; ii++ )
                    {
                        w(ii-1) = T( (value_type)(0.5 * ( 1 - std::cos( 2.0*M_PI*ii/(2*width+1) ) )) );
                    }
                }
            break;

            default:
                Gadgetron::fill(w, T(1.0));
            break;
        }

        if ( densityComp )
        {
            size_t startSym(0), endSym(len-1);
            GADGET_CHECK_RETURN_FALSE(findSymmetricSampledRegion(start, end, len/2, startSym, endSym));

            if ( start==0 && end==len-1 )
            {
                for ( ii=1; ii<=width; ii++ )
                {
                    filter(ii-1) = w(ii-1);
                    filter(len-ii) = filter(ii-1);
                }
            }

            if ( start==0 && end<len-1 )
            {
                for ( ii=0; ii<startSym; ii++ )
                {
                    filter(ii) = 2.0;
                }

                for ( ii=1; ii<=width; ii++ )
                {
                    filter(ii-1+startSym) = T(1.0) + w(width-ii);
                    filter(end-ii+1) = w(ii-1);
                }
            }

            if ( start>0 && end==len-1 )
            {
                for ( ii=endSym+1; ii<len; ii++ )
                {
                    filter(ii) = 2.0;
                }

                for ( ii=1; ii<=width; ii++ )
                {
                    filter(endSym-ii+1) = T(1.0) + w(width-ii);
                    filter(start+ii-1) = w(ii-1);
                }
            }

            if ( start>0 && end<len-1 )
            {
                if ( start==startSym && end==endSym )
                {
                    for ( ii=1; ii<=width; ii++ )
                    {
                        filter(start+ii-1) = w(ii-1);
                        filter(end-ii+1) = w(ii-1);
                    }
                }
                else if ( start==startSym && end>endSym )
                {
                    for ( ii=endSym+1; ii<=end; ii++ )
                    {
                        filter(ii) = 2.0;
                    }

                    for ( ii=1; ii<=width; ii++ )
                    {
                        filter(end-ii+1) = T(1.0) + w(ii-1);
                        filter(endSym-ii+1) = w(width-ii);
                        filter(start+ii-1) = w(ii-1);
                    }
                }
                else if ( start<startSym && end==endSym )
                {
                    for ( ii=start; ii<startSym; ii++ )
                    {
                        filter(ii) = 2.0;
                    }

                    for ( ii=1; ii<=width; ii++ )
                    {
                        filter(ii-1+start) = T(1.0) + w(ii-1);
                        filter(ii-1+startSym) = w(width-ii);
                        filter(end-ii+1) = w(ii-1);
                    }
                }
                else
                {
                    for ( ii=1; ii<=width; ii++ )
                    {
                        filter(start+ii-1) = w(ii-1);
                        filter(end-ii+1) = w(ii-1);
                    }
                }
            }
        }
        else
        {
            if ( start==0 && end==len-1 )
            {
                for ( ii=1; ii<=width; ii++ )
                {
                    filter(ii-1) = w(ii-1);
                    filter(len-ii) = filter(ii-1);
                }
            }

            if ( start==0 && end<len-1 )
            {
                for ( ii=1; ii<=width; ii++ )
                {
                    filter(end-ii+1) = w(ii-1);
                }
            }

            if ( start>0 && end==len-1 )
            {
                for ( ii=1; ii<=width; ii++ )
                {
                    filter(start+ii-1) = w(ii-1);
                }
            }

            if ( start>0 && end<len-1 )
            {
                for ( ii=1; ii<=width; ii++ )
                {
                    filter(start+ii-1) = w(ii-1);
                    filter(end-ii+1) = w(ii-1);
                }
            }
        }

        T sos = 0.0f;
        for ( ii=0; ii<len; ii++ )
        {
            sos += filter(ii)*filter(ii);
        }

        // T r = 1.0/std::sqrt( std::abs(sos)/len );
        T r = (value_type)( 1.0/std::sqrt( std::abs(sos)/(end-start+1) ) ); // SNR unit filter
        for ( ii=0; ii<len; ii++ )
        {
            filter(ii) *= r;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::generateAsymmetricFilter(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
generateSymmetricFilterForRef(size_t len, size_t start, size_t end, 
        hoNDArray<T>& filter, ISMRMRDKSPACEFILTER filterType, double sigma, size_t width)
{
    try
    {
        if ( len < 2 ) return true;

        GADGET_CHECK_RETURN_FALSE(start>=0&&end<=len-1&&start<=end);

        if ( start==0 && end==len-1 )
        {
            GADGET_CHECK_RETURN_FALSE(generateSymmetricFilter(len, 0, len-1, filter, filterType, sigma, width));
            return true;
        }

        size_t centerInd = len/2;

        size_t lenFilter(0); // make a symmetric filter with zero at the center
        size_t lenFilterEnd = 2*(end-centerInd)+1;
        size_t lenFilterStart = 2*(centerInd-start)+1;

        if ( start==0 && end<len-1 )
        {
            lenFilter = lenFilterEnd;
        }
        else if ( start>0 && end==len-1 )
        {
            lenFilter = lenFilterStart;
        }
        else if ( start>0 && end<len-1 )
        {
            lenFilter = ( (lenFilterStart<lenFilterEnd) ? lenFilterStart : lenFilterEnd );
        }
        else
        {
            GERROR_STREAM("Invalid inputs : start - end - len : " << start << " " << end << " " << len);
        }

        GADGET_CHECK_RETURN_FALSE(lenFilter>0);

        hoNDArray<T> filterSym(lenFilter);
        GADGET_CHECK_RETURN_FALSE(generateSymmetricFilter(lenFilter, 0, lenFilter-1, filterSym, filterType, sigma, width));

        filter.create(len);
        Gadgetron::clear(&filter);

        if ( start==0 && end<len-1 )
        {
            memcpy(filter.begin()+end-lenFilter+1, filterSym.begin(), filterSym.get_number_of_bytes());
            return true;
        }
        else if ( start>0 && end==len-1 )
        {
            memcpy(filter.begin()+start, filterSym.begin(), filterSym.get_number_of_bytes());
            return true;
        }
        else if ( start>0 && end<len-1 )
        {
            if ( lenFilter == lenFilterStart ) 
            {
                memcpy(filter.begin()+start, filterSym.begin(), filterSym.get_number_of_bytes());
            }
            else
            {
                memcpy(filter.begin()+end-lenFilter+1, filterSym.begin(), filterSym.get_number_of_bytes());
            }

            return true;
        }
        else
        {
            GERROR_STREAM("Invalid inputs : start - end - len : " << start << " " << end << " " << len);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::generateSymmetricFilterForRef(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
findSymmetricSampledRegion(size_t start, size_t end, size_t center, size_t& startSym, size_t& endSym)
{
    GADGET_CHECK_RETURN_FALSE(end>=start);
    GADGET_CHECK_RETURN_FALSE(center>=start);
    GADGET_CHECK_RETURN_FALSE(end>=center);

    size_t halfSizeStart = center - start;
    size_t halfSizeEnd =  end - center;

    if ( halfSizeStart > halfSizeEnd )
    {
        startSym = center - halfSizeEnd;
        endSym = center + halfSizeEnd;
    }
    else
    {
        startSym = center - halfSizeStart;
        endSym = center + halfSizeStart;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::computeFilterSNRUnitScaleFactor(const hoNDArray<T>& filter, T& scalFactor)
{
    size_t ii, len;

    len = filter.get_number_of_elements();
    if ( len == 0 )
    {
        scalFactor = T(1.0);
        return true;
    }

    T sos(0.0);
    for ( ii=0; ii<len; ii++ )
    {
        sos += filter(ii)*filter(ii);
    }

    scalFactor = (value_type)(1.0/std::sqrt( std::abs(sos)/len ));

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledRegion2D(const hoNDArray<T>& data, size_t& startRO, size_t& endRO, size_t& startE1, size_t& endE1)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();

        hoNDArray<typename realType<T>::Type> mag(data.get_dimensions()), magSum, magSumE1, magSumRO;
        Gadgetron::abs(data, mag);

        if ( NDim > 2 )
        {
            size_t ii;
            std::vector<size_t> dim;
            for ( ii=0; ii<NDim-2; ii++ )
            {
                mag.get_dimensions(dim);

                // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(mag, magSum));
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(mag, magSum, mag.get_number_of_dimensions() - 1));

                std::vector<size_t> dimSum(dim.size()-1);
                memcpy(&dimSum[0], &dim[0], sizeof(size_t)*dimSum.size());
                magSum.reshape(dimSum);

                mag = magSum;
            }
        }

        size_t RO = mag.get_size(0);
        size_t E1 = mag.get_size(1);

        startRO = RO-1;
        endRO = 0;

        startE1 = E1-1;
        endE1 = 0;

        size_t ro, e1;

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, magSumE1, 1));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, magSumRO, 0));

        for ( ro=0; ro<RO; ro++ )
        {
            if ( magSumE1(ro) > 0 )
            {
                if ( ro < startRO ) startRO = ro;
                if ( ro > endRO ) endRO = ro;
            }
        }

        for ( e1=0; e1<E1; e1++ )
        {
            if ( magSumRO(e1) > 0 )
            {
                if ( e1 < startE1 ) startE1 = e1;
                if ( e1 > endE1 ) endE1 = e1;
            }
        }

        if ( startRO > endRO )
        {
            startRO = 0;
            endRO = RO-1;
        }

        if ( startE1 > endE1 )
        {
            startE1 = 0;
            endE1 = E1-1;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledRegion2D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledRegion3D(const hoNDArray<T>& data, size_t& startRO, size_t& endRO, size_t& startE1, size_t& endE1, size_t& startE2, size_t& endE2)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();

        hoNDArray<typename realType<T>::Type> mag(data.get_dimensions()), magSum, magSum2, magSumRO, magSumE1, magSumE2;
        Gadgetron::abs(data, mag);

        if ( NDim > 5 )
        {
            std::vector<size_t> dim;

            size_t ii;
            for ( ii=0; ii<NDim-5; ii++ )
            {
                mag.get_dimensions(dim);

                // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(mag, magSum));
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(mag, magSum, mag.get_number_of_dimensions() - 1));

                std::vector<size_t> dimSum(dim.size() - 1);
                memcpy(&dimSum[0], &dim[0], sizeof(size_t)*dimSum.size());
                magSum.reshape(dimSum);

                mag = magSum;
            }
        }

        size_t RO = mag.get_size(0);
        size_t E1 = mag.get_size(1);
        size_t E2 = mag.get_size(4);

        startRO = RO-1;
        endRO = 0;

        startE1 = E1-1;
        endE1 = 0;

        startE2 = E2-1;
        endE2 = 0;

        size_t ro, e1, e2;

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, magSum2, 4));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSum, 3));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum, magSum2, 2));

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSumE1, 1));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSumRO, 0));

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, magSum2, 3));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSum, 2));

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum, magSum2, 1));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(magSum2, magSumE2, 0));

        for ( ro=0; ro<RO; ro++ )
        {
            if ( magSumE1(ro) > 0 )
            {
                if ( ro < startRO ) startRO = ro;
                if ( ro > endRO ) endRO = ro;
            }
        }

        for ( e1=0; e1<E1; e1++ )
        {
            if ( magSumRO(e1) > 0 )
            {
                if ( e1 < startE1 ) startE1 = e1;
                if ( e1 > endE1 ) endE1 = e1;
            }
        }

        for ( e2=0; e2<E2; e2++ )
        {
            if ( magSumE2(e2) > 0 )
            {
                if ( e2 < startE2 ) startE2 = e2;
                if ( e2 > endE2 ) endE2 = e2;
            }
        }

        if ( startRO > endRO )
        {
            startRO = 0;
            endRO = RO-1;
        }

        if ( startE1 > endE1 )
        {
            startE1 = 0;
            endE1 = E1-1;
        }

        if ( startE2 > endE2 )
        {
            startE2 = 0;
            endE2 = E2-1;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledRegion3D(...) ... ");
        return false;
    }

    return true;
}

// ------------------------------------------------------------------------
// coil sensitivity
// ------------------------------------------------------------------------
template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
averageKSpace4D(const hoNDArray<T>& data, hoNDArray<T>& ave)
{
    try
    {
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, ave, 3));
        Gadgetron::scal( (typename realType<T>::Type)(1.0/data.get_size(3)), ave);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::averageKSpace4D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
averageKSpace4D(const hoNDArray<T>& data, hoNDArray<T>& ave, std::vector<size_t>& sampledTimes)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        if ( NDim < 4 )
        {
            ave = data;
            GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1(data, sampledTimes));
            return true;
        }

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t CHA = data.get_size(2);
        size_t N = data.get_size(3);

        hoNDArray<T> data4D(RO, E1, CHA, N, const_cast<T*>(data.begin()));
        GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1(data4D, sampledTimes));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, ave, 3));

        boost::shared_ptr< std::vector<size_t> > dim = ave.get_dimensions();

        if ( dim->size() != NDim )
        {
            (*dim).insert((*dim).begin()+3, 1);
            ave.reshape(dim.get());
        }

        hoNDArray<T> sampledTimes2D(RO, E1);
        T* pTimes = sampledTimes2D.begin();
        size_t ro, e1;
        for ( e1=0; e1<E1; e1++ )
        {
            double t = (double)sampledTimes[e1];
            if ( t == 0 ) t = 1;

            for ( ro=0; ro<RO; ro++ )
            {
                pTimes[e1*RO+ro] = (value_type)(1.0/t);
            }
        }

        // GADGET_CHECK_RETURN_FALSE(multipleMultiply(sampledTimes2D, ave, ave));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(multiply(ave, sampledTimes2D, ave));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::averageKSpace4D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
averageKSpace5D(const hoNDArray<T>& data, hoNDArray<T>& ave)
{
    try
    {
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, ave, 4));
        Gadgetron::scal( (typename realType<T>::Type)(1.0/data.get_size(4)), ave);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::averageKSpace5D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
averageKSpace5D(const hoNDArray<T>& data, hoNDArray<T>& ave, hoNDArray<size_t>& sampledTimes)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        if ( NDim < 5 )
        {
            ave = data;
            GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1E2(data, sampledTimes));
            return true;
        }

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t CHA = data.get_size(3);
        size_t N = data.get_size(4);

        hoNDArray<T> data5D(RO, E1, E2, CHA, N, const_cast<T*>(data.begin()));
        GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1E2(data5D, sampledTimes));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, ave, 4));

        hoNDArray<T> sampledTimes3D(RO, E1, E2);
        T* pTimes = sampledTimes3D.begin();
        size_t ro, e1, e2;
        for ( e2=0; e2<E2; e2++ )
        {
            for ( e1=0; e1<E1; e1++ )
            {
                double t = (double)sampledTimes(e1+e2*E1);
                if ( t == 0 ) t = 1;

                for ( ro=0; ro<RO; ro++ )
                {
                    pTimes[e2*RO*E1+e1*RO+ro] = (value_type)(1.0/t);
                }
            }
        }

        // GADGET_CHECK_RETURN_FALSE(multipleMultiply(sampledTimes3D, ave, ave));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(multiply(ave, sampledTimes3D, ave));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::averageKSpace5D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledTimesE1(const hoNDArray<T>& data4D, std::vector<size_t>& sampledTimes)
{
    try
    {
        size_t NDim = data4D.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = data4D.get_size(0);
        size_t E1 = data4D.get_size(1);
        size_t CHA = data4D.get_size(2);
        size_t N = data4D.get_size(3);

        hoNDArray<typename realType<T>::Type> mag(data4D.get_dimensions());
        Gadgetron::abs(data4D, mag);

        hoNDArray<typename realType<T>::Type> mag3D(RO, E1, 1, N);
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, mag3D, 2));

        hoNDArray<typename realType<T>::Type> mag2D(1, E1, 1, N);
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag3D, mag2D, 0));
        typename realType<T>::Type* pMag2D = mag2D.begin();

        sampledTimes.resize(E1, 0);

        size_t e1, n;
        for ( e1=0; e1<E1; e1++ )
        {
            for ( n=0; n<N; n++ )
            {
                if ( pMag2D[e1+n*E1] > 0 )
                {
                    sampledTimes[e1]++;
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledTimesE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledRegionE1(const hoNDArray<T>& data, size_t& startE1, size_t& endE1)
{
    try
    {
        std::vector<size_t> sampledTimes;
        GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1(data, sampledTimes));

        size_t E1 = sampledTimes.size();

        startE1 = E1-1;
        endE1 = 0;

        for ( size_t e1=0; e1<E1; e1++ )
        {
            if ( sampledTimes[e1] > 0 )
            {
                if ( e1 > endE1 ) endE1 = e1;
                if ( e1 < startE1 ) startE1 = e1;
            }
        }

        if ( endE1 < startE1 )
        {
            startE1 = 0;
            endE1 = E1-1;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledRegionE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledTimesE1E2(const hoNDArray<T>& data5D, hoNDArray<size_t>& sampledTimes)
{
    try
    {
        size_t NDim = data5D.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        size_t RO = data5D.get_size(0);
        size_t E1 = data5D.get_size(1);
        size_t E2 = data5D.get_size(2);
        size_t CHA = data5D.get_size(3);
        size_t N = data5D.get_size(4);

        hoNDArray<typename realType<T>::Type> mag(RO, E1, E2);

        hoNDArray<T> dataFirstChannel(RO, E1, E2, const_cast<T*>(data5D.begin()));
        Gadgetron::abs(dataFirstChannel, mag);

        hoNDArray<typename realType<T>::Type> mag3D(1, E1, E2);
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(mag, mag3D, 0));

        typename realType<T>::Type* pMag3D = mag3D.begin();

        sampledTimes.create(E1, E2);
        Gadgetron::clear(sampledTimes);
        size_t* pTimes = sampledTimes.get_data_ptr();

        size_t e1, e2, n;
        for ( e2=0; e2<E2; e2++ )
        {
            for ( e1=0; e1<E1; e1++ )
            {
                for ( n=0; n<N; n++ )
                {
                    if ( pMag3D[e1+e2*E1+n*E1*E2] > 0 )
                    {
                        pTimes[e1+e2*E1]++;
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledTimesE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
detectSampledRegionE1E2(const hoNDArray<T>& data, size_t& startE1, size_t& endE1, size_t& startE2, size_t& endE2)
{
    try
    {
        hoNDArray<size_t> sampledTimes;
        GADGET_CHECK_RETURN_FALSE(detectSampledTimesE1E2(data, sampledTimes));

        size_t E1 = sampledTimes.get_size(0);
        size_t E2 = sampledTimes.get_size(1);

        startE1 = E1-1;
        endE1 = 0;

        startE2 = E2-1;
        endE2 = 0;

        size_t e1, e2;
        for ( e2=0; e2<E2; e2++ )
        {
            for ( e1=0; e1<E1; e1++ )
            {
                if ( sampledTimes(e1+e2*E1) > 0 )
                {
                    if ( e1 > endE1 ) endE1 = e1;
                    if ( e1 < startE1 ) startE1 = e1;

                    if ( e2 > endE2 ) endE2 = e2;
                    if ( e2 < startE2 ) startE2 = e2;
                }
            }
        }

        if ( endE1 < startE1 )
        {
            startE1 = 0;
            endE1 = E1-1;
        }

        if ( endE2 < startE2 )
        {
            startE2 = 0;
            endE2 = E2-1;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::detectSampledRegionE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongE1(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startE1, size_t endE1)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongE1(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        size_t N = dst.get_number_of_elements()/(RO*E1);

        size_t n, e1;
        for ( n=0; n<N; n++ )
        {
            for ( e1=startE1; e1<=endE1; e1++ )
            {
                memcpy(dst.begin()+n*RO*E1+e1*RO, src.begin()+n*RO*E1+e1*RO, sizeof(T)*RO);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongROE1(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startRO>=RO) || (endRO>=RO) || (startRO>endRO) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
            return true;
        }

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        size_t N = dst.get_number_of_elements()/(RO*E1);
        const T* pSrc = src.begin();
        T* pDst = dst.begin();

        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, pSrc, pDst, RO, E1, startRO, endRO, startE1, endE1)
        for ( n=0; n<(long long)N; n++ )
        {
            for ( size_t e1=startE1; e1<=endE1; e1++ )
            {
                size_t offset = n*RO*E1+e1*RO+startRO;
                memcpy(pDst+offset, pSrc+offset, sizeof(T)*(endRO-startRO+1));
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongROE1(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongROE1E2(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);
        size_t E2 = dst.get_size(2);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);
        size_t E2_src = src.get_size(2);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(E2==E2_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startRO>=RO) || (endRO>=RO) || (startRO>endRO) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1E2(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
            return true;
        }

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1E2(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        if ( (startE2>=E2) || (endE2>=E2) || (startE2>endE2) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1E2(...) : (startE2>=E2) || (endE2>=E2) || (startE2>endE2) ... ");
            return true;
        }

        size_t N = dst.get_number_of_elements()/(RO*E1*E2);
        const T* pSrc = src.begin();
        T* pDst = dst.begin();

        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, pSrc, pDst, RO, E1, E2, startRO, endRO, startE1, endE1, startE2, endE2)
        for ( n=0; n<(long long)N; n++ )
        {
            for ( size_t e2=startE2; e2<=endE2; e2++ )
            {
                for ( size_t e1=startE1; e1<=endE1; e1++ )
                {
                    size_t offset = n*RO*E1*E2+e2*E1*RO+e1*RO+startRO;
                    memcpy(pDst+offset, pSrc+offset, sizeof(T)*(endRO-startRO+1));
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongROE1E2(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongROE1TransitionBand(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, 
        size_t startE1, size_t endE1, size_t transBandRO, size_t transBandE1)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=2);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startRO>=RO) || (endRO>=RO) || (startRO>endRO) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1TransitionBand(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
            return true;
        }

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1TransitionBand(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        while ( transBandRO>1 && startRO+transBandRO > RO/2 )
        {
             transBandRO--;
        }

        while ( transBandRO>1 && endRO-transBandRO < RO/2 )
        {
             transBandRO--;
        }

        while ( transBandE1>1 && startE1+transBandE1 > E1/2 )
        {
             transBandE1--;
        }

        while ( transBandE1>1 && endE1-transBandE1 < E1/2 )
        {
             transBandE1--;
        }

        ISMRMRDKSPACEFILTER filterType = ISMRMRD_FILTER_TAPERED_HANNING;
        bool densityComp = false;

        hoNDArray<T> filter_src_RO, filter_src_E1;

        if ( startRO==0 && endRO==RO-1 )
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_NONE, transBandRO, densityComp));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(RO, startRO, endRO, filter_src_RO, filterType, transBandRO, densityComp));
        }

        if ( startE1==0 && endE1==E1-1 )
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_NONE, transBandE1, densityComp));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(E1, startE1, endE1, filter_src_E1, filterType, transBandE1, densityComp));
        }

        // in this way, the SNR unit scale property is perserved
        T midValue = filter_src_RO(RO/2);
        T scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_RO);

        midValue = filter_src_E1(E1/2);
        scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_E1);

        hoNDArray<T> filter_dst_RO(RO), filter_dst_E1(E1);

        size_t ii;
        for ( ii=0; ii<RO; ii++ )
        {
            filter_dst_RO(ii) = T(1.0) - filter_src_RO(ii);
        }

        for ( ii=0; ii<E1; ii++ )
        {
            filter_dst_E1(ii) = T(1.0) - filter_src_E1(ii);
        }

        hoNDArray<T> srcFiltered(src), dstFiltered(dst);
        if ( startRO==0 && endRO==RO-1 )
        {
            GADGET_CHECK_RETURN_FALSE(kspacefilterE1(src, filter_src_E1, srcFiltered));
            GADGET_CHECK_RETURN_FALSE(kspacefilterE1(dst, filter_dst_E1, dstFiltered));
        }
        else if ( startE1==0 && endE1==E1-1 )
        {
            GADGET_CHECK_RETURN_FALSE(kspacefilterRO(src, filter_src_RO, srcFiltered));
            GADGET_CHECK_RETURN_FALSE(kspacefilterRO(dst, filter_dst_RO, dstFiltered));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(kspacefilterROE1(src, filter_src_RO, filter_src_E1, srcFiltered));

            hoNDArray<T> fxy;
            GADGET_CHECK_RETURN_FALSE(compute2DFilterFromTwo1D(filter_src_RO, filter_src_E1, fxy));

            size_t Nxy = RO*E1;
            for ( ii=0; ii<Nxy; ii++ )
            {
                fxy(ii) = T(1.0) - fxy(ii);
            }

            GADGET_CHECK_RETURN_FALSE(kspacefilterROE1(dst, fxy, dstFiltered));
        }

        Gadgetron::add(srcFiltered, dstFiltered, dst);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongROE1TransitionBand(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
copyAlongROE1E2TransitionBand(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, 
                        size_t startE1, size_t endE1, size_t startE2, size_t endE2, 
                        size_t transBandRO, size_t transBandE1, size_t transBandE2)
{
    try
    {
        size_t NDim = src.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        size_t RO = dst.get_size(0);
        size_t E1 = dst.get_size(1);
        size_t E2 = dst.get_size(2);

        size_t RO_src = src.get_size(0);
        size_t E1_src = src.get_size(1);
        size_t E2_src = src.get_size(2);

        GADGET_CHECK_RETURN_FALSE(RO==RO_src);
        GADGET_CHECK_RETURN_FALSE(E1==E1_src);
        GADGET_CHECK_RETURN_FALSE(E2==E2_src);
        GADGET_CHECK_RETURN_FALSE(src.get_number_of_elements()==dst.get_number_of_elements());

        if ( (startRO>=RO) || (endRO>=RO) || (startRO>endRO) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1TransitionBand(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
            return true;
        }

        if ( (startE1>=E1) || (endE1>=E1) || (startE1>endE1) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1TransitionBand(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
            return true;
        }

        if ( (startE2>=E2) || (endE2>=E2) || (startE2>endE2) )
        {
            dst = src;
            GWARN_STREAM("copyAlongROE1E2TransitionBand(...) : (startE2>=E2) || (endE2>=E2) || (startE2>endE2) ... ");
            return true;
        }

        while ( transBandRO>1 && startRO+transBandRO > RO/2 )
        {
             transBandRO--;
        }

        while ( transBandRO>1 && endRO-transBandRO < RO/2 )
        {
             transBandRO--;
        }

        while ( transBandE1>1 && startE1+transBandE1 > E1/2 )
        {
             transBandE1--;
        }

        while ( transBandE1>1 && endE1-transBandE1 < E1/2 )
        {
             transBandE1--;
        }

        while ( transBandE2>1 && startE2+transBandE2 > E2/2 )
        {
             transBandE2--;
        }

        while ( transBandE2>1 && endE2-transBandE2 < E2/2 )
        {
             transBandE2--;
        }

        ISMRMRDKSPACEFILTER filterType = ISMRMRD_FILTER_TAPERED_HANNING;
        bool densityComp = false;

        hoNDArray<T> filter_src_RO, filter_src_E1, filter_src_E2;

        if ( startRO==0 && endRO==RO-1 )
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_NONE, transBandRO, densityComp));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(RO, startRO, endRO, filter_src_RO, filterType, transBandRO, densityComp));
        }

        if ( startE1==0 && endE1==E1-1 )
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_NONE, transBandE1, densityComp));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(E1, startE1, endE1, filter_src_E1, filterType, transBandE1, densityComp));
        }

        if ( startE2==0 && endE2==E2-1 )
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(E2, startE2, endE2, filter_src_E2, ISMRMRD_FILTER_NONE, transBandE2, densityComp));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(generateAsymmetricFilter(E2, startE2, endE2, filter_src_E2, filterType, transBandE2, densityComp));
        }

        // in this way, the SNR unit scale property is perserved
        T midValue = filter_src_RO(RO/2);
        T scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_RO);

        midValue = filter_src_E1(E1/2);
        scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_E1);

        midValue = filter_src_E2(E2/2);
        scalFactor = T(1.0)/midValue;
        Gadgetron::scal(scalFactor, filter_src_E2);

        hoNDArray<T> filter_dst_RO(RO), filter_dst_E1(E1), filter_dst_E2(E2);

        size_t ii;
        for ( ii=0; ii<RO; ii++ )
        {
            filter_dst_RO(ii) = T(1.0) - filter_src_RO(ii);
        }

        for ( ii=0; ii<E1; ii++ )
        {
            filter_dst_E1(ii) = T(1.0) - filter_src_E1(ii);
        }

        for ( ii=0; ii<E2; ii++ )
        {
            filter_dst_E2(ii) = T(1.0) - filter_src_E2(ii);
        }

        hoNDArray<T> srcFiltered(src), dstFiltered(dst);
        if ( startRO>=0 && endRO<=RO-1 && startE1==0 && endE1==E1-1 && startE2==0 && endE1==E2-1 )
        {
            GADGET_CHECK_RETURN_FALSE(kspacefilterRO(src, filter_src_E1, srcFiltered));
            GADGET_CHECK_RETURN_FALSE(kspacefilterRO(dst, filter_dst_E1, dstFiltered));
        }
        else if ( startRO==0 && endRO==RO-1 && startE1>=0 && endE1<=E1-1 && startE2==0 && endE1==E2-1 )
        {
            GADGET_CHECK_RETURN_FALSE(kspacefilterE1(src, filter_src_RO, srcFiltered));
            GADGET_CHECK_RETURN_FALSE(kspacefilterE1(dst, filter_dst_RO, dstFiltered));
        }
        else if ( startRO==0 && endRO==RO-1 && startE1==0 && endE1==E1-1 && startE2>=0 && endE1<=E2-1 )
        {
            GADGET_CHECK_RETURN_FALSE(kspace3DfilterE2(src, filter_src_RO, srcFiltered));
            GADGET_CHECK_RETURN_FALSE(kspace3DfilterE2(dst, filter_dst_RO, dstFiltered));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(kspace3DfilterROE1E2(src, filter_src_RO, filter_src_E1, filter_src_E2, srcFiltered));

            hoNDArray<T> fxyz;
            GADGET_CHECK_RETURN_FALSE(compute3DFilterFromThree1D(filter_src_RO, filter_src_E1, filter_src_E2, fxyz));

            size_t Nxyz = RO*E1*E2;
            for ( ii=0; ii<Nxyz; ii++ )
            {
                fxyz(ii) = T(1.0) - fxyz(ii);
            }

            GADGET_CHECK_RETURN_FALSE(kspace3DfilterROE1E2(dst, fxyz, dstFiltered));
        }

        Gadgetron::add(srcFiltered, dstFiltered, dst);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::copyAlongROE1E2TransitionBand(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
std::string gtPlusISMRMRDReconUtil<T>::getISMRMRDDimName(const ISMRMRDDIM& dim)
{
    std::ostringstream os;
    switch (dim)
    {
        case DIM_ReadOut:
            os << "DIM_ReadOut";
        break;

        case DIM_Encoding1:
            os << "DIM_Encoding1";
        break;

        case DIM_Channel:
            os << "DIM_Channel";
        break;

        case DIM_Slice:
            os << "DIM_Slice";
        break;

        case DIM_Encoding2:
            os << "DIM_Encoding2";
        break;

        case DIM_Contrast:
            os << "DIM_Contrast";
        break;

        case DIM_Phase:
            os << "DIM_Phase";
        break;

        case DIM_Repetition:
            os << "DIM_Repetition";
        break;

        case DIM_Set:
            os << "DIM_Set";
        break;

        case DIM_Segment:
            os << "DIM_Segment";
        break;

        case DIM_Average:
            os << "DIM_Average";
        break;

        case DIM_other1:
            os << "DIM_other1";
        break;

        case DIM_other2:
            os << "DIM_other2";
        break;

        case DIM_other3:
            os << "DIM_other3";
        break;

        default:
            os << "DIM_NONE";
    }

    std::string dimStr(os.str());
    return dimStr;
}

template <typename T> 
ISMRMRDDIM gtPlusISMRMRDReconUtil<T>::getISMRMRDDimFromName(const std::string& name)
{
    if ( name == "DIM_ReadOut" ) return DIM_ReadOut;
    if ( name == "DIM_Encoding1" ) return DIM_Encoding1;
    if ( name == "DIM_Channel" ) return DIM_Channel;
    if ( name == "DIM_Slice" ) return DIM_Slice;
    if ( name == "DIM_Encoding2" ) return DIM_Encoding2;
    if ( name == "DIM_Contrast" ) return DIM_Contrast;
    if ( name == "DIM_Phase" ) return DIM_Phase;
    if ( name == "DIM_Repetition" ) return DIM_Repetition;
    if ( name == "DIM_Set" ) return DIM_Set;
    if ( name == "DIM_Segment" ) return DIM_Segment;
    if ( name == "DIM_Average" ) return DIM_Average;
    if ( name == "DIM_other1" ) return DIM_other1;
    if ( name == "DIM_other2" ) return DIM_other2;
    if ( name == "DIM_other3" ) return DIM_other3;

    return DIM_NONE;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::getISMRMRDDimIndex(const ISMRMRDDIM& dim, long long& ind)
{
    switch (dim)
    {
        case Gadgetron::DIM_ReadOut:
            ind = 0;
        break;

        case Gadgetron::DIM_Encoding1:
            ind = 1;
        break;

        case Gadgetron::DIM_Channel:
            ind = 2;
        break;

        case Gadgetron::DIM_Slice:
            ind = 3;
        break;

        case Gadgetron::DIM_Encoding2:
            ind = 4;
        break;

        case Gadgetron::DIM_Contrast:
            ind = 5;
        break;

        case Gadgetron::DIM_Phase:
            ind = 6;
        break;

        case Gadgetron::DIM_Repetition:
            ind = 7;
        break;

        case Gadgetron::DIM_Set:
            ind = 8;
        break;

        case Gadgetron::DIM_Segment:
            ind = 9;
        break;

        case Gadgetron::DIM_Average:
            ind = 10;
        break;

        default:
            ind = -1;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::findDimIndex(const std::vector<DimensionRecordType>& dimStartingIndexes, ISMRMRDDIM dim, size_t ind)
{
    size_t N = dimStartingIndexes.size();

    size_t n;
    for ( n=0; n<N; n++ )
    {
        if ( dimStartingIndexes[n].first == dim )
        {
            ind = dimStartingIndexes[n].second;
            return true;
        }
    }

    return false;
}

template <typename T> 
ISMRMRDALGO gtPlusISMRMRDReconUtil<T>::getISMRMRDReconAlgoFromName(const std::string& name)
{
    if ( name == "ISMRMRD_GRAPPA" ) return ISMRMRD_GRAPPA;
    if ( name == "ISMRMRD_SENSE" ) return ISMRMRD_SENSE;
    if ( name == "ISMRMRD_SPIRIT" ) return ISMRMRD_SPIRIT;
    if ( name == "ISMRMRD_L1SPIRIT" ) return ISMRMRD_L1SPIRIT;
    if ( name == "ISMRMRD_SOFTSENSE" ) return ISMRMRD_SOFTSENSE;
    if ( name == "ISMRMRD_L1SOFTSENSE" ) return ISMRMRD_L1SOFTSENSE;
    if ( name == "ISMRMRD_2DTBINNING" ) return ISMRMRD_2DTBINNING;
    if ( name == "ISMRMRD_2DTBINNING_FLOW" ) return ISMRMRD_2DTBINNING_FLOW;
    if ( name == "ISMRMRD_L1SPIRIT_SLEP" ) return ISMRMRD_L1SPIRIT_SLEP;
    if ( name == "ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP" ) return ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP;

    return ISMRMRD_NONE;
}

template <typename T> 
ISMRMRDCOILMAPALGO gtPlusISMRMRDReconUtil<T>::getISMRMRDCoilMapAlgoFromName(const std::string& name)
{
    if ( name == "ISMRMRD_SOUHEIL" ) return ISMRMRD_SOUHEIL;
    if ( name == "ISMRMRD_SOUHEIL_ITER" ) return ISMRMRD_SOUHEIL_ITER;

    return ISMRMRD_SOUHEIL;
}

template <typename T> 
ISMRMRDPFALGO gtPlusISMRMRDReconUtil<T>::getISMRMRDPartialFourierReconAlgoFromName(const std::string& name)
{
    if ( name == "ISMRMRD_PF_HOMODYNE" ) return ISMRMRD_PF_HOMODYNE;
    if ( name == "ISMRMRD_PF_FENGHUANG" ) return ISMRMRD_PF_FENGHUANG;
    if ( name == "ISMRMRD_PF_POCS" ) return ISMRMRD_PF_POCS;
    if ( name == "ISMRMRD_PF_ZEROFILLING_FILTER" ) return ISMRMRD_PF_ZEROFILLING_FILTER;
    if ( name == "ISMRMRD_PF_ZEROFILLING" ) return ISMRMRD_PF_ZEROFILLING;

    return ISMRMRD_PF_NONE;
}

template <typename T> 
std::string gtPlusISMRMRDReconUtil<T>::getNameFromISMRMRDPartialFourierReconAlgo(ISMRMRDPFALGO algo)
{
    if ( algo == ISMRMRD_PF_HOMODYNE ) return std::string("ISMRMRD_PF_HOMODYNE");
    if ( algo == ISMRMRD_PF_FENGHUANG ) return std::string("ISMRMRD_PF_FENGHUANG");
    if ( algo == ISMRMRD_PF_ZEROFILLING_FILTER ) return std::string("ISMRMRD_PF_ZEROFILLING_FILTER");
    if ( algo == ISMRMRD_PF_POCS ) return std::string("ISMRMRD_PF_POCS");
    if ( algo == ISMRMRD_PF_ZEROFILLING ) return std::string("ISMRMRD_PF_ZEROFILLING");

    return std::string("ISMRMRD_PF_NONE");
}

template <typename T> 
ISMRMRDKSPACEFILTER gtPlusISMRMRDReconUtil<T>::
getISMRMRDKSpaceFilterFromName(const std::string& name)
{
    if ( name == "ISMRMRD_FILTER_GAUSSIAN" ) return ISMRMRD_FILTER_GAUSSIAN;
    if ( name == "ISMRMRD_FILTER_HANNING" ) return ISMRMRD_FILTER_HANNING;
    if ( name == "ISMRMRD_FILTER_TUKEY" ) return ISMRMRD_FILTER_TUKEY;
    if ( name == "ISMRMRD_FILTER_TAPERED_HANNING" ) return ISMRMRD_FILTER_TAPERED_HANNING;
    if ( name == "ISMRMRD_FILTER_NONE" ) return ISMRMRD_FILTER_NONE;

    return ISMRMRD_FILTER_NONE;
}

template <typename T> 
ISMRMRDINTERPRETROGATING gtPlusISMRMRDReconUtil<T>::getISMRMRDRetroGatingInterpFromName(const std::string& name)
{
    if ( name == "ISMRMRD_INTERP_RETRO_GATING_LINEAR" ) return ISMRMRD_INTERP_RETRO_GATING_LINEAR;
    if ( name == "ISMRMRD_INTERP_RETRO_GATING_CUBIC" ) return ISMRMRD_INTERP_RETRO_GATING_CUBIC;
    if ( name == "ISMRMRD_INTERP_RETRO_GATING_BSPLINE" ) return ISMRMRD_INTERP_RETRO_GATING_BSPLINE;

    return ISMRMRD_INTERP_RETRO_GATING_LINEAR;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
extractSubArrayForDim(const hoNDArray<T>& x, hoNDArray<T>& r, ISMRMRDDIM& dim, size_t value, bool lessEqual)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        long long dimInd;
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim, dimInd));

        GADGET_CHECK_RETURN_FALSE(value<(*dimX)[dimInd]);

        std::vector<size_t> crop_offset(11, 0);
        crop_offset[0] = 0;
        crop_offset[1] = 0;
        crop_offset[2] = 0;
        crop_offset[3] = 0;
        crop_offset[4] = 0;
        crop_offset[5] = 0;
        crop_offset[6] = 0;
        crop_offset[7] = 0;
        crop_offset[8] = 0;
        crop_offset[9] = 0;
        crop_offset[10] = 0;

        std::vector<size_t> crop_size(11, 0);
        crop_size[0] = (*dimX)[0];
        crop_size[1] = (*dimX)[1];
        crop_size[2] = (*dimX)[2];
        crop_size[3] = (*dimX)[3];
        crop_size[4] = (*dimX)[4];
        crop_size[5] = (*dimX)[5];
        crop_size[6] = (*dimX)[6];
        crop_size[7] = (*dimX)[7];
        crop_size[8] = (*dimX)[8];
        crop_size[9] = (*dimX)[9];
        crop_size[10] = (*dimX)[10];

        if ( lessEqual )
        {
            crop_size[dimInd] = value+1;
        }
        else
        {
            crop_offset[dimInd] = value;
            crop_size[dimInd] = 1;
        }

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(x, r, crop_offset, crop_size));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::extractSubArrayForDim(dim, value) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
extractSubArrayForDim(const hoNDArray<T>& x, hoNDArray<T>& r, ISMRMRDDIM& dim1, size_t value1, ISMRMRDDIM& dim2, size_t value2, bool lessEqual)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        long long dimInd1, dimInd2;
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim1, dimInd1));
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim2, dimInd2));

        GADGET_CHECK_RETURN_FALSE(value1<(*dimX)[dimInd1]);
        GADGET_CHECK_RETURN_FALSE(value2<(*dimX)[dimInd2]);

        std::vector<size_t> crop_offset(11, 0);
        crop_offset[0] = 0;
        crop_offset[1] = 0;
        crop_offset[2] = 0;
        crop_offset[3] = 0;
        crop_offset[4] = 0;
        crop_offset[5] = 0;
        crop_offset[6] = 0;
        crop_offset[7] = 0;
        crop_offset[8] = 0;
        crop_offset[9] = 0;
        crop_offset[10] = 0;

        std::vector<size_t> crop_size(11, 0);
        crop_size[0] = (*dimX)[0];
        crop_size[1] = (*dimX)[1];
        crop_size[2] = (*dimX)[2];
        crop_size[3] = (*dimX)[3];
        crop_size[4] = (*dimX)[4];
        crop_size[5] = (*dimX)[5];
        crop_size[6] = (*dimX)[6];
        crop_size[7] = (*dimX)[7];
        crop_size[8] = (*dimX)[8];
        crop_size[9] = (*dimX)[9];
        crop_size[10] = (*dimX)[10];

        if ( lessEqual )
        {
            crop_size[dimInd1] = value1+1;
            crop_size[dimInd2] = value2+1;
        }
        else
        {
            crop_offset[dimInd1] = value1;
            crop_size[dimInd1] = 1;

            crop_offset[dimInd2] = value2;
            crop_size[dimInd2] = 1;
        }

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(x, r, crop_offset, crop_size));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::extractSubArrayForDim(dim1, value1, dim2, value2) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
extractSubArrayForDim1LessEqualDim2Equal(const hoNDArray<T>& x, hoNDArray<T>& r, ISMRMRDDIM& dim1, size_t value1, ISMRMRDDIM& dim2, size_t value2)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        long long dimInd1, dimInd2;
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim1, dimInd1));
        GADGET_CHECK_RETURN_FALSE(getISMRMRDDimIndex(dim2, dimInd2));

        GADGET_CHECK_RETURN_FALSE(value1<(*dimX)[dimInd1]);
        GADGET_CHECK_RETURN_FALSE(value2<(*dimX)[dimInd2]);

        std::vector<size_t> crop_offset(11, 0);
        crop_offset[0] = 0;
        crop_offset[1] = 0;
        crop_offset[2] = 0;
        crop_offset[3] = 0;
        crop_offset[4] = 0;
        crop_offset[5] = 0;
        crop_offset[6] = 0;
        crop_offset[7] = 0;
        crop_offset[8] = 0;
        crop_offset[9] = 0;
        crop_offset[10] = 0;

        std::vector<size_t> crop_size(11, 0);
        crop_size[0] = (*dimX)[0];
        crop_size[1] = (*dimX)[1];
        crop_size[2] = (*dimX)[2];
        crop_size[3] = (*dimX)[3];
        crop_size[4] = (*dimX)[4];
        crop_size[5] = (*dimX)[5];
        crop_size[6] = (*dimX)[6];
        crop_size[7] = (*dimX)[7];
        crop_size[8] = (*dimX)[8];
        crop_size[9] = (*dimX)[9];
        crop_size[10] = (*dimX)[9];

        crop_size[dimInd1] = value1+1;

        crop_offset[dimInd2] = value2;
        crop_size[dimInd2] = 1;

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(x, r, crop_offset, crop_size));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::extractSubArrayForDim1LessEqualDim2Equal(dim1, value1, dim2, value2) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::
extractSubArrayForMaxEncodingCounters(const hoNDArray<T>& x, hoNDArray<T>& r, const ISMRMRD::EncodingCounters& maxIdx)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimX = x.get_dimensions();

        std::vector<size_t> crop_offset(11, 0);
        crop_offset[0] = 0;
        crop_offset[1] = 0;
        crop_offset[2] = 0;
        crop_offset[3] = 0;
        crop_offset[4] = 0;
        crop_offset[5] = 0;
        crop_offset[6] = 0;
        crop_offset[7] = 0;
        crop_offset[8] = 0;
        crop_offset[9] = 0;
        crop_offset[10] = 0;

        // [RO E1 Cha Slice E2 Contrast Phase Rep Set Seg Ave]
        std::vector<size_t> crop_size(11, 0);
        crop_size[0] = (*dimX)[0];
        crop_size[1] = (*dimX)[1]; if ( maxIdx.kspace_encode_step_1 < crop_size[1]-1 ) crop_size[1] = maxIdx.kspace_encode_step_1+1;
        crop_size[2] = (*dimX)[2]; 
        crop_size[3] = (*dimX)[3]; if ( maxIdx.slice                < crop_size[3]-1 ) crop_size[3] = maxIdx.slice+1;
        crop_size[4] = (*dimX)[4]; if ( maxIdx.kspace_encode_step_2 < crop_size[4]-1 ) crop_size[4] = maxIdx.kspace_encode_step_2+1;
        crop_size[5] = (*dimX)[5]; if ( maxIdx.contrast             < crop_size[5]-1 ) crop_size[5] = maxIdx.contrast+1;
        crop_size[6] = (*dimX)[6]; if ( maxIdx.phase                < crop_size[6]-1 ) crop_size[6] = maxIdx.phase+1;
        crop_size[7] = (*dimX)[7]; if ( maxIdx.repetition           < crop_size[7]-1 ) crop_size[7] = maxIdx.repetition+1;
        crop_size[8] = (*dimX)[8]; if ( maxIdx.set                  < crop_size[8]-1 ) crop_size[8] = maxIdx.set+1;
        crop_size[9] = (*dimX)[9]; if ( maxIdx.segment              < crop_size[9]-1 ) crop_size[9] = maxIdx.segment+1;
        crop_size[10] = (*dimX)[10]; if ( maxIdx.average            < crop_size[10]-1 ) crop_size[10] = maxIdx.average+1;

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(x, r, crop_offset, crop_size));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::extractSubArrayForMaxEncodingCounters(const hoNDArray<T>& x, hoNDArray<T>& r, const ISMRMRD::EncodingCounters& maxIdx) ... ");
        return false;
    }

    return true;
}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::clearAcquisitionHeaderISMRMRD(ISMRMRD::AcquisitionHeader& acqHeader)
{
    memset(&acqHeader, 0, sizeof(ISMRMRD::AcquisitionHeader));
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::hasIdenticalGeometryISMRMRD(const ISMRMRD::AcquisitionHeader& acqHeader1, const ISMRMRD::AcquisitionHeader& acqHeader2)
{
    long long ii;

    for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
    {
        if ( std::abs(acqHeader1.position[ii]-acqHeader2.position[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
        if ( std::abs(acqHeader1.patient_table_position[ii]-acqHeader2.patient_table_position[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
    }

    for ( ii=0; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
    {
        if ( std::abs(acqHeader1.read_dir[ii]-acqHeader2.read_dir[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
        if ( std::abs(acqHeader1.phase_dir[ii]-acqHeader2.phase_dir[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
        if ( std::abs(acqHeader1.slice_dir[ii]-acqHeader2.slice_dir[ii]) > GT_IMAGING_GEOMETRY_DELTA ) return false;
    }

    return true;
}

template <typename T> 
long long gtPlusISMRMRDReconUtil<T>::addPrePostZeros(size_t centre_column, size_t samples)
{
    // 1 : pre zeros
    // 2 : post zeros
    // 0 : no zeros
    if ( 2*centre_column == samples )
    {
        return 0;
    }

    if ( 2*centre_column < samples )
    {
        return 1;
    }

    if ( 2*centre_column > samples )
    {
        return 2;
    }

    return 0;
}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::findStartEndRO(size_t centre_column, size_t samples, long long& startRO, long long& endRO)
{
    long long zerosFlag = addPrePostZeros(centre_column, samples);

    if ( zerosFlag == 0 )
    {
        startRO = 0;
        endRO = (long long)samples-1;
    }

    if ( zerosFlag == 1 )
    {
        endRO = (long long)2*(samples-centre_column)-1;
        startRO = (long long)endRO-samples+1;
    }

    if ( zerosFlag == 2 )
    {
        startRO = 0;
        endRO = (long long)samples-1;
    }

    return;
}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::findStartEndROAfterZeroFilling(size_t centre_column, size_t samples_zerofilled, int& startRO, int& endRO)
{
    size_t num = samples_zerofilled/2;

    if ( centre_column == num )
    {
        startRO = 0;
        endRO = (int)samples_zerofilled-1;
    }

    if ( centre_column+num < samples_zerofilled ) // pre zeros
    {
        endRO = (int)samples_zerofilled-1;
        startRO = endRO-(int)(centre_column+num)+1;
    }

    if ( centre_column+num > samples_zerofilled ) // post zeros
    {
        startRO = 0;
        endRO = (int)samples_zerofilled-1;
    }

    return;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::setMetaAttributesFromImageHeaderISMRMRD(const ISMRMRD::ImageHeader& imgHeader, ISMRMRD::MetaContainer& attrib)
{
    try
    {
        unsigned int ii;

        attrib.set(ISMRMRD_IMAGE_version,                 (long)imgHeader.version);
        attrib.set(ISMRMRD_IMAGE_flags,                   (long)imgHeader.flags);
        attrib.set(ISMRMRD_IMAGE_measurement_uid,         (long)imgHeader.measurement_uid);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_matrix_size, (long)imgHeader.matrix_size[0]);
        attrib.append(ISMRMRD_IMAGE_matrix_size, (long)imgHeader.matrix_size[1]);
        attrib.append(ISMRMRD_IMAGE_matrix_size, (long)imgHeader.matrix_size[2]);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_field_of_view, (double)imgHeader.field_of_view[0]);
        attrib.append(ISMRMRD_IMAGE_field_of_view, (double)imgHeader.field_of_view[1]);
        attrib.append(ISMRMRD_IMAGE_field_of_view, (double)imgHeader.field_of_view[2]);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_channels, (long)imgHeader.channels);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_position, (double)imgHeader.position[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_position, (double)imgHeader.position[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_read_dir, (double)imgHeader.read_dir[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_read_dir, (double)imgHeader.read_dir[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_phase_dir, (double)imgHeader.phase_dir[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_phase_dir, (double)imgHeader.phase_dir[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_slice_dir, (double)imgHeader.slice_dir[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_slice_dir, (double)imgHeader.slice_dir[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_patient_table_position, (double)imgHeader.patient_table_position[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_patient_table_position, (double)imgHeader.patient_table_position[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_average,       (long)imgHeader.average);
        attrib.set(ISMRMRD_IMAGE_slice,         (long)imgHeader.slice);
        attrib.set(ISMRMRD_IMAGE_contrast,      (long)imgHeader.contrast);
        attrib.set(ISMRMRD_IMAGE_phase,         (long)imgHeader.phase);
        attrib.set(ISMRMRD_IMAGE_repetition,    (long)imgHeader.repetition);
        attrib.set(ISMRMRD_IMAGE_set,           (long)imgHeader.set);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_acquisition_time_stamp, (long)imgHeader.acquisition_time_stamp);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_physiology_time_stamp, (long)imgHeader.physiology_time_stamp[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_PHYS_STAMPS; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_physiology_time_stamp, (long)imgHeader.physiology_time_stamp[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_image_data_type,       (long)imgHeader.data_type);
        attrib.set(ISMRMRD_IMAGE_image_type,            (long)imgHeader.image_type);
        attrib.set(ISMRMRD_IMAGE_image_series_index,    (long)imgHeader.image_series_index);

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_user_int, (long)imgHeader.user_int[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_USER_INTS; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_user_int, (long)imgHeader.user_int[ii]);
        }

        // ----------------------------------------

        attrib.set(ISMRMRD_IMAGE_user_float, (double)imgHeader.user_float[0]);
        for ( ii=1; ii<ISMRMRD::ISMRMRD_USER_FLOATS; ii++ )
        {
            attrib.append(ISMRMRD_IMAGE_user_float, (double)imgHeader.user_float[ii]);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::setMetaAttributesFromImageHeaderISMRMRD(const ISMRMRD::ImageHeader& imgHeader, ISMRMRD::MetaContainer& attrib) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtil<T>::setImageHeaderISMRMRDFromMetaAttributes(const ISMRMRD::MetaContainer& attrib, ISMRMRD::ImageHeader& imgHeader)
{
    try
    {
        unsigned int ii;

        imgHeader.version = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_version, 0);
        imgHeader.flags = (uint64_t)attrib.as_long(ISMRMRD_IMAGE_flags, 0);
        imgHeader.measurement_uid = (uint32_t)attrib.as_long(ISMRMRD_IMAGE_measurement_uid, 0);

        // ----------------------------------------

        imgHeader.matrix_size[0] = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_matrix_size, 0);
        imgHeader.matrix_size[1] = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_matrix_size, 1);
        imgHeader.matrix_size[2] = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_matrix_size, 2);

        // ----------------------------------------

        imgHeader.field_of_view[0] = (float)attrib.as_double(ISMRMRD_IMAGE_field_of_view, 0);
        imgHeader.field_of_view[1] = (float)attrib.as_double(ISMRMRD_IMAGE_field_of_view, 1);
        imgHeader.field_of_view[2] = (float)attrib.as_double(ISMRMRD_IMAGE_field_of_view, 2);

        // ----------------------------------------

        imgHeader.channels = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_channels, 0);;

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
        {
            imgHeader.position[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_position, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            imgHeader.read_dir[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_read_dir, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            imgHeader.phase_dir[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_phase_dir, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_DIRECTION_LENGTH; ii++ )
        {
            imgHeader.slice_dir[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_slice_dir, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_POSITION_LENGTH; ii++ )
        {
            imgHeader.patient_table_position[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_patient_table_position, ii);
        }

        // ----------------------------------------

        imgHeader.average = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_average, 0);
        imgHeader.slice = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_slice, 0);
        imgHeader.contrast = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_contrast, 0);
        imgHeader.phase = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_phase, 0);
        imgHeader.repetition = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_repetition, 0);
        imgHeader.set = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_set, 0);

        // ----------------------------------------

        imgHeader.acquisition_time_stamp = (uint32_t)attrib.as_long(ISMRMRD_IMAGE_acquisition_time_stamp, 0);

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_PHYS_STAMPS; ii++ )
        {
            imgHeader.physiology_time_stamp[ii] = (uint32_t)attrib.as_long(ISMRMRD_IMAGE_physiology_time_stamp, ii);
        }

        // ----------------------------------------

        imgHeader.data_type = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_image_data_type, 0);
        imgHeader.image_type = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_image_type, 0);
        imgHeader.image_series_index = (uint16_t)attrib.as_long(ISMRMRD_IMAGE_image_series_index, 0);

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_USER_INTS; ii++ )
        {
            imgHeader.user_int[ii] = (int32_t)attrib.as_long(ISMRMRD_IMAGE_user_int, ii);
        }

        // ----------------------------------------

        for ( ii=0; ii<ISMRMRD::ISMRMRD_USER_FLOATS; ii++ )
        {
            imgHeader.user_float[ii] = (float)attrib.as_double(ISMRMRD_IMAGE_user_float, ii);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::setImageHeaderISMRMRDFromMetaAttributes(const ISMRMRD::MetaContainer& attrib, ISMRMRD::ImageHeader& imgHeader) ... ");
        return false;
    }

    return true;
}

//#ifdef USE_CUDA
//
//template <typename T> 
//bool gtPlusISMRMRDReconUtil<T>::
//cudaJobSplitter(const std::vector<unsigned int>& jobIDs, size_t jobSize, size_t minimalMemoryForValidDevice, 
//                std::vector< std::pair<unsigned int, std::vector<std::vector<unsigned int> > > >& jobSchedule)
//{
//    try
//    {
//        unsigned int numOfJobs = jobIDs.size();
//        if ( numOfJobs == 0 )
//        {
//            GWARN_STREAM("numOfJobs == 0");
//            return true;
//        }
//
//        // find valid device
//        int numOfDevices(0);
//        GADGET_CHECK_RETURN_FALSE(cudaGetDeviceCount( &numOfDevices )==cudaSuccess);
//
//        if ( numOfDevices == 0 )
//        {
//            GWARN_STREAM("numOfDevices == 0");
//            return true;
//        }
//
//        std::vector<unsigned int> validDevices;
//        int d;
//        for ( d=0; d<numOfDevices; d++ )
//        {
//            size_t totalMem = cudaDeviceManager::Instance()->total_global_mem(d);
//            if ( totalMem >= minimalMemoryForValidDevice )
//            {
//                validDevices.push_back(d);
//            }
//        }
//
//        if ( validDevices.empty() )
//        {
//            GERROR_STREAM("No valid device can be found : " << minimalMemoryForValidDevice);
//            return false;
//        }
//
//        std::vector<unsigned int> maxJobN(validDevices.size());
//        for ( d=0; d<validDevices.size(); d++ )
//        {
//            size_t totalMem = cudaDeviceManager::Instance()->total_global_mem(validDevices[d]);
//            maxJobN[d] = totalMem/jobSize;
//        }
//
//        jobSchedule.clear();
//
//        size_t job = 0;
//        unsigned int validDevice = 0;
//        while ( job < numOfJobs )
//        {
//            size_t start = job;
//            size_t end = job + maxJobN[validDevice] - 1;
//
//            if ( end >= numOfJobs ) end = numOfJobs - 1;
//
//            unsigned int deviceID = validDevices[validDevice];
//
//            unsigned int loc;
//            for ( loc=0; loc<jobSchedule.size(); loc++ )
//            {
//                if ( jobSchedule[loc].first == deviceID ) break;
//            }
//
//            if ( loc < jobSchedule.size() )
//            {
//                // insert a new job package
//                std::vector<unsigned int> jobPackage;
//                for ( unsigned int jj=start; jj<=end; jj++ )
//                {
//                    jobPackage.push_back(jobIDs[jj]);
//                }
//
//                jobSchedule[loc].second.push_back(jobPackage);
//            }
//            else
//            {
//                // create a new entry
//                std::pair<unsigned int, std::vector<std::vector<unsigned int> > > jobItem;
//                jobItem.first = deviceID;
//
//                std::vector<unsigned int> jobPackage;
//                for ( unsigned int jj=start; jj<=end; jj++ )
//                {
//                    jobPackage.push_back(jobIDs[jj]);
//                }
//                jobItem.second.push_back(jobPackage);
//
//                jobSchedule.push_back(jobItem);
//            }
//
//            job = end+1;
//            validDevice++;
//
//            if ( validDevice >= validDevices.size() )
//            {
//                validDevice = 0;
//            }
//        }
//    }
//    catch(...)
//    {
//        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtil<T>::cudaJobSplitter(...) ... ");
//        return false;
//    }
//
//    return true;
//}
//
//template <typename T> 
//bool gtPlusISMRMRDReconUtil<T>::
//cudaJobSplitter(unsigned int numOfJobs, size_t jobSize, size_t minimalMemoryForValidDevice, 
//            std::vector< std::pair<unsigned int, std::vector<std::vector<unsigned int> > > >& jobSchedule)
//{
//    if ( numOfJobs == 0 )
//    {
//        GWARN_STREAM("numOfJobs == 0");
//        return true;
//    }
//
//    std::vector<unsigned int> jobIDs(numOfJobs, 0);
//    unsigned int ii;
//    for ( ii=0; ii<numOfJobs; ii++ ) jobIDs[ii] = ii;
//    return cudaJobSplitter(jobIDs, jobSize, minimalMemoryForValidDevice, jobSchedule);
//}
//
//#endif // USE_CUDA

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::
compareAgainstGroundTruthArray(const std::string& gt_filename, const hoNDArray<T>& x, typename realType<T>::Type& normDiff, typename realType<T>::Type& maxNormDiff)
{
    hoNDArray<T> gt;

    gtPlusIOAnalyze gt_io;
    gt_io.importArray(gt, gt_filename);

    compareAgainstGroundTruthArray(gt, x, normDiff, maxNormDiff);
}

template <typename T> 
void gtPlusISMRMRDReconUtil<T>::
compareAgainstGroundTruthArray(const hoNDArray<T>& gt, const hoNDArray<T>& x, typename realType<T>::Type& normDiff, typename realType<T>::Type& maxNormDiff)
{
    hoNDArray<T> diff(x);
    Gadgetron::subtract(gt, x, diff);

    typename realType<T>::Type v;
    Gadgetron::norm2(diff, v);
    normDiff = v;

    T maxV;
    size_t ind;
    Gadgetron::maxAbsolute(diff, maxV, ind);
    maxNormDiff = std::abs(maxV);
}

// ========================================================================================== //

template <typename T> 
gtPlusISMRMRDReconUtilComplex<T>::gtPlusISMRMRDReconUtilComplex() {}

template <typename T> 
gtPlusISMRMRDReconUtilComplex<T>::~gtPlusISMRMRDReconUtilComplex() {}

template <typename T> 
void gtPlusISMRMRDReconUtilComplex<T>::printInfo(std::ostream& os)
{
    using namespace std;
    os << "-------------- GTPlus ISMRMRD Recon Util Complex -------------" << endl;
    os << "Implementation of recon utilities for ISMRMRD complex data type" << endl;
    os << "--------------------------------------------------------------" << endl;
}

// ------------------------------------------------------------------------
// noise prewhitening
// ------------------------------------------------------------------------
template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
computeNoisePrewhiteningMatrix(const hoNDArray<T>& noise, double noiseBandWidth, double receiverBWRatio, double ADCSamplingTimeinSecond, hoMatrix<T>& prewhiteningMatrix)
{
    try
    {
        size_t RO = noise.get_size(0);
        size_t E1 = noise.get_size(1);
        size_t CHA = noise.get_size(2);

        GADGET_CHECK_RETURN_FALSE(prewhiteningMatrix.createMatrix(CHA, CHA));
        Gadgetron::clear(prewhiteningMatrix);

        typedef typename realType<T>::Type ValueType;

        // noise sampling time in second
        ValueType noiseSamplingTimeinSecond = (ValueType)(1.0/(noiseBandWidth*RO));

        // scaling factor
        ValueType scaling = (ValueType)(noiseSamplingTimeinSecond/ADCSamplingTimeinSecond/receiverBWRatio);
        scaling /= (RO*E1-1);

        // compute the noise covariance matrix
        hoMatrix<T> R(RO*E1, CHA, const_cast<T*>(noise.begin()));

        // R'*R --> CHA by CHA covariance matrix
        Gadgetron::gemm(prewhiteningMatrix, R, true, R, false);
        Gadgetron::scal(scaling, prewhiteningMatrix);

        // 0.5*(R+R')
        hoMatrix<T> RH(prewhiteningMatrix);
        conjugatetrans(prewhiteningMatrix, RH);
        Gadgetron::add(prewhiteningMatrix, RH, prewhiteningMatrix);
        Gadgetron::scal( (ValueType)0.5, prewhiteningMatrix);

        Gadgetron::potrf(prewhiteningMatrix, 'U');
        Gadgetron::trtri(prewhiteningMatrix, 'U');
        Gadgetron::scal( (value_type)(std::sqrt((double)2.0)), prewhiteningMatrix);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::computeNoisePrewhiteningMatrix(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
performNoisePrewhitening(hoNDArray<T>& data, const hoMatrix<T>& prewhiteningMatrix)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t CHA = data.get_size(2);

        GADGET_CHECK_RETURN_FALSE(prewhiteningMatrix.rows()==CHA);
        GADGET_CHECK_RETURN_FALSE(prewhiteningMatrix.cols()==CHA);

        size_t N = data.get_number_of_elements()/(RO*E1*CHA);

        long long n;
        #pragma omp parallel default(none) private(n) shared(RO, E1, CHA, N, data, prewhiteningMatrix)
        {
            hoMatrix<T> tmp(RO*E1, CHA);
            Gadgetron::clear(tmp);

            #pragma omp for
            for ( n=0; n<(long long)N; n++ )
            {
                hoMatrix<T> D(RO*E1, CHA, data.begin()+n*RO*E1*CHA);
                Gadgetron::gemm(tmp, D, false, prewhiteningMatrix, false);
                memcpy(data.begin()+n*RO*E1*CHA, tmp.begin(), sizeof(T)*RO*E1*CHA);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::performNoisePrewhitening(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);

        if ( RO==sizeX && E1==sizeY )
        {
            dataResized = data;
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY )
        {
            dataResized.create(sizeX, sizeY);
        }

        Gadgetron::clear(&dataResized);

        hoNDArray<T> kspace(data);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(data, kspace);
        GADGET_CHECK_RETURN_FALSE(zpadResize2DOnKSpace(kspace, sizeX, sizeY, dataResized));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize2D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize2DOnKSpace(const hoNDArray<T>& kspace, size_t sizeX, size_t sizeY, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);

        if ( RO==sizeX && E1==sizeY )
        {
            dataResized = kspace;
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(dataResized);
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY )
        {
            dataResized.create(sizeX, sizeY);
        }

        Gadgetron::clear(&dataResized);

        // GADGET_CHECK_RETURN_FALSE(this->zeropad2D(kspace, sizeX, sizeY, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad2D(kspace, sizeX, sizeY, dataResized));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(dataResized);

        typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeX*sizeY)/std::sqrt((double)RO*E1));
        Gadgetron::scal(scaling, dataResized);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize2DOnKSpace(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);
        GADGET_CHECK_RETURN_FALSE(sizeZ>=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataResized = data;
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY || dataResized.get_size(2)!=sizeZ )
        {
            dataResized.create(sizeX, sizeY, sizeZ);
        }

        Gadgetron::clear(&dataResized);

        hoNDArray<T> kspace(data);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, kspace);
        GADGET_CHECK_RETURN_FALSE(zpadResize3DOnKSpace(kspace, sizeX, sizeY, sizeZ, dataResized));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize3DOnKSpace(const hoNDArray<T>& kspace, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);
        GADGET_CHECK_RETURN_FALSE(sizeZ>=E2);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataResized = kspace;
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataResized);
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY || dataResized.get_size(2)!=sizeZ )
        {
            dataResized.create(sizeX, sizeY, sizeZ);
        }

        Gadgetron::clear(&dataResized);

        // GADGET_CHECK_RETURN_FALSE(this->zeropad3D(kspace, sizeX, sizeY, sizeZ, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(kspace, sizeX, sizeY, sizeZ, dataResized));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataResized);

        typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeX*sizeY*sizeZ)/std::sqrt((double)RO*E1*E2));
        Gadgetron::scal(scaling, dataResized);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize3DOnKSpace(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize2DFilter(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, const hoNDArray<T>& filter2D, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);

        GADGET_CHECK_RETURN_FALSE(filter2D.get_size(0)==sizeX);
        GADGET_CHECK_RETURN_FALSE(filter2D.get_size(1)==sizeY);

        if ( RO==sizeX && E1==sizeY )
        {
            dataResized = data;
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY )
        {
            dataResized.create(sizeX, sizeY);
        }

        Gadgetron::clear(&dataResized);

        hoNDArray<T> kspace(data);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(data, kspace);
        // GADGET_CHECK_RETURN_FALSE(this->zeropad2D(kspace, sizeX, sizeY, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad2D(kspace, sizeX, sizeY, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(dataResized, filter2D, dataResized));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(dataResized);

        typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeX*sizeY)/std::sqrt((double)RO*E1));
        Gadgetron::scal(scaling, dataResized);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize2DFilter(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
zpadResize3DFilter(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, const hoNDArray<T>& filter3D, hoNDArray<T>& dataResized)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);

        GADGET_CHECK_RETURN_FALSE(sizeX>=RO);
        GADGET_CHECK_RETURN_FALSE(sizeY>=E1);
        GADGET_CHECK_RETURN_FALSE(sizeZ>=E2);

        GADGET_CHECK_RETURN_FALSE(filter3D.get_size(0)==sizeX);
        GADGET_CHECK_RETURN_FALSE(filter3D.get_size(1)==sizeY);
        GADGET_CHECK_RETURN_FALSE(filter3D.get_size(2)==sizeZ);

        if ( RO==sizeX && E1==sizeY && E2==sizeZ )
        {
            dataResized = data;
            return true;
        }

        if ( dataResized.get_size(0)!=sizeX || dataResized.get_size(1)!=sizeY || dataResized.get_size(2)!=sizeZ )
        {
            dataResized.create(sizeX, sizeY, sizeZ);
        }

        Gadgetron::clear(&dataResized);

        hoNDArray<T> kspace(data);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, kspace);
        // GADGET_CHECK_RETURN_FALSE(this->zeropad3D(kspace, sizeX, sizeY, sizeZ, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(kspace, sizeX, sizeY, sizeZ, dataResized));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(dataResized, filter3D, dataResized));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataResized);

        typename realType<T>::Type scaling = (typename realType<T>::Type)(std::sqrt((double)sizeX*sizeY*sizeZ)/std::sqrt((double)RO*E1*E2));
        Gadgetron::scal(scaling, dataResized);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::zpadResize3DFilter(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterROImage(hoNDArray<T>& data, const hoNDArray<T>& fRO)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(data);
        GADGET_CHECK_RETURN_FALSE(this->kspacefilterRO(data, fRO));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(data);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterROImage(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterROImage(const hoNDArray<T>& data, const hoNDArray<T>& fRO, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(data, dataFiltered);
        GADGET_CHECK_RETURN_FALSE(this->kspacefilterRO(dataFiltered, fRO));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterROImage(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterE1Image(const hoNDArray<T>& data, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(data, dataFiltered);
        GADGET_CHECK_RETURN_FALSE(this->kspacefilterRO(dataFiltered, fE1, dataFiltered));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterE1Image(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterE2Image(const hoNDArray<T>& data, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, dataFiltered);
        GADGET_CHECK_RETURN_FALSE(this->kspacefilterRO(dataFiltered, fE2, dataFiltered));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterE2Image(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterE1E2Image(const hoNDArray<T>& data, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, dataFiltered);
        GADGET_CHECK_RETURN_FALSE(this->kspacefilterE1E2(dataFiltered, fE1, fE2, dataFiltered));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterE1E2Image(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
kspacefilterROE1E2Image(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==fRO.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==fE1.get_number_of_elements());
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==fE2.get_number_of_elements());

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(data, dataFiltered);
        GADGET_CHECK_RETURN_FALSE(this->kspacefilterROE1E2(dataFiltered, fRO, fE1, fE2, dataFiltered));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(dataFiltered);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::kspacefilterROE1E2Image(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilMap2DNIHInner(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t power)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long CHA = data.get_size(2);

        long long N = data.get_number_of_elements()/(RO*E1*CHA);
        GADGET_CHECK_RETURN_FALSE(N==1);

        const T* pData = data.begin();

        if ( !data.dimensions_equal(&coilMap) )
        {
            coilMap = data;
        }
        T* pSen = coilMap.begin();

        if ( ks%2 != 1 )
        {
            ks++;
        }

        size_t kss = ks*ks;
        long long halfKs = (long long)ks/2;

        int e1;

        // #pragma omp parallel default(none) private(e1) shared(ks, RO, E1, CHA, pSen, pData, halfKs, power, kss)
        #pragma omp parallel private(e1) shared(ks, RO, E1, CHA, pSen, pData, halfKs, power, kss)
        {
            hoNDArray<T> D(ks*ks, CHA);
            T* pD = D.begin();

            hoNDArray<T> DC(ks*ks, CHA);
            T* pDC = DC.begin();

            hoNDArray<T> DH_D(CHA, CHA);
            Gadgetron::clear(DH_D);

            hoNDArray<T> U1(ks*ks, 1);
            T* pU1 = U1.begin();

            hoNDArray<T> V1(CHA, 1);
            T* pV1 = V1.begin();

            hoNDArray<T> V(CHA, 1);

            Gadgetron::clear(D);
            Gadgetron::clear(DC);
            Gadgetron::clear(DH_D);
            Gadgetron::clear(U1);
            Gadgetron::clear(V1);
            Gadgetron::clear(V);

            T phaseU1;

            value_type v1Norm(1), u1Norm(1);

            long long cha, ro, kro, ke1, de1, dro;
            size_t po;

            #pragma omp for
            for ( e1=0; e1<(int)E1; e1++ )
            {
                for ( ro=0; ro<(long long)RO; ro++ )
                {
                    // fill the data matrix D
                    if ( e1>=halfKs && e1<E1-halfKs && ro>=halfKs && ro<RO-halfKs )
                    {
                        for ( cha=0; cha<CHA; cha++ )
                        {
                            const T* pDataCurr = pData + cha*RO*E1;
                            int ind=0;
                            for ( ke1=-halfKs; ke1<=halfKs; ke1++ )
                            {
                                de1 = e1 + ke1;
                                for ( kro=-halfKs; kro<=halfKs; kro++ )
                                {
                                    // D(ind++, cha) = pDataCurr[de1*RO+ro+kro];
                                    pD[ind+cha*kss] = pDataCurr[de1*RO+ro+kro];
                                    ind++;
                                }
                            }
                        }
                    }
                    else
                    {
                        for ( cha=0; cha<CHA; cha++ )
                        {
                            const T* pDataCurr = pData + cha*RO*E1;
                            int ind=0;
                            for ( ke1=-halfKs; ke1<=halfKs; ke1++ )
                            {
                                de1 = e1 + ke1;
                                if ( de1 < 0 ) de1 += E1;
                                if ( de1 >= E1 ) de1 -= E1;

                                for ( kro=-halfKs; kro<=halfKs; kro++ )
                                {
                                    dro = ro + kro;
                                    if ( dro < 0 ) dro += RO;
                                    if ( dro >= RO ) dro -= RO;

                                    // D(ind++, cha) = pDataCurr[de1*RO+dro];
                                    pD[ind+cha*kss] = pDataCurr[de1*RO+dro];
                                    ind++;
                                }
                            }
                        }
                    }

                    // compute V1
                    // D.sumOverCol(V1);
                    T* pTmp;
                    for ( cha=0; cha<CHA; cha++ )
                    {
                        pTmp = pD + cha*kss;
                        pV1[cha] = pTmp[0];
                        for ( po=1; po<kss; po++ )
                        {
                            pV1[cha] += pTmp[po];
                        }
                    }

                    // norm2(V1, v1Norm);
                    // Gadgetron::math::norm2(CHA, V1.begin(), v1Norm);
                    value_type sum(0);
                    for ( cha=0; cha<CHA; cha++ )
                    {
                        const T& c = pV1[cha];
                        const value_type re = c.real();
                        const value_type im = c.imag();
                        sum += ( (re*re) + (im * im) );
                    }
                    v1Norm = std::sqrt(sum);

                    // scal( (value_type)1.0/v1Norm, V1);
                    value_type v1NormInv = (value_type)1.0/v1Norm;
                    for ( cha=0; cha<CHA; cha++ )
                    {
                        pV1[cha] *= v1NormInv;
                    }

                    memcpy(pDC, pD, sizeof(T)*ks*ks*CHA);
                    gemm(DH_D, DC, true, D, false);

                    for ( po=0; po<power; po++ )
                    {
                        gemm(V, DH_D, false, V1, false);
                        // V1 = V;
                        memcpy(V1.begin(), V.begin(), V.get_number_of_bytes());

                        // norm2(V1, v1Norm);

                        sum = 0;
                        for ( cha=0; cha<CHA; cha++ )
                        {
                            const T& c = pV1[cha];
                            const value_type re = c.real();
                            const value_type im = c.imag();
                            sum += ( (re*re) + (im * im) );
                        }
                        v1Norm = std::sqrt(sum);

                        // scal( (value_type)1.0/v1Norm, V1);

                        value_type v1NormInv = (value_type)1.0/v1Norm;
                        for ( cha=0; cha<CHA; cha++ )
                        {
                            pV1[cha] *= v1NormInv;
                        }
                    }

                    // compute U1
                    gemm(U1, D, false, V1, false);

                    //phaseU1 = U1(0, 0);
                    phaseU1 = pU1[0];
                    for ( po=1; po<kss; po++ )
                    {
                        //phaseU1 += U1(po, 0);
                        phaseU1 += pU1[po];
                    }
                    phaseU1 /= std::abs(phaseU1);

                    // put the mean object phase to coil map
                    // conjugate(V1, V1);
                    // scal(phaseU1, V1);

                    const value_type c = phaseU1.real();
                    const value_type d = phaseU1.imag();

                    for ( cha=0; cha<CHA; cha++ )
                    {
                        const T& v = pV1[cha];
                        const value_type a = v.real();
                        const value_type b = v.imag();

                        reinterpret_cast< value_type(&)[2] >(pV1[cha])[0] = a*c+b*d;
                        reinterpret_cast< value_type(&)[2] >(pV1[cha])[1] = a*d-b*c;
                    }

                    for ( cha=0; cha<CHA; cha++ )
                    {
                        // pSen[cha*RO*E1+e1*RO+ro] = pV1[cha];
                        pSen[cha*RO*E1+e1*RO+ro] = V1(cha, 0);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilMap2DNIHInner(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilMap3DNIHInner(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t power)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long E2 = data.get_size(2);
        long long CHA = data.get_size(3);

        long long N = data.get_number_of_elements()/(RO*E1*E2*CHA);
        GADGET_CHECK_RETURN_FALSE(N==1);

        const T* pData = data.begin();

        if ( !data.dimensions_equal(&coilMap) )
        {
            coilMap = data;
        }
        T* pSen = coilMap.begin();

        if ( ks%2 != 1 )
        {
            ks++;
        }

        size_t kss = ks*ks*ks;
        long long halfKs = (long long)ks/2;

        long long e2;

        #pragma omp parallel default(none) private(e2) shared(ks, RO, E1, E2, CHA, pSen, pData, halfKs, power, kss)
        {
            hoMatrix<T> D(kss, CHA);
            hoMatrix<T> DC(kss, CHA);
            hoMatrix<T> DH_D(CHA, CHA);

            hoMatrix<T> U1(kss, 1);
            hoMatrix<T> V1(CHA, 1);
            hoMatrix<T> V(CHA, 1);

            Gadgetron::clear(D);
            Gadgetron::clear(DC);
            Gadgetron::clear(DH_D);
            Gadgetron::clear(U1);
            Gadgetron::clear(V1);
            Gadgetron::clear(V);

            T phaseU1;

            value_type v1Norm(1);

            long long cha, ro, e1, kro, dro, ke1, de1, ke2, de2;
            size_t po;

            #pragma omp for
            for ( e2=0; e2<(long long)E2; e2++ )
            {
                for ( e1=0; e1<(long long)E1; e1++ )
                {
                    for ( ro=0; ro<(long long)RO; ro++ )
                    {
                        // fill the data matrix D
                        if ( e2>=halfKs && e2<E2-halfKs && e1>=halfKs && e1<E1-halfKs && ro>=halfKs && ro<RO-halfKs )
                        {
                            for ( cha=0; cha<CHA; cha++ )
                            {
                                const T* pDataCurr = pData + cha*RO*E1*E2;
                                long long ind=0;
                                for ( ke2=-halfKs; ke2<=halfKs; ke2++ )
                                {
                                    de2 = e2 + ke2;
                                    for ( ke1=-halfKs; ke1<=halfKs; ke1++ )
                                    {
                                        de1 = e1 + ke1;
                                        for ( kro=-halfKs; kro<=halfKs; kro++ )
                                        {
                                            D(ind++, cha) = pDataCurr[de2*RO*E1+de1*RO+ro+kro];
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            for ( cha=0; cha<CHA; cha++ )
                            {
                                const T* pDataCurr = pData + cha*RO*E1*E2;
                                long long ind=0;
                                for ( ke2=-halfKs; ke2<=halfKs; ke2++ )
                                {
                                    de2 = e2 + ke2;
                                    if ( de2 < 0 ) de2 += E2;
                                    if ( de2 >= E2 ) de2 -= E2;

                                    for ( ke1=-halfKs; ke1<=halfKs; ke1++ )
                                    {
                                        de1 = e1 + ke1;
                                        if ( de1 < 0 ) de1 += E1;
                                        if ( de1 >= E1 ) de1 -= E1;

                                        for ( kro=-halfKs; kro<=halfKs; kro++ )
                                        {
                                            dro = ro + kro;
                                            if ( dro < 0 ) dro += RO;
                                            if ( dro >= RO ) dro -= RO;

                                            D(ind++, cha) = pDataCurr[de2*RO*E1+de1*RO+dro];
                                        }
                                    }
                                }
                            }
                        }

                        // compute V1
                        D.sumOverCol(V1);
                        norm2(V1, v1Norm);
                        scal( (value_type)1.0/v1Norm, V1);

                        memcpy(DC.begin(), D.begin(), sizeof(T)*kss*CHA);
                        gemm(DH_D, DC, true, D, false);
                        // gemm(DH_D, D, true, D, false);

                        for ( po=0; po<power; po++ )
                        {
                            gemm(V, DH_D, false, V1, false);
                            V1 = V;
                            norm2(V1, v1Norm);
                            scal( (value_type)1.0/v1Norm, V1);
                        }

                        // compute U1
                        gemm(U1, D, false, V1, false);

                        phaseU1 = U1(0, 0);
                        for ( po=1; po<kss; po++ )
                        {
                            phaseU1 += U1(po, 0);
                        }
                        phaseU1 /= std::abs(phaseU1);

                        // put the mean object phase to coil map
                        conjugate(V1, V1);
                        scal(phaseU1, V1);

                        for ( cha=0; cha<CHA; cha++ )
                        {
                            pSen[cha*RO*E1*E2+e2*RO*E1+e1*RO+ro] = V1(cha, 0);
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilMap3DNIHInner(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilMap2DNIH2Inner(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t iterNum, typename realType<T>::Type thres)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long CHA = data.get_size(2);

        long long N = data.get_number_of_elements()/(RO*E1*CHA);
        GADGET_CHECK_RETURN_FALSE(N==1);

        const T* pData = data.begin();

        if ( !data.dimensions_equal(&coilMap) )
        {
            coilMap = data;
        }

        // create convolution kernel
        hoNDArray<T> ker(ks, ks);
        Gadgetron::fill( ker, T( (value_type)1.0/(ks*ks)) );

        hoNDArray<T> prevR(RO, E1, 1), R(RO, E1, 1), imT(RO, E1, 1), magT(RO, E1, 1), diffR(RO, E1, 1);
        hoNDArray<T> coilMapConv(RO, E1, CHA);
        hoNDArray<T> D(RO, E1, CHA);
        hoNDArray<T> D_sum(1, E1, CHA);
        hoNDArray<T> D_sum_1st_2nd(1, 1, CHA);
        typename realType<T>::Type v, vR, vDiffR;
        T vCha;
        size_t iter;
        long long cha;

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, D_sum, 0));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D_sum, D_sum_1st_2nd, 1));
        Gadgetron::norm2(D_sum_1st_2nd, v);
        Gadgetron::scal( (value_type)1.0/v, D_sum_1st_2nd);

        Gadgetron::clear(R);
        for ( cha=0; cha<CHA; cha++ )
        {
            hoNDArray<T> dataCHA(RO, E1, const_cast<T*>(data.begin())+cha*RO*E1);
            vCha = D_sum_1st_2nd(cha);
            Gadgetron::axpy( std::conj(vCha), dataCHA, R, R);
        }

        for ( iter=0; iter<iterNum; iter++ )
        {
            prevR = R;

            Gadgetron::conjugate(R, R);

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, R, coilMap));

            Gadgetron::conv2(coilMap, ker, coilMapConv);

            Gadgetron::multiplyConj(coilMapConv, coilMapConv, D);

            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, R, 2));

            Gadgetron::sqrt(R, R);

            Gadgetron::addEpsilon(R);
            Gadgetron::inv(R, R);

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(coilMapConv, R, coilMap));

            Gadgetron::multiplyConj(data, coilMap, D);
            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, R, 2));

            //if ( iter < iterNum - 1 )
            //{
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(coilMap, R, D));
            //}
            //else
            //{
            //    D = coilMap;
            //}

            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, D_sum, 0));
            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D_sum, D_sum_1st_2nd, 1));

            Gadgetron::norm2(D_sum_1st_2nd, v);
            Gadgetron::scal( (value_type)1.0/v, D_sum_1st_2nd);

            Gadgetron::clear(imT);
            for ( cha=0; cha<CHA; cha++ )
            {
                hoNDArray<T> coilMapCHA(RO, E1, coilMap.begin()+cha*RO*E1);
                vCha = D_sum_1st_2nd(cha);
                Gadgetron::axpy( std::conj(vCha), coilMapCHA, imT, imT);
            }

            Gadgetron::abs(imT, magT);
            Gadgetron::divide(imT, magT, imT);

            Gadgetron::multiply(R, imT, R);
            Gadgetron::conjugate(imT, imT);
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(coilMap, imT, coilMap));

            Gadgetron::subtract(prevR, R, diffR);
            Gadgetron::norm2(diffR, vDiffR);
            Gadgetron::norm2(R, vR);

            // GDEBUG_STREAM("coilMap2DNIH2Inner - iter : " << iter << " - norm(prevR-R)/norm(R) : " << vDiffR/vR);

            if ( vDiffR/vR < thres ) break;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilMap2DNIH2Inner(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilMap2DNIH(const hoNDArray<T>& data, hoNDArray<T>& coilMap, ISMRMRDCOILMAPALGO algo, size_t ks, size_t power, size_t iterNum, typename realType<T>::Type thres, bool useGPU)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long CHA = data.get_size(2);

        size_t N = data.get_number_of_elements()/(RO*E1*CHA);
        size_t num = RO*E1*CHA;

        if ( !data.dimensions_equal(&coilMap) )
        {
            coilMap = data;
        }

        if ( ks%2 != 1 )
        {
            ks++;
        }

        long long n;

        if ( N >= 16 )
        {
            #pragma omp parallel default(none) private(n) shared(ks, RO, E1, CHA, num, algo, N, data, coilMap, power, iterNum, thres)
            {
                #pragma omp for
                for ( n=0; n<(long long)N; n++ )
                {
                    hoNDArray<T> dataCurr(RO, E1, CHA, const_cast<T*>(data.begin()+n*num));
                    hoNDArray<T> coilMapCurr(RO, E1, CHA, coilMap.begin()+n*num);

                    if ( algo == ISMRMRD_SOUHEIL_ITER )
                    {
                        coilMap2DNIH2Inner(dataCurr, coilMapCurr, ks, iterNum, thres);
                    }
                    else
                    {
                        coilMap2DNIHInner(dataCurr, coilMapCurr, ks, power);
                    }
                }
            }
        }
        else if ( N == 1 )
        {
            if ( algo == ISMRMRD_SOUHEIL_ITER )
            {
                GADGET_CHECK_RETURN_FALSE(coilMap2DNIH2Inner(data, coilMap, ks, iterNum, thres));
            }
            else
            {
                GADGET_CHECK_RETURN_FALSE(coilMap2DNIHInner(data, coilMap, ks, power));
            }
        }
        else
        {
            for ( n=0; n<(long long)N; n++ )
            {
                hoNDArray<T> dataCurr(RO, E1, CHA, const_cast<T*>(data.begin()+n*num));
                hoNDArray<T> coilMapCurr(RO, E1, CHA, coilMap.begin()+n*num);
                if ( algo == ISMRMRD_SOUHEIL_ITER )
                {
                    GADGET_CHECK_RETURN_FALSE(coilMap2DNIH2Inner(dataCurr, coilMapCurr, ks, iterNum, thres));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(coilMap2DNIHInner(dataCurr, coilMapCurr, ks, power));
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilMap2DNIH(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilMap3DNIH2Inner(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t kz, size_t iterNum, typename realType<T>::Type thres)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t CHA = data.get_size(3);

        size_t N = data.get_number_of_elements()/(RO*E1*E2*CHA);
        GADGET_CHECK_RETURN_FALSE(N==1);

        const T* pData = data.begin();

        if ( !data.dimensions_equal(&coilMap) )
        {
            coilMap = data;
        }

        // create convolution kernel
        hoNDArray<T> ker(ks, ks, kz);
        Gadgetron::fill( &ker, T( (value_type)1.0/(ks*ks*kz)) );

        hoNDArray<T> R(RO, E1, E2, 1), imT(RO, E1, E2, 1), magT(RO, E1, E2, 1);
        hoNDArray<T> coilMapConv(RO, E1, E2, CHA);
        hoNDArray<T> D(RO, E1, E2, CHA);
        hoNDArray<T> D_sum(1, CHA);
        typename realType<T>::Type v;
        T vCha;
        size_t iter, cha;

        hoNDArray<T> dataByCha(RO*E1*E2, CHA, const_cast<T*>(data.begin()));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, D_sum, 0));
        Gadgetron::norm2(D_sum, v);
        Gadgetron::scal( (value_type)1.0/v, D_sum);

        Gadgetron::clear(R);
        for ( cha=0; cha<CHA; cha++ )
        {
            hoNDArray<T> dataCHA(RO, E1, E2, const_cast<T*>(data.begin())+cha*RO*E1*E2);
            vCha = D_sum(cha);
            Gadgetron::axpy( std::conj(vCha), dataCHA, R, R);
        }

        for ( iter=0; iter<iterNum; iter++ )
        {
            Gadgetron::conjugate(R, R);

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(data, R, coilMap));

            Gadgetron::conv2(coilMap, ker, coilMapConv);

            Gadgetron::multiplyConj(coilMapConv, coilMapConv, D);

            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, R, 3));

            Gadgetron::sqrt(R, R);

            Gadgetron::addEpsilon(R);
            Gadgetron::inv(R, R);

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(coilMapConv, R, coilMap));

            Gadgetron::multiplyConj(data, coilMap, D);
            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, R, 3));

            //if ( iter < iterNum - 1 )
            //{
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(coilMap, R, D));
            //}
            //else
            //{
            //    D = coilMap;
            //}

            hoNDArray<T> DByCha(RO*E1*E2, CHA, D.begin());
            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(DByCha, D_sum, 0));

            Gadgetron::norm2(D_sum, v);
            Gadgetron::scal( (value_type)1.0/v, D_sum);

            Gadgetron::clear(imT);
            for ( cha=0; cha<CHA; cha++ )
            {
                hoNDArray<T> coilMapCHA(RO, E1, E2, 1, coilMap.begin()+cha*RO*E1*E2);
                vCha = D_sum(cha);
                Gadgetron::axpy( std::conj(vCha), coilMapCHA, imT, imT);
            }

            Gadgetron::abs(imT, magT);
            Gadgetron::divide(imT, magT, imT);

            Gadgetron::multiply(R, imT, R);
            Gadgetron::conjugate(imT, imT);
            Gadgetron::multiply(coilMap, imT, coilMap);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilMap3DNIH2Inner(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilMap3DNIH(const hoNDArray<T>& data, hoNDArray<T>& coilMap, ISMRMRDCOILMAPALGO algo, size_t ks, size_t power, size_t iterNum, typename realType<T>::Type thres, bool true3D)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t CHA = data.get_size(3);

        size_t N = data.get_number_of_elements()/(RO*E1*E2*CHA);

        if ( !data.dimensions_equal(&coilMap) )
        {
            coilMap = data;
        }

        if ( ks%2 != 1 )
        {
            ks++;
        }

        //std::string debugFolder = "D:/software/Gadgetron/20130114/install/gadgetron/DebugOutput/";
        //gtPlusIOAnalyze gt_io;

        hoNDArray<T> data2D, coilMap2D;

        if ( algo == ISMRMRD_SOUHEIL )
        {
            data2D.create(RO, E1, CHA);
            coilMap2D.create(RO, E1, CHA);
        }

        int n, e2;
        for ( n=0; n<(long long)N; n++ )
        {
            if ( algo == ISMRMRD_SOUHEIL_ITER )
            {
                GDEBUG_STREAM("calling 3D version of Souhiel iterative coil map estimation ... ");
                GADGET_CHECK_RETURN_FALSE(this->coilMap3DNIH2Inner(data, coilMap, ks, ks, iterNum, thres));
            }
            else if ( algo==ISMRMRD_SOUHEIL && E2>5*ks && true3D )
            {
                GDEBUG_STREAM("calling 3D version of Souhiel coil map estimation ... ");
                GADGET_CHECK_RETURN_FALSE(this->coilMap3DNIHInner(data, coilMap, ks, power));
            }
            else
            {
                hoNDArray<T> dataCurr(RO, E1, E2, CHA, const_cast<T*>(data.begin()+n*RO*E1*E2*CHA));
                hoNDArray<T> coilMapCurr(RO, E1, E2, CHA, coilMap.begin()+n*RO*E1*E2*CHA);

                #pragma omp parallel default(none) private(e2) shared(dataCurr, coilMapCurr, RO, E1, E2, CHA, algo, ks, power, iterNum, thres) if (E2>12)
                {
                    hoNDArray<T> data2D(RO, E1, CHA);
                    hoNDArray<T> coilMap2D(RO, E1, CHA);

                    #pragma omp for
                    for ( e2=0; e2<(int)E2; e2++ )
                    {
                        long long cha;

                        for ( cha=0; cha<(long long)CHA; cha++ )
                        {
                            memcpy(data2D.begin()+cha*RO*E1, dataCurr.begin()+cha*RO*E1*E2+e2*RO*E1, sizeof(T)*RO*E1);
                        }

                        if ( algo == ISMRMRD_SOUHEIL_ITER )
                        {
                            coilMap2DNIH2Inner(data2D, coilMap2D, ks, iterNum, thres);
                        }
                        else
                        {
                            coilMap2DNIHInner(data2D, coilMap2D, ks, power);
                        }

                        for ( cha=0; cha<(long long)CHA; cha++ )
                        {
                            memcpy(coilMapCurr.begin()+cha*RO*E1*E2+e2*RO*E1, coilMap2D.begin()+cha*RO*E1, sizeof(T)*RO*E1);
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilMap3DNIH(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
sumOfSquare(const hoNDArray<T>& data, hoNDArray<T>& sos)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim>=3);

        hoNDArray<T> tmp(data);
        Gadgetron::multiplyConj(data, data, tmp);

        if ( NDim == 3 )
        {
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(tmp, sos));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(tmp, sos, 2));

            std::vector<size_t> dim(2);
            dim[0] = sos.get_size(0);
            dim[1] = sos.get_size(1);

            sos.reshape(dim);
        }
        else if ( NDim == 4 )
        {
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverSecondLastDimension(tmp, sos));

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(tmp, sos, 2));

            std::vector<size_t> dim(3);
            dim[0] = sos.get_size(0);
            dim[1] = sos.get_size(1);
            dim[2] = sos.get_size(3);

            sos.reshape(dim);
        }
        else
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(tmp, sos, 2));
        }

        Gadgetron::sqrt(sos, sos);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::sumOfSquare(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilCombine(const hoNDArray<T>& data, const hoNDArray<T>& coilMap, hoNDArray<T>& combined)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        size_t NDimCoil = coilMap.get_number_of_dimensions();

        // GADGET_CHECK_RETURN_FALSE(NDimCoil<=NDim);
        GADGET_CHECK_RETURN_FALSE(data.get_number_of_elements()>=coilMap.get_number_of_elements());

        size_t n;
        for ( n=0; n<NDimCoil; n++ )
        {
            if ( n<NDim && coilMap.get_size(n)>1 )
            {
                GADGET_CHECK_RETURN_FALSE(data.get_size(n)==coilMap.get_size(n));
            }
        }

        boost::shared_ptr< std::vector<size_t> > dim = data.get_dimensions();
        boost::shared_ptr< std::vector<size_t> > dimCoil = coilMap.get_dimensions();

        std::vector<size_t> dimCombined(*dim);
        dimCombined.erase(dimCombined.begin()+2);
        combined.create(&dimCombined);

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t CHA = data.get_size(2);
        size_t N = data.get_size(3);

        size_t coilN = coilMap.get_size(3);

        if ( coilN < N )
        {
            size_t NCombined = data.get_number_of_elements()/(RO*E1*CHA);

            std::vector<size_t> dataInd, coilMapInd(NDimCoil, 0), coimbinedInd(dimCombined.size(), 0);

            size_t nn;
            size_t d;
            hoNDArray<T> dataTmp(RO, E1, CHA);
            hoNDArray<T> combinedCurr(RO, E1, 1);

            for ( nn=0; nn<NCombined; nn++ )
            {
                size_t offsetData = nn*RO*E1*CHA;
                dataInd = data.calculate_index(offsetData);

                for ( d=0; d<NDimCoil; d++ )
                {
                    if ( dataInd[d]<coilMap.get_size(d) )
                    {
                        coilMapInd[d] = dataInd[d];
                    }
                    else
                    {
                        coilMapInd[d] = 0;
                    }
                }

                for ( d=3; d<NDim; d++ )
                {
                    coimbinedInd[d-1] = dataInd[d];
                }

                size_t offsetCoilMap = coilMap.calculate_offset(coilMapInd);
                size_t offsetCombined = combined.calculate_offset(coimbinedInd);

                hoNDArray<T> dataCurr(RO, E1, CHA, const_cast<T*>(data.begin())+offsetData);
                hoNDArray<T> coilMapCurr(RO, E1, CHA, const_cast<T*>(coilMap.begin())+offsetCoilMap);

                Gadgetron::multiplyConj(dataCurr, coilMapCurr, dataTmp);
                GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(dataTmp, combinedCurr, 2));

                memcpy(combined.begin()+offsetCombined, combinedCurr.begin(), sizeof(T)*RO*E1);
            }
        }
        else
        {
            size_t NCombined = data.get_number_of_elements()/(RO*E1*CHA*N);

            std::vector<size_t> dataInd, coilMapInd(NDimCoil, 0), coimbinedInd(dimCombined.size(), 0);

            size_t nn;
            size_t d;
            hoNDArray<T> dataTmp(RO, E1, CHA, N);
            hoNDArray<T> combinedCurr(RO, E1, 1, N);

            for ( nn=0; nn<NCombined; nn++ )
            {
                size_t offsetData = nn*RO*E1*CHA*N;
                dataInd = data.calculate_index(offsetData);

                for ( d=0; d<NDimCoil; d++ )
                {
                    if ( dataInd[d]<coilMap.get_size(d) )
                    {
                        coilMapInd[d] = dataInd[d];
                    }
                    else
                    {
                        coilMapInd[d] = 0;
                    }
                }

                for ( d=3; d<NDim; d++ )
                {
                    coimbinedInd[d-1] = dataInd[d];
                }

                size_t offsetCoilMap = coilMap.calculate_offset(coilMapInd);
                size_t offsetCombined = combined.calculate_offset(coimbinedInd);

                hoNDArray<T> dataCurr(RO, E1, CHA, N, const_cast<T*>(data.begin())+offsetData);
                hoNDArray<T> coilMapCurr(RO, E1, CHA, N, const_cast<T*>(coilMap.begin())+offsetCoilMap);

                Gadgetron::multiplyConj(dataCurr, coilMapCurr, dataTmp);
                GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(dataTmp, combinedCurr, 2));

                memcpy(combined.begin()+offsetCombined, combinedCurr.begin(), sizeof(T)*RO*E1*N);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilCombine(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
coilCombine3D(const hoNDArray<T>& data, const hoNDArray<T>& coilMap, hoNDArray<T>& combined)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        size_t NDimCoil = coilMap.get_number_of_dimensions();

        // GADGET_CHECK_RETURN_FALSE(NDimCoil<=NDim);
        GADGET_CHECK_RETURN_FALSE(data.get_number_of_elements()>=coilMap.get_number_of_elements());

        /*size_t n;
        for ( n=0; n<NDimCoil; n++ )
        {
            GADGET_CHECK_RETURN_FALSE(data.get_size(n)==coilMap.get_size(n));
        }*/

        GADGET_CHECK_RETURN_FALSE(data.get_size(0)==coilMap.get_size(0));
        GADGET_CHECK_RETURN_FALSE(data.get_size(1)==coilMap.get_size(1));
        GADGET_CHECK_RETURN_FALSE(data.get_size(2)==coilMap.get_size(2));
        GADGET_CHECK_RETURN_FALSE(data.get_size(3)==coilMap.get_size(3));

        boost::shared_ptr< std::vector<size_t> > dim = data.get_dimensions();
        boost::shared_ptr< std::vector<size_t> > dimCoil = coilMap.get_dimensions();

        size_t N = coilMap.get_number_of_elements();
        size_t num = data.get_number_of_elements()/coilMap.get_number_of_elements();

        std::vector<size_t> dimCombined(*dim);
        dimCombined.erase(dimCombined.begin()+3);
        combined.create(&dimCombined);

        std::vector<size_t> dimCombinedCurr(*dimCoil);
        dimCombinedCurr[3] = 1;

        size_t NCombined = combined.get_number_of_elements()/num;

        long long nn;
        #pragma omp parallel default(none) private(nn) shared(data, coilMap, num, dimCoil, dimCombinedCurr, combined, N, NCombined) if (num>=6)
        {
            hoNDArray<T> dataTmp(coilMap);

            #pragma omp for
            for ( nn=0; nn<(long long)num; nn++ )
            {
                hoNDArray<T> dataCurr(dimCoil.get(), const_cast<T*>(data.begin()+nn*N));
                Gadgetron::multiplyConj(dataCurr, coilMap, dataTmp);

                hoNDArray<T> dataCombinedCurr(&dimCombinedCurr, const_cast<T*>(combined.begin()+nn*NCombined));
                Gadgetron::sum_over_dimension(dataTmp, dataCombinedCurr, 3);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::coilCombine3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
conjugateSymmetry2D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj)
{
    try
    {
        if ( !kspaceConj.dimensions_equal(&kspace) )
        {
            kspaceConj.create(kspace.get_dimensions());
        }

        long long RO = kspace.get_size(0);
        long long E1 = kspace.get_size(1);
        long long num = kspace.get_number_of_elements()/(RO*E1);

        long long centerRO = RO/2;
        long long centerE1 = E1/2;

        long long ii;

        #pragma omp parallel for default(none) private(ii) shared(RO, E1, num, centerRO, centerE1, kspace, kspaceConj)
        for ( ii=0; ii<num; ii++ )
        {
            ho2DArray<T> src(RO, E1, const_cast<T*>(kspace.begin()+ii*RO*E1));
            ho2DArray<T> dst(RO, E1, const_cast<T*>(kspaceConj.begin()+ii*RO*E1));

            long long ro, e1;
            long long cro, ce1;

            for ( e1=0; e1<E1; e1++ )
            {
                ce1 = 2*centerE1-e1;
                if ( ce1 > E1-1 )
                {
                    ce1 -= E1;
                }

                for ( ro=0; ro<RO; ro++ )
                {
                    cro = 2*centerRO-ro;
                    if ( cro > RO-1 )
                    {
                        cro -= RO;
                    }

                    dst(ro, e1) = std::conj(src(cro, ce1));
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::conjugateSymmetry2D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconUtilComplex<T>::
conjugateSymmetry3D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj)
{
    try
    {
        if ( !kspaceConj.dimensions_equal(&kspace) )
        {
            kspaceConj.create(kspace.get_dimensions());
        }

        long long RO = kspace.get_size(0);
        long long E1 = kspace.get_size(1);
        long long E2 = kspace.get_size(2);
        long long num = kspace.get_number_of_elements()/(RO*E1*E2);

        long long centerRO = RO/2;
        long long centerE1 = E1/2;
        long long centerE2 = E2/2;

        long long ii;

        #pragma omp parallel for default(none) private(ii) shared(RO, E1, E2, num, centerRO, centerE1, centerE2, kspace, kspaceConj)
        for ( ii=0; ii<num; ii++ )
        {
            ho3DArray<T> src(RO, E1, E2, const_cast<T*>(kspace.begin()+ii*RO*E1*E2));
            ho3DArray<T> dst(RO, E1, E2, const_cast<T*>(kspaceConj.begin()+ii*RO*E1*E2));

            long long ro, e1, e2;
            long long cro, ce1, ce2;

            for ( e2=0; e2<E2; e2++ )
            {
                ce2 = 2*centerE2-e2;
                if ( ce2 > E2-1 )
                {
                    ce2 -= E2;
                }

                for ( e1=0; e1<E1; e1++ )
                {
                    ce1 = 2*centerE1-e1;
                    if ( ce1 > E1-1 )
                    {
                        ce1 -= E1;
                    }

                    for ( ro=0; ro<RO; ro++ )
                    {
                        cro = 2*centerRO-ro;
                        if ( cro > RO-1 )
                        {
                            cro -= RO;
                        }

                        dst(ro, e1, e2) = std::conj(src(cro, ce1, ce2));
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconUtilComplex<T>::conjugateSymmetry3D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj) ... ");
        return false;
    }

    return true;
}

}}
