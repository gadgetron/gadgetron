
/** \file   gtPlusGRAPPA.h
    \brief  GRAPPA implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue

    References to the implementation can be found in:

    Griswold MA, Jakob PM, Heidemann RM, Nittka M, Jellus V, Wang J, Kiefer B, Haase A. 
    Generalized autocalibrating partially parallel acquisitions (GRAPPA). 
    Magnetic Resonance in Medicine 2002;47(6):1202-1210.

    Kellman P, Epstein FH, McVeigh ER. 
    Adaptive sensitivity encoding incorporating temporal filtering (TSENSE). 
    Magnetic Resonance in Medicine 2001;45(5):846-852.

    Breuer FA, Kellman P, Griswold MA, Jakob PM. .
    Dynamic autocalibrated parallel imaging using temporal GRAPPA (TGRAPPA). 
    Magnetic Resonance in Medicine 2005;53(4):981-985.

    Saybasili H., Kellman P., Griswold MA., Derbyshire JA. Guttman, MA. 
    HTGRAPPA: Real-time B1-weighted image domain TGRAPPA reconstruction. 
    Magnetic Resonance in Medicine 2009;61(6): 1425-1433. 
*/

#pragma once

#include "gtPlusAlgorithmBase.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusGRAPPA : public gtPlusAlgorithmBase<T>
{
public:

    typedef gtPlusAlgorithmBase<T> BaseClass;

    gtPlusGRAPPA() : calib_use_gpu_(true), BaseClass() {}
    virtual ~gtPlusGRAPPA() {}

    virtual void printInfo(std::ostream& os);

    // get the kernel pattern, given the acceleration factor and kernel size
    bool kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, size_t accelFactor, size_t kNE1, bool fitItself);

    // grappa calibration for 2D case
    // acsSrc : [RO E1 srcCHA]
    // acsDst : [RO E1 dstCHA]
    // ker : [kRO kE1 srcCHA dstCHA oE1]
    bool calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, 
            size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray<T>& ker);

    // image domain kernel for 2D kernel
    // kIm: image domain kernel [RO E1 srcCHA dstCHA]
    bool imageDomainKernel(const ho5DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t ro, size_t e1, hoNDArray<T>& kIm);

    // grappa calibration for 3D case
    // acsSrc : [RO E1 E2 srcCHA]
    // acsDst : [RO E1 E2 dstCHA]
    // ker : [kRO kE1 kE2 srcCHA dstCHA oE1 oE2]
    bool calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, double thres, double overDetermineRatio, 
            size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, ho7DArray<T>& ker);

    // image domain kernel for 3D kernel
    // kIm: image domain kernel [RO E1 E2 srcCHA dstCHA]
    bool imageDomainKernel3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, size_t ro, size_t e1, size_t e2, hoNDArray<T>& kIm);

    // convert the calibrated kernel to the convlution kernel in kspace
    // if ROis3rdDim == true, the kernel dimension is [E1 E2 RO], otherwise [RO E1 E2]
    bool kspaceDomainConvKernel3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, ho5DArray<T>& convKerFlip, bool ROis3rdDim=true);

    // image domain kernel for 3D kernel, only RO direction is converted to image domain
    // E1 and E2 stays in the kspace domain
    // kImRO: kspace-image hybrid kernel [convE1 convE2 RO srcCHA dstCHA]
    bool imageDomainKernelRO3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, size_t ro, hoNDArray<T>& kImRO);

    // image domain kernel for 3D kernel, E1 and E2 directions are converted to image domain
    // kImRO : kspace-image hybrid kernel where first two dimensions are E1 and E2 and in kspace
    bool imageDomainKernelE1E2RO(const hoNDArray<T>& kImRO, size_t e1, size_t e2, hoNDArray<T>& kImE1E2RO);

    // use gpu in the kernel calibration
    bool calib_use_gpu_;

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;
    using BaseClass::gtPlus_mem_manager_;
};

template <typename T> 
void gtPlusGRAPPA<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD GRAPPA reconstruction ------------------" << endl;
    os << "Implementation of GRAPPA algorithms for ISMRMRD package" << endl;
    os << "Both 2D and 3D version are implemented" << endl;
    os << "Algorithms are published at:" << endl;
    os << "Generalized autocalibrating partially parallel acquisitions (GRAPPA), Magnetic Resonance in Medicine, Volume 47, Issue 6, pages 1202�1210, June 2002" << endl;
    os << "HTGRAPPA: Real-time B1-weighted image domain TGRAPPA reconstruction, Magnetic Resonance in Medicine, Volume 61, Issue 6, pages 1425�1433, June 2009" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, size_t accelFactor, size_t kNE1, bool fitItself)
{
    if ( accelFactor == 1 )
    {
        kE1.resize(1, 0);
        oE1.resize(1, 0);
        return true;
    }

    kE1.resize(kNE1, 0);
    if ( kNE1%2 == 0 )
    {
        long long k;
        for ( k=-((long long)kNE1/2-1); k<=(long long)kNE1/2; k++ )
        {
            kE1[k+kNE1/2-1] = (int)(k*accelFactor);
        }
    }
    else
    {
        long long k;
        for ( k=-(long long)kNE1/2; k<=(long long)kNE1/2; k++ )
        {
            kE1[k+kNE1/2] = (int)(k*accelFactor);
        }
    }

    if ( fitItself )
    {
        oE1.resize(accelFactor);
        for ( long long a=0; a<(long long)accelFactor; a++ )
        {
            oE1[a] = (int)a;
        }
    }
    else
    {
        oE1.resize(accelFactor-1);
        for ( long long a=1; a<(long long)accelFactor; a++ )
        {
            oE1[a-1] = (int)a;
        }
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, 
    size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray<T>& ker)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(0)==acsDst.get_size(0));
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(1)==acsDst.get_size(1));
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(2)>=acsDst.get_size(2));

        size_t RO = acsSrc.get_size(0);
        size_t E1 = acsSrc.get_size(1);
        size_t srcCHA = acsSrc.get_size(2);
        size_t dstCHA = acsDst.get_size(2);

        const T* pSrc = acsSrc.begin();
        const T* pDst = acsDst.begin();

        long long kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GWARN_STREAM("gtPlusGRAPPA<T>::calib(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        size_t kNE1 = kE1.size();
        size_t oNE1 = oE1.size();

        // allocate kernel
        GADGET_CHECK_RETURN_FALSE(ker.createArray(kRO, kNE1, srcCHA, dstCHA, oNE1));

        // loop over the calibration region and assemble the equation
        // Ax = b

        size_t eRO = RO - kROhalf -1;
        size_t sE1 = std::abs(kE1[0]);
        size_t eE1 = E1 -1 - kE1[kNE1-1];

        size_t lenRO = eRO-kROhalf+1;

        size_t rowA = (eE1-sE1+1)*lenRO;
        size_t colA = kRO*kNE1*srcCHA;
        size_t colB = dstCHA*oNE1;

        hoMatrix<T> A;
        hoMatrix<T> B;
        hoMatrix<T> x( colA, colB );

        if ( performTiming_ ) { gt_timer3_.start("grappa 2D calibration - allocate matrix storage ... "); }
        hoNDArrayMemoryManaged<T> A_mem(rowA, colA, gtPlus_mem_manager_);
        A.createMatrix( rowA, colA, A_mem.begin() );
        T* pA = A.begin();

        hoNDArrayMemoryManaged<T> B_mem(rowA, colB, gtPlus_mem_manager_);
        B.createMatrix( A.rows(), colB, B_mem.begin() );
        T* pB = B.begin();
        if ( performTiming_ ) { gt_timer3_.stop(); }

        long long e1;
        for ( e1=(long long)sE1; e1<=(long long)eE1; e1++ )
        {
            for ( long long ro=kROhalf; ro<=(long long)eRO; ro++ )
            {
                long long rInd = (e1-sE1)*lenRO+ro-kROhalf;

                size_t src, dst, ke1, oe1;
                long long kro;

                // fill matrix A
                size_t col = 0;
                size_t offset = 0;
                for ( src=0; src<srcCHA; src++ )
                {
                    for ( ke1=0; ke1<kNE1; ke1++ )
                    {
                        offset = src*RO*E1 + (e1+kE1[ke1])*RO;
                        for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                        {
                            // A(rInd, col++) = acsSrc(ro+kro, e1+kE1[ke1], src);
                            pA[rInd + col*rowA] = pSrc[ro+kro+offset];
                            col++;
                        }
                    }
                }

                // fill matrix B
                col = 0;
                for ( oe1=0; oe1<oNE1; oe1++ )
                {
                    for ( dst=0; dst<dstCHA; dst++ )
                    {
                        B(rInd, col++) = acsDst(ro, e1+oE1[oe1], dst);
                    }
                }
            }
        }

        //typename realType<T>::Type v;

        //Gadgetron::norm2(A, v);
        //GDEBUG_STREAM("A = " << v);

        //Gadgetron::norm2(B, v);
        //GDEBUG_STREAM("B = " << v);

        //if ( performTiming_ ) { gt_timer2_.start("SolveLinearSystem_Tikhonov"); }
        //#ifdef USE_CUDA
        //    // go to device
        //    try
        //    {
        //        if ( typeid(typename realType<T>::Type)==typeid(float) && calib_use_gpu_ )
        //        {
        //            GDEBUG_STREAM("grappa 2D - calling GPU kernel estimation ... ");
        //            hoNDArray<float_complext> A_tmp(A.get_dimensions(), reinterpret_cast<float_complext*>(A.begin()));
        //            hoNDArray<float_complext> B_tmp(B.get_dimensions(), reinterpret_cast<float_complext*>(B.begin()));

        //            int ret(0);
        //            boost::shared_ptr< hoNDArray<complext<float> > > host_x;

        //            #pragma omp critical(inverse)
        //            {
        //                cuNDArray<float_complext> device_A(A_tmp);
        //                cuNDArray<float_complext> device_B(B_tmp);
        //                cuNDArray<float_complext> device_x;

        //                ret = Gadgetron::inverse_clib_matrix(&device_A, &device_B, &device_x, thres);
        //                if ( ret == 0 )
        //                {
        //                    host_x = device_x.to_host();
        //                }
        //            }

        //            if ( ret != 0 )
        //            {
        //                GERROR_STREAM("failed in Gadgetron::inverse_clib_matrix(&device_A, &device_B, &device_x, thres) ... ");
        //                SolveLinearSystem_Tikhonov(A, B, x, thres);
        //            }
        //            else
        //            {
        //                memcpy(x.begin(), host_x->begin(), x.get_number_of_bytes());
        //            }
        //        }
        //        else
        //        {
        //            if ( calib_use_gpu_ )
        //            {
        //                GWARN_STREAM("GPU inverse_clib_matrix for grappa is only available for single-precision, calling the CPU version ... ");
        //            }

        //            GADGET_CHECK_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));
        //        }
        //    }
        //    catch(...)
        //    {
        //        GERROR_STREAM("failed in GPU inverse_clib_matrix for grappa, calling the CPU version ... ");
        //        GADGET_CHECK_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));
        //    }

        //#else
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));
        //#endif // USE_CUDA

        // GADGET_CHECK_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));
        //if ( performTiming_ ) { gt_timer2_.stop(); }

        //Gadgetron::norm2(x, v);
        //GDEBUG_STREAM("x = " << v);

        // the matrix dimension just matches
        // hoMatrix<T> xt(x.cols(), x.rows(), ker.begin());
        // GADGET_CHECK_RETURN_FALSE(Gadgetron::trans(x, xt));
        memcpy(ker.begin(), x.begin(), ker.get_number_of_bytes());

        //Gadgetron::norm2(ker, v);
        //GDEBUG_STREAM("ker = " << v);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusGRAPPA<T>::calib(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
imageDomainKernel(const ho5DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t ro, size_t e1, hoNDArray<T>& kIm)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(2));
        long long dstCHA = (long long)(ker.get_size(3));
        long long kNE1 = (long long)(kE1.size());
        long long oNE1 = (long long)(oE1.size());

        long long kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GWARN_STREAM("gtPlusGRAPPA<T>::imageDomainKernel(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        // allocate image domain kernel
        kIm.create(ro, e1, srcCHA, dstCHA);

        /// fill the convolution kernels
        long long convKRO = 2*kRO+3;

        long long maxKE1 = std::abs(kE1[0]);
        if ( std::abs(kE1[kNE1-1]) > maxKE1 )
        {
            maxKE1 = std::abs(kE1[kNE1-1]);
        }
        long long convKE1 = 2*maxKE1+1;

        /// allocate the convolution kernel
        hoNDArray<T> convKer(convKRO, convKE1, srcCHA, dstCHA);
        Gadgetron::clear(&convKer);

        /// index
        long long oe1, kro, ke1, src, dst;

        /// fill the convolution kernel and sum up multiple kernels
        for ( oe1=0; oe1<oNE1; oe1++ )
        {
            for ( ke1=0; ke1<kNE1; ke1++ )
            {
                for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                {
                    for ( dst=0; dst<dstCHA; dst++ )
                    {
                        for ( src=0; src<srcCHA; src++ )
                        {
                            convKer(-kro+kRO+1, oE1[oe1]-kE1[ke1]+maxKE1, src, dst) = ker(kro+kROhalf, ke1, src, dst, oe1);
                        }
                    }

                }
            }
        }

        if ( (oE1[0]!=0) && (srcCHA==dstCHA) )
        {
            for ( dst=0; dst<dstCHA; dst++ )
            {
                convKer(kRO+1, maxKE1, dst, dst) = 1.0;
            }
        }

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1)) ), convKer ));
        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad2D(convKer, ro, e1, kIm));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad2D(convKer, ro, e1, kIm));

        //uint64d<2>::Type sizeArray;
        //sizeArray[0] = ro;
        //sizeArray[1] = e1;
        //Gadgetron::pad<T, 2>(sizeArray, &convKer, &kIm);

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kIm);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusGRAPPA<T>::imageDomainKernel(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, 
        double thres, double overDetermineRatio, 
        size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, 
        const std::vector<int>& oE1, const std::vector<int>& oE2, 
        ho7DArray<T>& ker)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(0)==acsDst.get_size(0));
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(1)==acsDst.get_size(1));
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(2)>=acsDst.get_size(2));
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(3)>=acsDst.get_size(3));

        size_t RO = acsSrc.get_size(0);
        size_t E1 = acsSrc.get_size(1);
        size_t E2 = acsSrc.get_size(2);
        size_t srcCHA = acsSrc.get_size(3);
        size_t dstCHA = acsDst.get_size(3);

        const T* pSrc = acsSrc.begin();
        const T* pDst = acsDst.begin();

        long long kROhalf = (long long)kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GWARN_STREAM("gtPlusGRAPPA<T>::calib3D(...) - 2*kROhalf == kRO " << kRO);
        }

        kRO = 2*kROhalf + 1;

        size_t kNE1 = kE1.size();
        size_t oNE1 = oE1.size();

        size_t kNE2 = kE2.size();
        size_t oNE2 = oE2.size();

        // allocate kernel
        GADGET_CHECK_RETURN_FALSE(ker.createArray(kRO, kNE1, kNE2, srcCHA, dstCHA, oNE1, oNE2));

        // loop over the calibration region and assemble the equation
        // Ax = b

        size_t sRO = kROhalf;
        size_t eRO = RO - kROhalf -1;

        size_t sE1 = std::abs(kE1[0]);
        size_t eE1 = E1 -1 - kE1[kNE1-1];

        size_t sE2 = std::abs(kE2[0]);
        size_t eE2 = E2 -1 - kE2[kNE2-1];

        size_t lenRO = eRO-kROhalf+1;
        size_t lenE1 = eE1-sE1+1;
        size_t lenE2 = eE2-sE2+1;

        size_t colA = kRO*kNE1*kNE2*srcCHA;
        size_t colB = dstCHA*oNE1*oNE2;

        if ( overDetermineRatio > 1.0 )
        {
            size_t maxRowA = (size_t)std::ceil(overDetermineRatio*colA);
            size_t maxROUsed = maxRowA/(lenE1*lenE2);
            if ( maxROUsed < lenRO )
            {
                // find the peak signal of acsSrc
                hoNDArray<T> acsSrc1stCha(RO, E1, E2, const_cast<T*>(acsSrc.begin()));
                hoNDArray<T> acsSrc1stChaSumE2(RO, E1, 1), acsSrc1stChaSumE2E1(RO, 1, 1);

                try
                {
                    Gadgetron::sum_over_dimension(acsSrc1stCha, acsSrc1stChaSumE2, 2);
                    Gadgetron::sum_over_dimension(acsSrc1stChaSumE2, acsSrc1stChaSumE2E1, 1);

                    T maxSignal;
                    size_t roInd(0);
                    try
                    {
                        Gadgetron::maxAbsolute(acsSrc1stChaSumE2E1, maxSignal, roInd);

                        if ( roInd > maxROUsed/2+kROhalf )
                        {
                            sRO = roInd - maxROUsed/2;
                        }
                        else
                        {
                            sRO = kROhalf;
                        }

                        if( sRO+maxROUsed-1 <= RO-kROhalf-1 )
                        {
                            eRO = sRO + maxROUsed - 1;
                        }
                        else
                        {
                            eRO = RO - kROhalf -1;
                        }

                        lenRO = eRO-sRO+1;
                        GDEBUG_STREAM("gtPlusGRAPPA<T>::calib3D(...) - overDetermineRatio = " << overDetermineRatio << " ; RO data range used : [" << sRO << " " << eRO << "] ...");
                    }
                    catch(...)
                    {
                        GWARN_STREAM("gtPlusGRAPPA<T>::calib3D(...) - overDetermineRatio is ignored ... ");
                        throw;
                    }
                }
                catch (...)
                {
                    GWARN_STREAM("gtPlusGRAPPA<T>::calib3D(...) - overDetermineRatio is ignored ... ");
                    throw;
                }
            }
        }

        size_t rowA = lenRO*lenE1*lenE2;

        hoMatrix<T> A, B, x( colA, colB );

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - allocate matrix storage ... "); }
        hoNDArrayMemoryManaged<T> A_mem(rowA, colA, gtPlus_mem_manager_);
        A.createMatrix( rowA, colA, A_mem.begin() );
        T* pA = A.begin();

        hoNDArrayMemoryManaged<T> B_mem(rowA, colB, gtPlus_mem_manager_);
        B.createMatrix( rowA, colB, B_mem.begin() );
        if ( performTiming_ ) { gt_timer3_.stop(); }
        T* pB = B.begin();

        long long e2;

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - fill calib matrixes ... "); }
        #pragma omp parallel for default(none) private(e2) shared(sE2, eE2, sE1, eE1, kROhalf, sRO, eRO, lenRO, lenE1, srcCHA, kNE2, kNE1, A, rowA, pA, acsSrc, kE1, kE2, oNE2, oNE1, dstCHA, B, pB, acsDst, oE1, oE2)
        for ( e2=(long long)sE2; e2<=(long long)eE2; e2++ )
        {
            long long e1;
            for ( e1=(long long)sE1; e1<=(long long)eE1; e1++ )
            {
                for ( long long ro=(long long)sRO; ro<=(long long)eRO; ro++ )
                {
                    size_t rInd = (e2-sE2)*lenRO*lenE1 + (e1-sE1)*lenRO + ro-sRO;

                    size_t src, dst, ke1, ke2, oe1, oe2;
                    long long kro;

                    // fill matrix A
                    size_t col = 0;
                    for ( src=0; src<srcCHA; src++ )
                    {
                        for ( ke2=0; ke2<kNE2; ke2++ )
                        {
                            for ( ke1=0; ke1<kNE1; ke1++ )
                            {
                                for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                                {
                                    // A(rInd, col++) = acsSrc(ro+kro, e1+kE1[ke1], e2+kE2[ke2], src);
                                    pA[rInd + col*rowA] = acsSrc(ro+kro, e1+kE1[ke1], e2+kE2[ke2], src);
                                    col++;
                                }
                            }
                        }
                    }

                    // fill matrix B
                    col = 0;
                    for ( oe2=0; oe2<oNE2; oe2++ )
                    {
                        for ( oe1=0; oe1<oNE1; oe1++ )
                        {
                            for ( dst=0; dst<dstCHA; dst++ )
                            {
                                // B(rInd, col++) = acsDst(ro, e1+oE1[oe1], e2+oE2[oe2], dst);
                                pB[rInd + col*rowA] = acsDst(ro, e1+oE1[oe1], e2+oE2[oe2], dst);
                                col++;
                            }
                        }
                    }
                }
            }
        }
        if ( performTiming_ ) { gt_timer3_.stop(); }

        //typename realType<T>::Type v;

        //Gadgetron::norm2(A, v);
        //GDEBUG_STREAM("A = " << v);

        //Gadgetron::norm2(B, v);
        //GDEBUG_STREAM("B = " << v);

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - solve linear system ... "); }
        //#ifdef USE_CUDA
        //    // go to device
        //    try
        //    {
        //        if ( typeid(typename realType<T>::Type)==typeid(float) && calib_use_gpu_ )
        //        {
        //            GDEBUG_STREAM("grappa 3D - calling GPU kernel estimation ... ");
        //            //hoNDArray<float_complext> A_tmp(A.get_dimensions(), reinterpret_cast<float_complext*>(A.begin()));
        //            //hoNDArray<float_complext> B_tmp(B.get_dimensions(), reinterpret_cast<float_complext*>(B.begin()));

        //            //cuNDArray<float_complext> device_A(A_tmp);
        //            //cuNDArray<float_complext> device_B(B_tmp);
        //            //cuNDArray<float_complext> device_x;
        //            //if ( Gadgetron::inverse_clib_matrix(&device_A, &device_B, &device_x, thres) != 0 )
        //            //{
        //            //    GERROR_STREAM("failed in Gadgetron::inverse_clib_matrix(&device_A, &device_B, &device_x, thres) ... ");
        //            //    SolveLinearSystem_Tikhonov(A, B, x, thres);
        //            //}
        //            //else
        //            //{
        //            //    // go back to host
        //            //    boost::shared_ptr< hoNDArray<complext<float> > > host_x = device_x.to_host();
        //            //    memcpy(x.begin(), host_x->begin(), x.get_number_of_bytes());
        //            //}

        //            hoNDArray<float_complext> A_tmp(A.get_dimensions(), reinterpret_cast<float_complext*>(A.begin()));
        //            hoNDArray<float_complext> B_tmp(B.get_dimensions(), reinterpret_cast<float_complext*>(B.begin()));

        //            int ret(0);
        //            boost::shared_ptr< hoNDArray<complext<float> > > host_x;

        //            #pragma omp critical(inverse3D)
        //            {
        //                cuNDArray<float_complext> device_A(A_tmp);
        //                cuNDArray<float_complext> device_B(B_tmp);
        //                cuNDArray<float_complext> device_x;

        //                ret = Gadgetron::inverse_clib_matrix(&device_A, &device_B, &device_x, thres);
        //                if ( ret == 0 )
        //                {
        //                    host_x = device_x.to_host();
        //                }
        //            }

        //            if ( ret != 0 )
        //            {
        //                GERROR_STREAM("failed in Gadgetron::inverse_clib_matrix(&device_A, &device_B, &device_x, thres) ... ");
        //                SolveLinearSystem_Tikhonov(A, B, x, thres);
        //            }
        //            else
        //            {
        //                memcpy(x.begin(), host_x->begin(), x.get_number_of_bytes());
        //            }
        //        }
        //        else
        //        {
        //            GWARN_STREAM("GPU inverse_clib_matrix for grappa is only available for single-precision, calling the CPU version ... ");
        //            GADGET_CHECK_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));
        //        }
        //    }
        //    catch(...)
        //    {
        //        GERROR_STREAM("failed in GPU inverse_clib_matrix for grappa, calling the CPU version ... ");
        //        GADGET_CHECK_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));
        //    }
        //#else
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));
        //#endif // USE_CUDA
        if ( performTiming_ ) { gt_timer3_.stop(); }

        //Gadgetron::norm2(x, v);
        //GDEBUG_STREAM("x = " << v);

        // the matrix dimension just matches
        //hoMatrix<T> xt(x.cols(), x.rows(), ker.begin());
        //GADGET_CHECK_RETURN_FALSE(Gadgetron::trans(x, xt));
        memcpy(ker.begin(), x.begin(), ker.get_number_of_bytes());

        //Gadgetron::norm2(ker, v);
        //GDEBUG_STREAM("ker = " << v);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusGRAPPA<T>::calib3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
kspaceDomainConvKernel3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, ho5DArray<T>& convKer, bool ROis3rdDim)
{
try
    {
        long long srcCHA = (long long)(ker.get_size(3));
        long long dstCHA = (long long)(ker.get_size(4));

        long long kNE1 = (long long)(kE1.size());
        long long oNE1 = (long long)(oE1.size());

        long long kNE2 = (long long)(kE2.size());
        long long oNE2 = (long long)(oE2.size());

        long long kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GWARN_STREAM("gtPlusGRAPPA<T>::imageDomainKernel(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - convert to conv kernel ... "); }
        /// fill the convolution kernels
        long long convKRO = 2*kRO+3;

        long long maxKE1 = std::abs(kE1[0]);
        if ( std::abs(kE1[kNE1-1]) > maxKE1 )
        {
            maxKE1 = std::abs(kE1[kNE1-1]);
        }
        long long convKE1 = 2*maxKE1+1;

        long long maxKE2 = std::abs(kE2[0]);
        if ( std::abs(kE2[kNE2-1]) > maxKE2 )
        {
            maxKE2 = std::abs(kE2[kNE2-1]);
        }
        long long convKE2 = 2*maxKE2+1;

        /// allocate the convolution kernel
        if ( ROis3rdDim )
        {
            convKer.createArray(convKE1, convKE2, convKRO, srcCHA, dstCHA);
        }
        else
        {
            convKer.createArray(convKRO, convKE1, convKE2, srcCHA, dstCHA);
        }
        Gadgetron::clear(&convKer);

        /// index
        long long oe1, oe2, kro, ke1, ke2, src, dst;

        /// fill the convolution kernel and sum up multiple kernels
        for ( oe2=0; oe2<oNE2; oe2++ )
        {
            for ( oe1=0; oe1<oNE1; oe1++ )
            {
                for ( ke2=0; ke2<kNE2; ke2++ )
                {
                    for ( ke1=0; ke1<kNE1; ke1++ )
                    {
                        for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                        {
                            for ( dst=0; dst<dstCHA; dst++ )
                            {
                                if ( ROis3rdDim )
                                {
                                    for ( src=0; src<srcCHA; src++ )
                                    {
                                        convKer(oE1[oe1]-kE1[ke1]+maxKE1, oE2[oe2]-kE2[ke2]+maxKE2, -kro+kRO+1, src, dst) = ker(kro+kROhalf, ke1, ke2, src, dst, oe1, oe2);
                                    }
                                }
                                else
                                {
                                    for ( src=0; src<srcCHA; src++ )
                                    {
                                        convKer(-kro+kRO+1, oE1[oe1]-kE1[ke1]+maxKE1, oE2[oe2]-kE2[ke2]+maxKE2, src, dst) = ker(kro+kROhalf, ke1, ke2, src, dst, oe1, oe2);
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }

        if ( (oE1[0]!=0) && (oE2[0]!=0) && (srcCHA==dstCHA) )
        {
            for ( dst=0; dst<dstCHA; dst++ )
            {
                if ( ROis3rdDim )
                {
                    for ( src=0; src<srcCHA; src++ )
                    {
                        if ( src == dst )
                        {
                            convKer(maxKE1, maxKE2, kRO+1, src, dst) = 1.0;
                        }
                        else
                        {
                            convKer(maxKE1, maxKE2, kRO+1, src, dst) = 0.0;
                        }
                    }
                }
                else
                {
                    for ( src=0; src<srcCHA; src++ )
                    {
                        if ( src == dst )
                        {
                            convKer(kRO+1, maxKE1, maxKE2, src, dst) = 1.0;
                        }
                        else
                        {
                            convKer(kRO+1, maxKE1, maxKE2, src, dst) = 0.0;
                        }
                    }
                }
            }
        }
        if ( performTiming_ ) { gt_timer3_.stop(); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusGRAPPA<T>::kspaceDomainConvKernel3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
imageDomainKernel3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, size_t ro, size_t e1, size_t e2, hoNDArray<T>& kIm)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(3));
        long long dstCHA = (long long)(ker.get_size(4));

        long long kNE1 = (long long)(kE1.size());
        long long oNE1 = (long long)(oE1.size());

        long long kNE2 = (long long)(kE2.size());
        long long oNE2 = (long long)(oE2.size());

        long long kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GWARN_STREAM("gtPlusGRAPPA<T>::imageDomainKernel(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        // allocate image domain kernel
        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - create kIm array ... "); }
        if ( kIm.get_number_of_elements() < (size_t)ro*e1*e2*srcCHA*dstCHA )
        {
            kIm.create(ro, e1, e2, srcCHA, dstCHA);
        }
        if ( performTiming_ ) { gt_timer3_.stop(); }

        ho5DArray<T> convKer;
        bool ROis3rdDim = false;
        GADGET_CHECK_RETURN_FALSE(this->kspaceDomainConvKernel3D(ker, kRO, kE1, kE2, oE1, oE2, convKer, ROis3rdDim));

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - SNR unit scaling ... "); }
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1*e2)) ), convKer ));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - zero padding ... "); }
        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3DNoPresetZeros(convKer, ro, e1, e2, kIm));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(convKer, ro, e1, e2, kIm, false));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - conver to image domain ... "); }
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kIm);
        if ( performTiming_ ) { gt_timer3_.stop(); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusGRAPPA<T>::imageDomainKernel3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
imageDomainKernelRO3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, size_t ro, hoNDArray<T>& kImRO)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(3));
        long long dstCHA = (long long)(ker.get_size(4));

        GADGET_CHECK_RETURN_FALSE(kRO==ker.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kE1.size()==ker.get_size(1));
        GADGET_CHECK_RETURN_FALSE(kE2.size()==ker.get_size(2));
        GADGET_CHECK_RETURN_FALSE(oE1.size()==ker.get_size(5));
        GADGET_CHECK_RETURN_FALSE(oE2.size()==ker.get_size(6));

        bool ROat3rdDim = false;
        ho5DArray<T> convKer;
        GADGET_CHECK_RETURN_FALSE(this->kspaceDomainConvKernel3D(ker, kRO, kE1,  kE2, oE1, oE2, convKer, ROat3rdDim));

        // allocate image domain kernel
        size_t kConvE1 = convKer.get_size(1);
        size_t kConvE2 = convKer.get_size(2);

        kImRO.create(kConvE1, kConvE2, ro, srcCHA, dstCHA);

        hoNDArray<T> kImROTemp(ro, kConvE1, kConvE2, srcCHA, dstCHA);
        Gadgetron::clear(kImROTemp);

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - SNR unit scaling ... "); }
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro)) ), convKer ));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(convKer, debugFolder_+"convKer_scal_RO"); }

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - zero padding only for RO ... "); }
        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3DNoPresetZeros(convKer, ro, kConvE1, kConvE2, kImROTemp));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(convKer, ro, kConvE1, kConvE2, kImROTemp, false));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kImROTemp, debugFolder_+"convKer_scal_RO_zeropadded"); }

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - conver to image domain only for RO ... "); }
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(kImROTemp);
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - permute kernel dimensions to be [kE1 kE2 RO ...]  ... "); }

        std::vector<size_t> dim_order(3);
        dim_order[0] = 1;
        dim_order[1] = 2;
        dim_order[2] = 0;

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&kImROTemp, &kImRO, &dim_order));

        if ( performTiming_ ) { gt_timer3_.stop(); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusGRAPPA<T>::imageDomainKernelRO3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
imageDomainKernelE1E2RO(const hoNDArray<T>& kImRO, size_t e1, size_t e2, hoNDArray<T>& kImE1E2RO)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dim = kImRO.get_dimensions();

        std::vector<size_t> dimR(*dim);
        dimR[0] = e1;
        dimR[1] = e2;

        kImE1E2RO.create(&dimR);
        Gadgetron::clear(kImE1E2RO);

        hoNDArray<T> kImROScaled(kImRO);

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - SNR unit scaling for E1 and E2 ... "); }
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(e1*e2)) ), kImROScaled ));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kImROScaled, debugFolder_+"kImROScaledE1E2"); }

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - zero padding for E1 and E2 ... "); }
        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3DNoPresetZeros(kImROScaled, e1, e2, dimR[2], kImE1E2RO));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(kImROScaled, e1, e2, dimR[2], kImE1E2RO, false));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kImE1E2RO, debugFolder_+"kImE1E2RO_zeropadded_E1E2"); }

        if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - conver to image domain for E1 and E2 ... "); }
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kImE1E2RO);
        if ( performTiming_ ) { gt_timer3_.stop(); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusGRAPPA<T>::imageDomainKernelE1E2RO(...) ... ");
        return false;
    }

    return true;
}

}}
