
/** \file   mri_core_grappa.cpp
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

#include "mri_core_grappa.h"
#include "ismrmrd/ismrmrd.h"
#include "GadgetronTimer.h"
#include "mri_core_utility.h"
#include "mri_core_utility.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

template <typename T> 
void grappa<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "---------------------------- GRAPPA reconstruction ------------------" << endl;
    os << "Implementation of GRAPPA algorithms for ISMRMRD package" << endl;
    os << "Both 2D and 3D version are implemented" << endl;
    os << "Algorithms are published at:" << endl;
    os << "Generalized autocalibrating partially parallel acquisitions (GRAPPA), Magnetic Resonance in Medicine, Volume 47, Issue 6, pages 1202�1210, June 2002" << endl;
    os << "HTGRAPPA: Real-time B1-weighted image domain TGRAPPA reconstruction, Magnetic Resonance in Medicine, Volume 61, Issue 6, pages 1425�1433, June 2009" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
void grappa<T>::
kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, size_t accelFactor, size_t kNE1, bool fitItself)
{
    if ( accelFactor == 1 )
    {
        kE1.resize(1, 0);
        oE1.resize(1, 0);
        return;
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

    return;
}

template <typename T> 
void grappa<T>::
calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, 
    size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray<T>& ker)
{
    try
    {
        GADGET_CHECK_THROW(acsSrc.get_size(0)==acsDst.get_size(0));
        GADGET_CHECK_THROW(acsSrc.get_size(1)==acsDst.get_size(1));
        GADGET_CHECK_THROW(acsSrc.get_size(2)>=acsDst.get_size(2));

        size_t RO = acsSrc.get_size(0);
        size_t E1 = acsSrc.get_size(1);
        size_t srcCHA = acsSrc.get_size(2);
        size_t dstCHA = acsDst.get_size(2);

        const T* pSrc = acsSrc.begin();
        const T* pDst = acsDst.begin();

        long long kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GADGET_WARN_MSG("grappa<T>::calib(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        size_t kNE1 = kE1.size();
        size_t oNE1 = oE1.size();

        /// allocate kernel
        GADGET_CHECK_THROW(ker.createArray(kRO, kNE1, srcCHA, dstCHA, oNE1));

        /// loop over the calibration region and assemble the equation
        /// Ax = b

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

        hoNDArray<T> A_mem(rowA, colA);
        A.createMatrix( rowA, colA, A_mem.begin() );
        T* pA = A.begin();

        hoNDArray<T> B_mem(rowA, colB);
        B.createMatrix( A.rows(), colB, B_mem.begin() );
        T* pB = B.begin();

        long long e1;
        for ( e1=(long long)sE1; e1<=(long long)eE1; e1++ )
        {
            for ( long long ro=kROhalf; ro<=(long long)eRO; ro++ )
            {
                long long rInd = (e1-sE1)*lenRO+ro-kROhalf;

                size_t src, dst, ke1, oe1;
                long long kro;

                /// fill matrix A
                size_t col = 0;
                size_t offset = 0;
                for ( src=0; src<srcCHA; src++ )
                {
                    for ( ke1=0; ke1<kNE1; ke1++ )
                    {
                        offset = src*RO*E1 + (e1+kE1[ke1])*RO;
                        for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                        {
                            /// A(rInd, col++) = acsSrc(ro+kro, e1+kE1[ke1], src);
                            pA[rInd + col*rowA] = pSrc[ro+kro+offset];
                            col++;
                        }
                    }
                }

                /// fill matrix B
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

        SolveLinearSystem_Tikhonov(A, B, x, thres);
        memcpy(ker.begin(), x.begin(), ker.get_number_of_bytes());
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa<T>::calib(...) ... ");
    }

    return;
}

template <typename T> 
void grappa<T>::
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
            GADGET_WARN_MSG("grappa<T>::imageDomainKernel(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        /// allocate image domain kernel
        kIm.create(ro, e1, srcCHA, dstCHA);

        //// fill the convolution kernels
        long long convKRO = 2*kRO+3;

        long long maxKE1 = std::abs(kE1[0]);
        if ( std::abs(kE1[kNE1-1]) > maxKE1 )
        {
            maxKE1 = std::abs(kE1[kNE1-1]);
        }
        long long convKE1 = 2*maxKE1+1;

        //// allocate the convolution kernel
        ho4DArray<T> convKer(convKRO, convKE1, srcCHA, dstCHA);
        Gadgetron::clear(&convKer);

        //// index
        long long oe1, kro, ke1, src, dst;

        //// fill the convolution kernel and sum up multiple kernels
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

        if ( oE1[0] != 0 )
        {
            for ( dst=0; dst<dstCHA; dst++ )
            {
                convKer(kRO+1, maxKE1, dst, dst) = 1.0;
            }
        }

        Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1)) ), convKer );
        Gadgetron::zeropad2D(convKer, ro, e1, kIm);
        GADGET_CHECK_THROW(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kIm));
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa<T>::imageDomainKernel(...) ... ");
    }

    return;
}

template <typename T> 
void grappa<T>::
calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, 
        double thres, double overDetermineRatio, 
        size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, 
        const std::vector<int>& oE1, const std::vector<int>& oE2, 
        ho7DArray<T>& ker)
{
    try
    {
        GADGET_CHECK_THROW(acsSrc.get_size(0)==acsDst.get_size(0));
        GADGET_CHECK_THROW(acsSrc.get_size(1)==acsDst.get_size(1));
        GADGET_CHECK_THROW(acsSrc.get_size(2)>=acsDst.get_size(2));
        GADGET_CHECK_THROW(acsSrc.get_size(3)>=acsDst.get_size(3));

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
            GADGET_WARN_MSG("grappa<T>::calib3D(...) - 2*kROhalf == kRO " << kRO);
        }

        kRO = 2*kROhalf + 1;

        size_t kNE1 = kE1.size();
        size_t oNE1 = oE1.size();

        size_t kNE2 = kE2.size();
        size_t oNE2 = oE2.size();

        /// allocate kernel
        GADGET_CHECK_THROW(ker.createArray(kRO, kNE1, kNE2, srcCHA, dstCHA, oNE1, oNE2));

        /// loop over the calibration region and assemble the equation
        /// Ax = b

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
                /// find the peak signal of acsSrc
                hoNDArray<T> acsSrc1stCha(RO, E1, E2, const_cast<T*>(acsSrc.begin()));
                hoNDArray<T> acsSrc1stChaSumE2(RO, E1, 1), acsSrc1stChaSumE2E1(RO, 1, 1);

                if ( Gadgetron::sumOver3rdDimension(acsSrc1stCha, acsSrc1stChaSumE2) )
                {
                    if ( Gadgetron::sumOver2ndDimension(acsSrc1stChaSumE2, acsSrc1stChaSumE2E1) )
                    {
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
                            GADGET_MSG("grappa<T>::calib3D(...) - overDetermineRatio = " << overDetermineRatio << " ; RO data range used : [" << sRO << " " << eRO << "] ...");
                        }
                        catch(...)
                        {
                            GADGET_WARN_MSG("grappa<T>::calib3D(...) - overDetermineRatio is ignored ... ");
                        }
                    }
                }
                else
                {
                    GADGET_WARN_MSG("grappa<T>::calib3D(...) - overDetermineRatio is ignored ... ");
                }
            }
        }

        size_t rowA = lenRO*lenE1*lenE2;

        hoMatrix<T> A, B, x( colA, colB );

        hoNDArray<T> A_mem(rowA, colA);
        A.createMatrix( rowA, colA, A_mem.begin() );
        T* pA = A.begin();

        hoNDArray<T> B_mem(rowA, colB);
        B.createMatrix( rowA, colB, B_mem.begin() );
        T* pB = B.begin();

        long long e2;

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

                    /// fill matrix A
                    size_t col = 0;
                    for ( src=0; src<srcCHA; src++ )
                    {
                        for ( ke2=0; ke2<kNE2; ke2++ )
                        {
                            for ( ke1=0; ke1<kNE1; ke1++ )
                            {
                                for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                                {
                                    /// A(rInd, col++) = acsSrc(ro+kro, e1+kE1[ke1], e2+kE2[ke2], src);
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

        SolveLinearSystem_Tikhonov(A, B, x, thres);
        memcpy(ker.begin(), x.begin(), ker.get_number_of_bytes());
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa<T>::calib3D(...) ... ");
    }

    return;
}

template <typename T> 
void grappa<T>::
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
            GADGET_WARN_MSG("grappa<T>::imageDomainKernel(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

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

        if ( (oE1[0]!=0) && (oE2[0]!=0) )
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
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa<T>::kspaceDomainConvKernel3D(...) ... ");
    }

    return;
}

template <typename T> 
void grappa<T>::
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
            GADGET_WARN_MSG("grappa<T>::imageDomainKernel(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        // allocate image domain kernel
        if ( kIm.get_number_of_elements() < (size_t)ro*e1*e2*srcCHA*dstCHA )
        {
            kIm.create(ro, e1, e2, srcCHA, dstCHA);
        }

        ho5DArray<T> convKer;
        bool ROis3rdDim = false;
        this->kspaceDomainConvKernel3D(ker, kRO, kE1, kE2, oE1, oE2, convKer, ROis3rdDim);

        Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1*e2)) ), convKer );

        Gadgetron::zeropad3DNoPresetZeros(convKer, ro, e1, e2, kIm);

        GADGET_CHECK_THROW(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kIm));
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa<T>::imageDomainKernel3D(...) ... ");
    }

    return;
}

template <typename T> 
void grappa<T>::
imageDomainKernelRO3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, size_t ro, hoNDArray<T>& kImRO)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(3));
        long long dstCHA = (long long)(ker.get_size(4));

        GADGET_CHECK_THROW(kRO==ker.get_size(0));
        GADGET_CHECK_THROW(kE1.size()==ker.get_size(1));
        GADGET_CHECK_THROW(kE2.size()==ker.get_size(2));
        GADGET_CHECK_THROW(oE1.size()==ker.get_size(5));
        GADGET_CHECK_THROW(oE2.size()==ker.get_size(6));

        bool ROat3rdDim = false;
        ho5DArray<T> convKer;
        this->kspaceDomainConvKernel3D(ker, kRO, kE1,  kE2, oE1, oE2, convKer, ROat3rdDim);

        // allocate image domain kernel
        size_t kConvE1 = convKer.get_size(1);
        size_t kConvE2 = convKer.get_size(2);

        kImRO.create(kConvE1, kConvE2, ro, srcCHA, dstCHA);

        hoNDArray<T> kImROTemp(ro, kConvE1, kConvE2, srcCHA, dstCHA);
        Gadgetron::clear(kImROTemp);

        Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro)) ), convKer );

        Gadgetron::zeropad3DNoPresetZeros(convKer, ro, kConvE1, kConvE2, kImROTemp);

        GADGET_CHECK_THROW(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(kImROTemp));

        GADGET_CHECK_THROW(Gadgetron::permuteROTo3rdDimensionFor3DRecon(kImROTemp, kImRO));
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa<T>::imageDomainKernelRO3D(...) ... ");
    }

    return;
}

template <typename T> 
void grappa<T>::
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

        Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(e1*e2)) ), kImROScaled );

        Gadgetron::zeropad3DNoPresetZeros(kImROScaled, e1, e2, dimR[2], kImE1E2RO);

        GADGET_CHECK_THROW(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kImE1E2RO));
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa<T>::imageDomainKernelE1E2RO(...) ... ");
    }

    return;
}

// 
// Template instantiation
//

template class EXPORTMRICORE grappa< std::complex<float> >;
template class EXPORTMRICORE grappa< std::complex<double> >;

}
