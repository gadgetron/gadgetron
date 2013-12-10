
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
    Magnetic Resonance in Medicine 2009;61(6): 1425�1433. 
*/

#pragma once

#include "gtPlusAlgorithmBase.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusGRAPPA : public gtPlusAlgorithmBase<T>
{
public:

    typedef gtPlusAlgorithmBase<T> BaseClass;

    gtPlusGRAPPA() : BaseClass() {}
    virtual ~gtPlusGRAPPA() {}

    virtual void printInfo(std::ostream& os);

    // get the kernel pattern, given the acceleration factor and kernel size
    bool kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, int accelFactor, size_t kNE1, bool fitItself);

    // grappa calibration for 2D case
    // acsSrc : [RO E1 srcCHA]
    // acsDst : [RO E1 dstCHA]
    // ker : [kRO kE1 srcCHA dstCHA oE1]
    bool calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, 
            int kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray<T>& ker);

    // image domain kernel for 2D kernel
    // kIm: image domain kernel [RO E1 srcCHA dstCHA]
    bool imageDomainKernel(const ho5DArray<T>& ker, int kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, int ro, int e1, hoNDArray<T>& kIm);

    // grappa calibration for 3D case
    // acsSrc : [RO E1 E2 srcCHA]
    // acsDst : [RO E1 E2 dstCHA]
    // ker : [kRO kE1 kE2 srcCHA dstCHA oE1 oE2]
    bool calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, double thres, 
            int kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, ho7DArray<T>& ker);

    // image domain kernel for 3D kernel
    // kIm: image domain kernel [RO E1 E2 srcCHA dstCHA]
    bool imageDomainKernel3D(const ho7DArray<T>& ker, int kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, int ro, int e1, int e2, hoNDArray<T>& kIm);

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
kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, int accelFactor, size_t kNE1, bool fitItself)
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
        int k;
        for ( k=-((int)kNE1/2-1); k<=(int)kNE1/2; k++ )
        {
            kE1[k+kNE1/2-1] = k*accelFactor;
        }
    }
    else
    {
        int k;
        for ( k=-(int)kNE1/2; k<=(int)kNE1/2; k++ )
        {
            kE1[k+kNE1/2] = k*accelFactor;
        }
    }

    if ( fitItself )
    {
        oE1.resize(accelFactor);
        for ( int a=0; a<accelFactor; a++ )
        {
            oE1[a] = a;
        }
    }
    else
    {
        oE1.resize(accelFactor-1);
        for ( int a=1; a<accelFactor; a++ )
        {
            oE1[a-1] = a;
        }
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, 
    int kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray<T>& ker)
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

        int kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GADGET_WARN_MSG("gtPlusGRAPPA<T>::calib(...) - 2*kROhalf == kRO " << kRO);
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
        hoMatrix<T> x( A.cols(), B.cols() );

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("grappa 2D calibration - allocate matrix storage ... "));
        hoNDArrayMemoryManaged<T> A_mem(colA, rowA, gtPlus_mem_manager_);
        A.createMatrix( rowA, colA, A_mem.begin() );

        hoNDArrayMemoryManaged<T> B_mem(colB, rowA, gtPlus_mem_manager_);
        B.createMatrix( A.rows(), colB, B_mem.begin() );
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        int e1;
        for ( e1=(int)sE1; e1<=(int)eE1; e1++ )
        {
            for ( int ro=kROhalf; ro<=(int)eRO; ro++ )
            {
                int rInd = (e1-sE1)*lenRO+ro-kROhalf;

                size_t src, dst, ke1, oe1;
                int kro;

                // fill matrix A
                size_t col = 0;
                for ( src=0; src<srcCHA; src++ )
                {
                    for ( ke1=0; ke1<kNE1; ke1++ )
                    {
                        for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                        {
                            A(rInd, col++) = acsSrc(ro+kro, e1+kE1[ke1], src);
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
        //GADGET_MSG("A = " << v);

        //Gadgetron::norm2(B, v);
        //GADGET_MSG("B = " << v);

        //GADGET_CHECK_PERFORM(performTiming_, gt_timer2_.start("SolveLinearSystem_Tikhonov"));
        GADGET_CHECK_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));
        //GADGET_CHECK_PERFORM(performTiming_, gt_timer2_.stop());

        //Gadgetron::norm2(x, v);
        //GADGET_MSG("x = " << v);

        // the matrix dimension just matches
        hoMatrix<T> xt(x.cols(), x.rows(), ker.begin());
        GADGET_CHECK_RETURN_FALSE(Gadgetron::trans(x, xt));

        //Gadgetron::norm2(ker, v);
        //GADGET_MSG("ker = " << v);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusGRAPPA<T>::calib(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
imageDomainKernel(const ho5DArray<T>& ker, int kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, int ro, int e1, hoNDArray<T>& kIm)
{
    try
    {
        int srcCHA = (int)(ker.get_size(2));
        int dstCHA = (int)(ker.get_size(3));
        int kNE1 = (int)(kE1.size());
        int oNE1 = (int)(oE1.size());

        int kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GADGET_WARN_MSG("gtPlusGRAPPA<T>::imageDomainKernel(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        // allocate image domain kernel
        kIm.create(ro, e1, srcCHA, dstCHA);

        /// fill the convolution kernels
        int convKRO = 2*kRO+3;

        int maxKE1 = std::abs(kE1[0]);
        if ( std::abs(kE1[kNE1-1]) > maxKE1 )
        {
            maxKE1 = std::abs(kE1[kNE1-1]);
        }
        int convKE1 = 2*maxKE1+1;

        /// allocate the convolution kernel
        ho4DArray<T> convKer(convKRO, convKE1, srcCHA, dstCHA);
        Gadgetron::clear(&convKer);

        /// index
        int oe1, kro, ke1, src, dst;

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

        GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1)) ), convKer ));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad2D(convKer, ro, e1, kIm));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kIm));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusGRAPPA<T>::imageDomainKernel(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, double thres, 
            int kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, 
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

        int kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GADGET_WARN_MSG("gtPlusGRAPPA<T>::calib3D(...) - 2*kROhalf == kRO " << kRO);
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

        size_t eRO = RO - kROhalf -1;

        size_t sE1 = std::abs(kE1[0]);
        size_t eE1 = E1 -1 - kE1[kNE1-1];

        size_t sE2 = std::abs(kE2[0]);
        size_t eE2 = E2 -1 - kE2[kNE2-1];

        size_t lenRO = eRO-kROhalf+1;
        size_t lenE1 = eE1-sE1+1;

        size_t rowA = (eE2-sE2+1)*lenE1*lenRO;
        size_t colA = kRO*kNE1*kNE2*srcCHA;
        size_t colB = dstCHA*oNE1*oNE2;

        hoMatrix<T> A, B, x( A.cols(), B.cols() );

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("grappa 3D calibration - allocate matrix storage ... "));
        hoNDArrayMemoryManaged<T> A_mem(colA, rowA, gtPlus_mem_manager_);
        A.createMatrix( rowA, colA, A_mem.begin() );

        hoNDArrayMemoryManaged<T> B_mem(colB, rowA, gtPlus_mem_manager_);
        B.createMatrix( A.rows(), colB, B_mem.begin() );
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        int e2;

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("grappa 3D calibration - fill calib matrixes ... "));
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel for default(none) private(e2) shared(sE2, eE2, sE1, eE1, kROhalf, eRO, lenRO, lenE1, srcCHA, kNE2, kNE1, A, oNE2, oNE1, dstCHA, B)
        #else
            #pragma omp parallel for default(none) private(e2) shared(sE2, eE2, sE1, eE1, kROhalf, eRO, lenRO, lenE1, srcCHA, kNE2, kNE1, A, acsSrc, kE1, kE2, oNE2, oNE1, dstCHA, B, acsDst, oE1, oE2)
        #endif
        for ( e2=(int)sE2; e2<=(int)eE2; e2++ )
        {
            int e1;
            for ( e1=(int)sE1; e1<=(int)eE1; e1++ )
            {
                for ( int ro=kROhalf; ro<=(int)eRO; ro++ )
                {
                    int rInd = (e2-sE2)*lenRO*lenE1 + (e1-sE1)*lenRO + ro-kROhalf;

                    size_t src, dst, ke1, ke2, oe1, oe2;
                    int kro;

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
                                    A(rInd, col++) = acsSrc(ro+kro, e1+kE1[ke1], e2+kE2[ke2], src);
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
                                B(rInd, col++) = acsDst(ro, e1+oE1[oe1], e2+oE2[oe2], dst);
                            }
                        }
                    }
                }
            }
        }
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        //typename realType<T>::Type v;

        //Gadgetron::norm2(A, v);
        //GADGET_MSG("A = " << v);

        //Gadgetron::norm2(B, v);
        //GADGET_MSG("B = " << v);

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("grappa 3D calibration - solve linear system ... "));
        GADGET_CHECK_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        //Gadgetron::norm2(x, v);
        //GADGET_MSG("x = " << v);

        // the matrix dimension just matches
        hoMatrix<T> xt(x.cols(), x.rows(), ker.begin());
        GADGET_CHECK_RETURN_FALSE(Gadgetron::trans(x, xt));

        //Gadgetron::norm2(ker, v);
        //GADGET_MSG("ker = " << v);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusGRAPPA<T>::calib3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusGRAPPA<T>::
imageDomainKernel3D(const ho7DArray<T>& ker, int kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, int ro, int e1, int e2, hoNDArray<T>& kIm)
{
    try
    {
        int srcCHA = (int)(ker.get_size(3));
        int dstCHA = (int)(ker.get_size(4));

        int kNE1 = (int)(kE1.size());
        int oNE1 = (int)(oE1.size());

        int kNE2 = (int)(kE2.size());
        int oNE2 = (int)(oE2.size());

        int kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GADGET_WARN_MSG("gtPlusGRAPPA<T>::imageDomainKernel(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        // allocate image domain kernel
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("grappa 3D calibration - create kIm array ... "));
        if ( kIm.get_number_of_elements() < (size_t)ro*e1*e2*srcCHA*dstCHA )
        {
            kIm.create(ro, e1, e2, srcCHA, dstCHA);
        }
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("grappa 3D calibration - convert to conv kernel ... "));
        /// fill the convolution kernels
        int convKRO = 2*kRO+3;

        int maxKE1 = std::abs(kE1[0]);
        if ( std::abs(kE1[kNE1-1]) > maxKE1 )
        {
            maxKE1 = std::abs(kE1[kNE1-1]);
        }
        int convKE1 = 2*maxKE1+1;

        int maxKE2 = std::abs(kE2[0]);
        if ( std::abs(kE2[kNE2-1]) > maxKE2 )
        {
            maxKE2 = std::abs(kE2[kNE2-1]);
        }
        int convKE2 = 2*maxKE2+1;

        /// allocate the convolution kernel
        ho5DArray<T> convKer(convKRO, convKE1, convKE2, srcCHA, dstCHA);
        Gadgetron::clear(&convKer);

        /// index
        int oe1, oe2, kro, ke1, ke2, src, dst;

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

        if ( (oE1[0]!=0) && (oE2[0]!=0) && (srcCHA==dstCHA) )
        {
            for ( dst=0; dst<dstCHA; dst++ )
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
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("grappa 3D calibration - SNR unit scaling ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1*e2)) ), convKer ));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("grappa 3D calibration - zero padding ... "));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3D(convKer, ro, e1, e2, kIm));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("grappa 3D calibration - conver to image domain ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kIm));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusGRAPPA<T>::imageDomainKernel3D(...) ... ");
        return false;
    }

    return true;
}

}}
