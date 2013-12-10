
/** \file   gtPlusSPIRIT.h
    \brief  SPIRIT kernel estimation for 2D and 3D MRI parallel imaging
    \author Hui Xue

    References to the implementation can be found in:

    Lustig M, Pauly JM. 
    SPIRiT: Iterative self-consistent parallel imaging reconstruction from arbitrary k-space. 
    Magnetic Resonance in Medicine 2010;64(2):457-471.

    ISMRM 2013 sunrise course on Parallel Imaging
    Michael S. Hansen, Philip Beatty
    http://gadgetron.sourceforge.net/sunrise/
    http://cds.ismrm.org/protected/13MPresentations/abstracts/7059.pdf
*/

#pragma once

#include "gtPlusAlgorithmBase.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusSPIRIT : public gtPlusAlgorithmBase<T>
{
public:

    typedef gtPlusAlgorithmBase<T> BaseClass;

    typedef typename realType<T>::Type ValueType;

    gtPlusSPIRIT() : BaseClass() {}
    virtual ~gtPlusSPIRIT() {}

    virtual void printInfo(std::ostream& os);

    // SPIRIT calibration for 2D case
    // acsSrc : [RO E1 srcCHA]
    // acsDst : [RO E1 dstCHA]
    // ker : [kRO kE1 srcCHA dstCHA 1 1]
    bool calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, 
            int kRO, int kE1, int a, int b, ho6DArray<T>& ker);

    // image domain kernel for 2D kernel
    // kIm: image domain kernel [RO E1 srcCHA dstCHA]
    // if minusI==true, compute image domain G-I kernel
    bool imageDomainKernel(const ho6DArray<T>& ker, int kRO, int kE1, int a, int b, int ro, int e1, hoNDArray<T>& kIm, bool minusI=false);

    // SPIRIT calibration for 3D case
    // acsSrc : [RO E1 E2 srcCHA]
    // acsDst : [RO E1 E2 dstCHA]
    // ker : [kRO kE1 kE2 srcCHA dstCHA 1 1 1]
    // overDetermineRatio : over determine ratio of calib matrix, if < 1, all data are used
    bool calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, double thres, double overDetermineRatio, 
            int kRO, int kE1, int kE2, int a, int b, int c, hoNDArray<T>& ker);

    // convert the calibrated kernel to the convlution kernel in kspace
    // if ROis3rdDim == true, the kernel dimension is [E1 E2 RO], otherwise [RO E1 E2]
    bool kspaceDomainConvKernel3D(const hoNDArray<T>& ker, int kRO, int kE1, int kE2, int a, int b, int c, ho5DArray<T>& convKerFlip, bool minusI=true, bool ROis3rdDim=true);

    // image domain kernel for 3D kernel
    // kIm: image domain kernel [E1 E2 RO srcCHA dstCHA]
    // if minusI==true, compute image domain G-I kernel
    bool imageDomainKernel3D(const hoNDArray<T>& ker, int kRO, int kE1, int kE2, 
        int a, int b, int c, int ro, int e1, int e2, hoNDArray<T>& kIm, bool minusI=false);

    // image domain kernel for 3D kernel, only RO direction is converted to image domain
    // E1 and E2 stays in the kspace domain
    // kImRO: kspace-image hybrid kernel [convE1 convE2 RO srcCHA dstCHA]
    bool imageDomainKernelRO3D(const hoNDArray<T>& ker, int kRO, int kE1, int kE2, 
        int a, int b, int c, int ro, hoNDArray<T>& kImRO, bool minusI=false);

    // image domain kernel for 3D kernel, E1 and E2 directions are converted to image domain
    // kImRO : kspace-image hybrid kernel where first two dimensions are E1 and E2 and in kspace
    bool imageDomainKernelE1E2RO(const hoNDArray<T>& kImRO, int e1, int e2, hoNDArray<T>& kImE1E2RO);

    // compute the image domain adjoint kernel
    bool imageDomainAdjointKernel(const hoNDArray<T>& kIm, hoNDArray<T>& adjkIm);

    // compute the (G-I)'*(G-I)
    bool AdjointForwardKernel(const hoNDArray<T>& kImS2D, const hoNDArray<T>& kImD2S, hoNDArray<T>& kIm);

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
void gtPlusSPIRIT<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT reconstruction ------------------" << endl;
    os << "Implementation of SPIRIT algorithms for ISMRMRD package" << endl;
    os << "Both 2D and 3D version are implemented" << endl;
    os << "Algorithms are published at:" << endl;
    os << "Lustig, M. and Pauly, J. M. (2010), SPIRiT: Iterative self-consistent parallel imaging reconstruction from arbitrary k-space. Magn Reson Med, 64: 457-471. doi: 10.1002/mrm.22428" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, int kRO, int kE1, int a, int b, ho6DArray<T>& ker)
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

        int kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GADGET_WARN_MSG("gtPlusSPIRIT<T>::calib(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        int kE1half = kE1/2;
        if ( 2*kE1half == kE1 )
        {
            GADGET_WARN_MSG("gtPlusSPIRIT<T>::calib(...) - 2*kE1half == kE1 " << kE1);
        }
        kE1 = 2*kE1half + 1;

        // allocate kernel
        GADGET_CHECK_RETURN_FALSE(ker.createArray(kRO, kE1, srcCHA, dstCHA, 1, 1));

        // loop over the calibration region and assemble the equation
        // Ax = b

        size_t sRO = kROhalf;
        size_t eRO = RO - kROhalf -1;

        size_t sE1 = kE1half;
        size_t eE1 = E1 - kE1half -1;

        size_t lenRO = eRO-sRO+1;
        size_t lenE1 = eE1-sE1+1;

        size_t rowA = lenE1*lenRO;
        size_t colA = (kRO*kE1-1)*srcCHA;
        size_t colB = dstCHA;

        hoMatrix<T> A;
        hoMatrix<T> B;
        hoMatrix<T> x( A.cols(), B.cols() );

        hoNDArrayMemoryManaged<T> A_mem(colA, rowA, gtPlus_mem_manager_);
        A.createMatrix( rowA, colA, A_mem.begin() );

        hoNDArrayMemoryManaged<T> B_mem(colB, rowA, gtPlus_mem_manager_);
        B.createMatrix( A.rows(), colB, B_mem.begin() );

        int dRO, dE1;

        for ( int e1=(int)sE1; e1<=(int)eE1; e1++ )
        {
            dE1 = e1;

            for ( int ro=sRO; ro<=(int)eRO; ro++ )
            {
                dRO = ro;

                int rInd = (e1-sE1)*lenRO+ro-sRO;

                // fill matrix A
                size_t col = 0;
                for ( size_t src=0; src<srcCHA; src++ )
                {
                    for ( int ke1=-kE1half; ke1<=kE1half; ke1++ )
                    {
                        for ( int kro=-kROhalf; kro<=kROhalf; kro++ )
                        {
                            if ( kro!=0 || ke1!=0 )
                            {
                                A(rInd, col++) = acsSrc(ro+kro, e1+ke1, src);
                            }
                        }
                    }
                }

                // fill matrix B
                for ( size_t dst=0; dst<dstCHA; dst++ )
                {
                    B(rInd, dst) = acsDst(dRO, dE1, dst);
                }
            }
        }

        SolveLinearSystem_Tikhonov(A, B, x, thres);

        int ind(0);
        for ( size_t src=0; src<srcCHA; src++ )
        {
            for ( int ke1=-kE1half; ke1<=kE1half; ke1++ ) 
            {
                for ( int kro=-kROhalf; kro<=kROhalf; kro++ ) 
                {
                    if ( kro!=0 || ke1!=0 )
                    {
                        for ( size_t dst=0; dst<dstCHA; dst++ )
                        {
                            ker(kro+kROhalf, ke1+kE1half, src, dst, 0, 0) = x(ind, dst);
                        }
                        ind++;
                    }
                    else
                    {
                        for ( size_t dst=0; dst<dstCHA; dst++ )
                        {
                            ker(kro+kROhalf, ke1+kE1half, src, dst, 0, 0) = 0;
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT<T>::calib(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
imageDomainKernel(const ho6DArray<T>& ker, int kRO, int kE1, int a, int b, int ro, int e1, hoNDArray<T>& kIm, bool minusI)
{
    try
    {
        int srcCHA = (int)(ker.get_size(2));
        int dstCHA = (int)(ker.get_size(3));

        GADGET_CHECK_RETURN_FALSE(kRO==ker.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kE1==ker.get_size(1));

        int kROhalf = kRO/2;
        int kE1half = kE1/2;

        // allocate image domain kernel
        kIm.create(ro, e1, srcCHA, dstCHA);

        /// fill the convolution kernels
        int convKRO = 2*kRO-1;
        int convKE1 = 2*kE1-1;

        /// fill in convolution kernel
        ho6DArray<T> convKer(convKRO, convKE1, srcCHA, dstCHA, 1, 1);
        Gadgetron::clear(&convKer);

        int kro, ke1, src, dst;
        for ( ke1=-kE1half; ke1<=kE1half; ke1++ )
        {
            for ( kro=-kROhalf; kro<=kROhalf; kro++ )
            {
                int iro = kro + kRO -1;
                int ie1 = ke1 + kE1 -1;

                for ( dst=0; dst<dstCHA; dst++ )
                {
                    for ( src=0; src<srcCHA; src++ )
                    {
                        convKer(iro, ie1, src, dst, 0, 0) = ker(kro+kROhalf, ke1+kE1half, src, dst, 0, 0);
                    }
                }
            }
        }

        hoNDArray<T> convKer2;
        ho4DArray<T> conKerMean(convKRO, convKE1, srcCHA, dstCHA);
        GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer, convKer2));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer2, conKerMean));

        // flip the kernel
        ho4DArray<T> convKerFlip(convKRO, convKE1, srcCHA, dstCHA);
        Gadgetron::clear(&convKerFlip);
        for ( ke1=0; ke1<convKE1; ke1++ )
        {
            for ( kro=0; kro<convKRO; kro++ )
            {
                for ( dst=0; dst<dstCHA; dst++ )
                {
                    for ( src=0; src<srcCHA; src++ )
                    {
                        convKerFlip( kro, ke1, src, dst) = conKerMean(convKRO-1-kro, convKE1-1-ke1, src, dst);
                    }
                }
            }
        }

        // minus I
        if ( minusI )
        {
            for ( dst=0; dst<dstCHA; dst++ )
            {
                T value = convKerFlip(kRO -1, kE1 -1, dst, dst);
                convKerFlip(kRO -1, kE1 -1, dst, dst) = value - T(1.0);
            }
        }

        GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1)) ), convKerFlip ));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad2D(convKerFlip, ro, e1, kIm));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kIm));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT<T>::imageDomainKernel(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, double thres, double overDetermineRatio, 
            int kRO, int kE1, int kE2, int a, int b, int c, hoNDArray<T>& ker)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(0)==acsDst.get_size(0));
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(1)==acsDst.get_size(1));
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(2)==acsDst.get_size(2));
        GADGET_CHECK_RETURN_FALSE(acsSrc.get_size(3)>=acsDst.get_size(3));

        size_t RO = acsSrc.get_size(0);
        size_t E1 = acsSrc.get_size(1);
        size_t E2 = acsSrc.get_size(2);
        size_t srcCHA = acsSrc.get_size(3);
        size_t dstCHA = acsDst.get_size(3);

        int kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GADGET_WARN_MSG("gtPlusSPIRIT<T>::calib3D(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        int kE1half = kE1/2;
        if ( 2*kE1half == kE1 )
        {
            GADGET_WARN_MSG("gtPlusSPIRIT<T>::calib3D(...) - 2*kE1half == kE1 " << kE1);
        }
        kE1 = 2*kE1half + 1;

        int kE2half = kE2/2;
        if ( 2*kE2half == kE2 )
        {
            GADGET_WARN_MSG("gtPlusSPIRIT<T>::calib3D(...) - 2*kE2half == kE2 " << kE2);
        }
        kE2 = 2*kE2half + 1;

        // allocate kernel
        ker.create(kRO, kE1, kE2, srcCHA, dstCHA, 1, 1, 1);

        // loop over the calibration region and assemble the equation
        // Ax = b

        size_t sRO = kROhalf;
        size_t eRO = RO - kROhalf -1;
        size_t lenRO = eRO-sRO+1;

        size_t sE1 = kE1half;
        size_t eE1 = E1 - kE1half -1;
        size_t lenE1 = eE1-sE1+1;

        size_t sE2 = kE2half;
        size_t eE2 = E2 - kE2half -1;
        size_t lenE2 = eE2-sE2+1;

        size_t colA = (kRO*kE1*kE2-1)*srcCHA;
        if ( overDetermineRatio > 1.0 )
        {
            size_t maxRowA = std::ceil(overDetermineRatio*colA);
            size_t maxROUsed = maxRowA/(lenE1*lenE2);
            if ( maxROUsed < lenRO )
            {
                // find the peak signal of acsSrc
                hoNDArray<T> acsSrc1stCha(RO, E1, E2, const_cast<T*>(acsSrc.begin()));
                hoNDArray<T> acsSrc1stChaSumE2(RO, E1, 1), acsSrc1stChaSumE2E1(RO, 1, 1);

                if ( Gadgetron::sumOver3rdDimension(acsSrc1stCha, acsSrc1stChaSumE2) )
                {
                    if ( Gadgetron::sumOver2ndDimension(acsSrc1stChaSumE2, acsSrc1stChaSumE2E1) )
                    {
                        T maxSignal;
                        size_t roInd;
                        if ( Gadgetron::maxAbsolute(acsSrc1stChaSumE2E1, maxSignal, roInd) )
                        {
                            sRO = roInd - maxROUsed/2;
                            eRO = sRO + maxROUsed - 1;
                            lenRO = eRO-sRO+1;
                            GADGET_MSG("gtPlusSPIRIT<T>::calib3D(...) - overDetermineRatio = " << overDetermineRatio << " ; RO data range used : [" << sRO << " " << eRO << "] ...");
                        }
                        else
                        {
                            GADGET_WARN_MSG("gtPlusSPIRIT<T>::calib3D(...) - overDetermineRatio is ignored ... ");
                        }
                    }
                }
                else
                {
                    GADGET_WARN_MSG("gtPlusSPIRIT<T>::calib3D(...) - overDetermineRatio is ignored ... ");
                }
            }
        }

        size_t rowA = lenRO*lenE1*lenE2;
        size_t colB = dstCHA;

        hoMatrix<T> A(rowA, colA);
        hoMatrix<T> B(rowA, colB);
        hoMatrix<T> x( A.cols(), B.cols() );

        int dRO, dE1, dE2;

        for ( int e2=(int)sE2; e2<=(int)eE2; e2++ )
        {
            dE2 = e2;

            for ( int e1=(int)sE1; e1<=(int)eE1; e1++ )
            {
                dE1 = e1;

                for ( int ro=sRO; ro<=(int)eRO; ro++ )
                {
                    dRO = ro;

                    int rInd = (e2-sE2)*lenRO*lenE1 + (e1-sE1)*lenRO + ro-sRO;

                    // fill matrix A
                    size_t col = 0;
                    for ( size_t src=0; src<srcCHA; src++ )
                    {
                        for ( int ke2=-kE2half; ke2<=kE2half; ke2++ )
                        {
                            for ( int ke1=-kE1half; ke1<=kE1half; ke1++ )
                            {
                                for ( int kro=-kROhalf; kro<=kROhalf; kro++ )
                                {
                                    if ( kro!=0 || ke1!=0 || ke2!=0 )
                                    {
                                        A(rInd, col++) = acsSrc(ro+kro, e1+ke1, e2+ke2, src);
                                    }
                                }
                            }
                        }
                    }

                    // fill matrix B
                    for ( size_t dst=0; dst<dstCHA; dst++ )
                    {
                        B(rInd, dst) = acsDst(dRO, dE1, dE2, dst);
                    }
                }
            }
        }

        SolveLinearSystem_Tikhonov(A, B, x, thres);

        int ind(0);

        std::vector<size_t> kerInd(8);
        kerInd[7] = 0;
        kerInd[6] = 0;
        kerInd[5] = 0;

        for ( size_t src=0; src<srcCHA; src++ )
        {
            kerInd[3] = src;
            for ( int ke2=-kE2half; ke2<=kE2half; ke2++ ) 
            {
                kerInd[2] = ke2+kE2half;
                for ( int ke1=-kE1half; ke1<=kE1half; ke1++ ) 
                {
                    kerInd[1] = ke1+kE1half;
                    for ( int kro=-kROhalf; kro<=kROhalf; kro++ ) 
                    {
                        kerInd[0] = kro+kROhalf;

                        if ( kro!=0 || ke1!=0 || ke2!=0 )
                        {
                            for ( size_t dst=0; dst<dstCHA; dst++ )
                            {
                                kerInd[4] = dst;
                                size_t offset = ker.calculate_offset(kerInd);
                                ker(offset) = x(ind, dst);
                            }
                            ind++;
                        }
                        else
                        {
                            for ( size_t dst=0; dst<dstCHA; dst++ )
                            {
                                kerInd[4] = dst;
                                size_t offset = ker.calculate_offset(kerInd);
                                ker(offset) = 0;
                            }
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT<T>::calib3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
kspaceDomainConvKernel3D(const hoNDArray<T>& ker, int kRO, int kE1, int kE2, int a, int b, int c, ho5DArray<T>& convKerFlip, bool minusI, bool ROis3rdDim)
{
    try
    {
        int srcCHA = (int)(ker.get_size(3));
        int dstCHA = (int)(ker.get_size(4));

        GADGET_CHECK_RETURN_FALSE(kRO==ker.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kE1==ker.get_size(1));
        GADGET_CHECK_RETURN_FALSE(kE2==ker.get_size(2));

        int kROhalf = kRO/2;
        int kE1half = kE1/2;
        int kE2half = kE2/2;

        /// fill the convolution kernels
        int convKRO = 2*kRO-1;
        int convKE1 = 2*kE1-1;
        int convKE2 = 2*kE2-1;

        /// fill in convolution kernel
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - convert to conv kernel ... "));

        hoNDArray<T> convKer(convKRO, convKE1, convKE2, srcCHA, dstCHA, 1, 1, 1);
        Gadgetron::clear(&convKer);

        int kro, ke1, ke2, src, dst;
        std::vector<size_t> kerInd(8), convKerInd(8);

        kerInd[7] = 0;
        convKerInd[7] = 0;

        kerInd[6] = 0;
        convKerInd[6] = 0;

        kerInd[5] = 0;
        convKerInd[5] = 0;

        for ( ke2=-kE2half; ke2<=kE2half; ke2++ )
        {
            kerInd[2] = ke2+kE2half;

            for ( ke1=-kE1half; ke1<=kE1half; ke1++ )
            {
                kerInd[1] = ke1+kE1half;

                for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                {
                    int iro = kro + kRO -1;
                    int ie1 = ke1 + kE1 -1;
                    int ie2 = ke2 + kE2 -1;

                    kerInd[0] = kro+kROhalf;

                    convKerInd[0] = iro;
                    convKerInd[1] = ie1;
                    convKerInd[2] = ie2;

                    for ( dst=0; dst<dstCHA; dst++ )
                    {
                        kerInd[4] = dst;
                        convKerInd[4] = dst;

                        for ( src=0; src<srcCHA; src++ )
                        {
                            kerInd[3] = src;
                            convKerInd[3] = src;

                            size_t offsetKer = ker.calculate_offset(kerInd);
                            size_t offsetConvKer = convKer.calculate_offset(convKerInd);

                            convKer(offsetConvKer) = ker(offsetKer);
                        }
                    }
                }
            }
        }

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - sum over output dimensions ... "));
        hoNDArray<T> convKer2, convKer3;
        ho5DArray<T> convKernMean(convKRO, convKE1, convKE2, srcCHA, dstCHA);
        GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer, convKer2));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer2, convKer3));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer3, convKernMean));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - flip along dimensions ... "));

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, convKernMean, "convKernMean");

        // flip the kernel
        if ( ROis3rdDim ) // E1, E2, RO
        {
            convKerFlip.createArray(convKE1, convKE2, convKRO, srcCHA, dstCHA);
            Gadgetron::clear(&convKerFlip);

            for ( ke2=0; ke2<convKE2; ke2++ )
            {
                for ( ke1=0; ke1<convKE1; ke1++ )
                {
                    for ( kro=0; kro<convKRO; kro++ )
                    {
                        for ( dst=0; dst<dstCHA; dst++ )
                        {
                            for ( src=0; src<srcCHA; src++ )
                            {
                                T value = convKernMean(convKRO-1-kro, convKE1-1-ke1, convKE2-1-ke2, src, dst);
                                convKerFlip(ke1, ke2, kro, src, dst) = value;
                            }
                        }
                    }
                }
            }
            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, convKerFlip, "convKerFlip");

            // minus I
            if ( minusI )
            {
                for ( dst=0; dst<dstCHA; dst++ )
                {
                    T value = convKerFlip(kE1 -1, kE2 -1, kRO -1, dst, dst);
                    convKerFlip(kE1 -1, kE2 -1, kRO -1, dst, dst) = value - T(1.0);
                }

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, convKerFlip, "convKerFlip_minusI");
            }
        }
        else
        {
            // RO, E1, E2
            convKerFlip.createArray(convKRO, convKE1, convKE2, srcCHA, dstCHA);
            Gadgetron::clear(&convKerFlip);

            for ( ke2=0; ke2<convKE2; ke2++ )
            {
                for ( ke1=0; ke1<convKE1; ke1++ )
                {
                    for ( kro=0; kro<convKRO; kro++ )
                    {
                        for ( dst=0; dst<dstCHA; dst++ )
                        {
                            for ( src=0; src<srcCHA; src++ )
                            {
                                T value = convKernMean(convKRO-1-kro, convKE1-1-ke1, convKE2-1-ke2, src, dst);
                                convKerFlip(kro, ke1, ke2, src, dst) = value;
                            }
                        }
                    }
                }
            }
            GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, convKerFlip, "convKerFlip");

            // minus I
            if ( minusI )
            {
                for ( dst=0; dst<dstCHA; dst++ )
                {
                    T value = convKerFlip(kRO -1, kE1 -1, kE2 -1, dst, dst);
                    convKerFlip(kRO -1, kE1 -1, kE2 -1, dst, dst) = value - T(1.0);
                }

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, convKerFlip, "convKerFlip_minusI");
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT<T>::kspaceDomainConvKernel3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
imageDomainKernel3D(const hoNDArray<T>& ker, int kRO, int kE1, int kE2, int a, int b, int c, int ro, int e1, int e2, hoNDArray<T>& kIm, bool minusI)
{
    try
    {
        int srcCHA = (int)(ker.get_size(3));
        int dstCHA = (int)(ker.get_size(4));

        GADGET_CHECK_RETURN_FALSE(kRO==ker.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kE1==ker.get_size(1));
        GADGET_CHECK_RETURN_FALSE(kE2==ker.get_size(2));

        // allocate image domain kernel
        kIm.create(e1, e2, ro, srcCHA, dstCHA);

        bool ROat3rdDim = true;
        ho5DArray<T> convKerFlip;
        GADGET_CHECK_RETURN_FALSE(this->kspaceDomainConvKernel3D(ker, kRO, kE1,  kE2, a, b, c, convKerFlip, minusI, ROat3rdDim));

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - SNR unit scaling ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1*e2)) ), convKerFlip ));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, convKerFlip, "convKerFlip_scal");

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - zero padding ... "));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3DNoPresetZeros(convKerFlip, e1, e2, ro, kIm));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kIm, "convKerFlip_scal_zeropadded");

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - conver to image domain ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kIm));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT<T>::imageDomainKernel3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
imageDomainKernelRO3D(const hoNDArray<T>& ker, int kRO, int kE1, int kE2, int a, int b, int c, int ro, hoNDArray<T>& kImRO, bool minusI)
{
    try
    {
        int srcCHA = (int)(ker.get_size(3));
        int dstCHA = (int)(ker.get_size(4));

        GADGET_CHECK_RETURN_FALSE(kRO==ker.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kE1==ker.get_size(1));
        GADGET_CHECK_RETURN_FALSE(kE2==ker.get_size(2));

        bool ROat3rdDim = false;
        ho5DArray<T> convKerFlip;
        GADGET_CHECK_RETURN_FALSE(this->kspaceDomainConvKernel3D(ker, kRO, kE1,  kE2, a, b, c, convKerFlip, minusI, ROat3rdDim));

        // allocate image domain kernel
        size_t kConvE1 = convKerFlip.get_size(1);
        size_t kConvE2 = convKerFlip.get_size(2);

        kImRO.create(kConvE1, kConvE2, ro, srcCHA, dstCHA);

        hoNDArray<T> kImROTemp(ro, kConvE1, kConvE2, srcCHA, dstCHA);
        Gadgetron::clear(kImROTemp);

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - SNR unit scaling ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro)) ), convKerFlip ));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, convKerFlip, "convKerFlip_scal_RO");

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - zero padding only for RO ... "));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3DNoPresetZeros(convKerFlip, ro, kConvE1, kConvE2, kImROTemp));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kImROTemp, "convKerFlip_scal_RO_zeropadded");

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - conver to image domain only for RO ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(kImROTemp));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - permute kernel dimensions to be [kE1 kE2 RO ...]  ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteROTo3rdDimensionFor3DRecon(kImROTemp, kImRO));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT<T>::imageDomainKernel3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
imageDomainKernelE1E2RO(const hoNDArray<T>& kImRO, int e1, int e2, hoNDArray<T>& kImE1E2RO)
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

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - SNR unit scaling for E1 and E2 ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(e1*e2)) ), kImROScaled ));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kImROScaled, "kImROScaledE1E2");

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - zero padding for E1 and E2 ... "));
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3DNoPresetZeros(kImROScaled, e1, e2, dimR[2], kImE1E2RO));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, kImE1E2RO, "kImE1E2RO_zeropadded_E1E2");

        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.start("spirit 3D calibration - conver to image domain for E1 and E2 ... "));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kImE1E2RO));
        GADGET_CHECK_PERFORM(performTiming_, gt_timer3_.stop());
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT<T>::imageDomainKernelE1E2RO(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
imageDomainAdjointKernel(const hoNDArray<T>& kIm, hoNDArray<T>& adjkIm)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteLastTwoDimensions(kIm, adjkIm));

        GADGET_CHECK_RETURN_FALSE(Gadgetron::conjugate(adjkIm, adjkIm));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT<T>::imageDomainAdjointKernel(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::AdjointForwardKernel(const hoNDArray<T>& kImS2D, const hoNDArray<T>& kImD2S, hoNDArray<T>& kIm)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dimS2D = kImS2D.get_dimensions();

        size_t NDim = kImS2D.get_number_of_dimensions();

        long long srcCHA = (*dimS2D)[NDim-2];
        long long dstCHA = (*dimS2D)[NDim-1];

        GADGET_CHECK_RETURN_FALSE(kImD2S.get_number_of_dimensions()==NDim);
        GADGET_CHECK_RETURN_FALSE(kImD2S.get_number_of_elements()==kImS2D.get_number_of_elements());

        std::vector<size_t> dimRes(*dimS2D);
        dimRes[NDim-2] = dstCHA;

        kIm.create(&dimRes);
        Gadgetron::clear(&kIm);

        size_t N = kImS2D.get_number_of_elements()/srcCHA/dstCHA;

        long long d;
        #ifdef GCC_OLD_FLAG
            #pragma omp parallel default(none) private(d) shared(N, dstCHA, srcCHA) num_threads(dstCHA) if (dstCHA > 4)
        #else
            #pragma omp parallel default(none) private(d) shared(N, dstCHA, srcCHA, kIm, kImS2D, kImD2S) num_threads(dstCHA) if (dstCHA > 4)
        #endif
        {
            hoNDArray<T> ker(N);

            #pragma omp for
            for ( d=0; d<dstCHA; d++ )
            {
                for ( size_t dprime=0; dprime<dstCHA; dprime++ )
                {
                    hoNDArray<T> dKer(N, kIm.begin()+d*N+dprime*N*dstCHA);

                    for ( size_t s=0; s<srcCHA; s++ )
                    {
                        hoNDArray<T> kerS2D(N, const_cast<T*>(kImS2D.begin())+s*N+dprime*N*srcCHA);
                        hoNDArray<T> kerD2S(N, const_cast<T*>(kImD2S.begin())+d*N+s*N*dstCHA);

                        Gadgetron::multiply(kerS2D, kerD2S, ker);
                        Gadgetron::add(dKer, ker, dKer);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT<T>::AdjointForwardKernel(...) ... ");
        return false;
    }

    return true;
}

}}
