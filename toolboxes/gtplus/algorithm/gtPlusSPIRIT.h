
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

    gtPlusSPIRIT() : calib_use_gpu_(true), BaseClass() {}
    virtual ~gtPlusSPIRIT() {}

    virtual void printInfo(std::ostream& os);

    // SPIRIT calibration for 2D case
    // acsSrc : [RO E1 srcCHA]
    // acsDst : [RO E1 dstCHA]
    // ker : [kRO kE1 srcCHA dstCHA oRO oE1]
    bool calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, 
            size_t kRO, size_t kE1, size_t oRO, size_t oE1, ho6DArray<T>& ker);

    // image domain kernel for 2D kernel
    // kIm: image domain kernel [RO E1 srcCHA dstCHA]
    // if minusI==true, compute image domain G-I kernel
    bool imageDomainKernel(const ho6DArray<T>& ker, size_t kRO, size_t kE1, 
        size_t oRO, size_t oE1, size_t ro, size_t e1, hoNDArray<T>& kIm, bool minusI=false);

    // SPIRIT calibration for 3D case
    // acsSrc : [RO E1 E2 srcCHA]
    // acsDst : [RO E1 E2 dstCHA]
    // ker : [kRO kE1 kE2 srcCHA dstCHA oRO oE1 oE2]
    // overDetermineRatio : over determine ratio of calib matrix, if < 1, all data are used
    bool calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, double thres, double overDetermineRatio, 
            size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray<T>& ker);

    // convert the calibrated kernel to the convlution kernel in kspace
    // if ROis3rdDim == true, the kernel dimension is [E1 E2 RO], otherwise [RO E1 E2]
    bool kspaceDomainConvKernel3D(const hoNDArray<T>& ker, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, ho5DArray<T>& convKerFlip, bool minusI=true, bool ROis3rdDim=true);

    // image domain kernel for 3D kernel
    // kIm: image domain kernel [E1 E2 RO srcCHA dstCHA]
    // if minusI==true, compute image domain G-I kernel
    bool imageDomainKernel3D(const hoNDArray<T>& ker, size_t kRO, size_t kE1, size_t kE2, 
        size_t oRO, size_t oE1, size_t oE2, size_t ro, size_t e1, size_t e2, hoNDArray<T>& kIm, bool minusI=false);

    // image domain kernel for 3D kernel, only RO direction is converted to image domain
    // E1 and E2 stays in the kspace domain
    // kImRO: kspace-image hybrid kernel [convE1 convE2 RO srcCHA dstCHA]
    bool imageDomainKernelRO3D(const hoNDArray<T>& ker, size_t kRO, size_t kE1, size_t kE2, 
        size_t oRO, size_t oE1, size_t oE2, size_t ro, hoNDArray<T>& kImRO, bool minusI=false);

    // image domain kernel for 3D kernel, E1 and E2 directions are converted to image domain
    // kImRO : kspace-image hybrid kernel where first two dimensions are E1 and E2 and in kspace
    bool imageDomainKernelE1E2RO(const hoNDArray<T>& kImRO, size_t e1, size_t e2, hoNDArray<T>& kImE1E2RO);

    // compute the image domain adjoint kernel
    bool imageDomainAdjointKernel(const hoNDArray<T>& kIm, hoNDArray<T>& adjkIm);

    // compute the (G-I)'*(G-I)
    bool AdjointForwardKernel(const hoNDArray<T>& kImS2D, const hoNDArray<T>& kImD2S, hoNDArray<T>& kIm);

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
void gtPlusSPIRIT<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT reconstruction ------------------" << endl;
    os << "Implementation of SPIRIT algorithms for ISMRMRD package" << endl;
    os << "Both 2D and 3D version are implemented" << endl;
    os << "Algorithms are published at:" << endl;
    os << "Lustig, M. and Pauly, J. M. (2010), SPIRiT: Iterative self-consistent parallel imaging reconstruction from arbitrary k-space. Magn Reson Med, 64: 457�471. doi: 10.1002/mrm.22428" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, 
            size_t kRO, size_t kE1, size_t oRO, size_t oE1, ho6DArray<T>& ker)
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

        long long kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        long long kE1half = kE1/2;
        if ( 2*kE1half == kE1 )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib(...) - 2*kE1half == kE1 " << kE1);
        }
        kE1 = 2*kE1half + 1;

        if ( oRO > kRO ) oRO = kRO;
        if ( oE1 > kE1 ) oE1 = kE1;

        long long oROhalf = oRO/2;
        if ( 2*oROhalf == oRO )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib(...) - 2*oROhalf == oRO " << oRO);
        }
        oRO = 2*oROhalf + 1;

        long long oE1half = oE1/2;
        if ( 2*oE1half == oE1 )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib(...) - 2*oE1half == oE1 " << oE1);
        }
        oE1 = 2*oE1half + 1;

        // allocate kernel
        GADGET_CHECK_RETURN_FALSE(ker.createArray(kRO, kE1, srcCHA, dstCHA, oRO, oE1));

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

        bool useGPU = (typeid(typename realType<T>::Type)==typeid(float) && calib_use_gpu_);
        //if ( useGPU )
        //{
        //    GDEBUG_STREAM("spirit 2D - calling GPU kernel estimation ... "); 
        //}

        const T* pAcsSrc = acsSrc.begin();

        #pragma omp parallel default(none) shared(RO, E1, sRO, eRO, sE1, eE1, oRO, oE1, lenRO, lenE1, rowA, colA, colB, kRO, kE1, kROhalf, kE1half, oROhalf, oE1half, pAcsSrc, acsSrc, acsDst, srcCHA, dstCHA, thres, ker, useGPU, std::cout) num_threads( (int)(oRO*oE1) ) if (oRO*oE1>=3)
        {
            hoMatrix<T> A(rowA, colA);
            T* pA = A.begin();

            hoMatrix<T> B(rowA, colB);
            T* pB = B.begin();

            hoMatrix<T> x( A.cols(), B.cols() );

            long long kInd = 0;
            #pragma omp for
            for ( kInd=0; kInd<(long long)(oRO*oE1); kInd++ )
            {
                long long oe1 = kInd/oRO;
                long long oro = kInd - oe1*oRO;

                oe1 -=oE1half;
                oro -=oROhalf;

                long long dRO, dE1;

                for ( long long e1=(long long)sE1; e1<=(long long)eE1; e1++ )
                {
                    dE1 = e1 + oe1;

                    for ( long long ro=sRO; ro<=(long long)eRO; ro++ )
                    {
                        dRO = ro + oro;

                        long long rInd = (e1-sE1)*lenRO+ro-sRO;

                        // fill matrix A
                        size_t col = 0;
                        for ( size_t src=0; src<srcCHA; src++ )
                        {
                            for ( long long ke1=-kE1half; ke1<=kE1half; ke1++ )
                            {
                                for ( long long kro=-kROhalf; kro<=kROhalf; kro++ )
                                {
                                    if ( kro!=oro || ke1!=oe1 )
                                    {
                                        //A(rInd, col++) = acsSrc(ro+kro, e1+ke1, src);
                                        // pA[rInd + col*rowA] = acsSrc(ro+kro, e1+ke1, src);
                                        pA[rInd + col*rowA] = pAcsSrc[ro+kro + (e1+ke1)*RO + src*RO*E1];
                                        col++;
                                    }
                                }
                            }
                        }

                        // fill matrix B
                        for ( size_t dst=0; dst<dstCHA; dst++ )
                        {
                            //B(rInd, dst) = acsDst(dRO, dE1, dst);
                            pB[rInd+dst*rowA] = acsDst(dRO, dE1, dst);
                        }
                    }
                }

                // GADGET_CHECK_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));

                //if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - solve linear system ... "); }
                //#ifdef USE_CUDA
                //    // go to device
                //    try
                //    {
                //        if ( useGPU )
                //        {
                //            hoNDArray<float_complext> A_tmp(A.get_dimensions(), reinterpret_cast<float_complext*>(A.begin()));
                //            hoNDArray<float_complext> B_tmp(B.get_dimensions(), reinterpret_cast<float_complext*>(B.begin()));

                //            int ret(0);
                //            boost::shared_ptr< hoNDArray<complext<float> > > host_x;

                //            #pragma omp critical(inverse_spirit)
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
                //                memcpy(x.begin(), host_x->begin(), host_x->get_number_of_bytes());
                //            }
                //        }
                //        else
                //        {
                //            GWARN_STREAM("GPU inverse_clib_matrix is only available for single-precision, calling the CPU version ... ");
                //            SolveLinearSystem_Tikhonov(A, B, x, thres);
                //        }
                //    }
                //    catch(...)
                //    {
                //        GERROR_STREAM("failed in GPU inverse_clib_matrix for grappa, calling the CPU version ... ");
                //        SolveLinearSystem_Tikhonov(A, B, x, thres);
                //    }
                //#else
                    SolveLinearSystem_Tikhonov(A, B, x, thres);
                //#endif // USE_CUDA
                //if ( performTiming_ ) { gt_timer3_.stop(); }

                //SolveLinearSystem_Tikhonov(A, B, x, thres);

                long long ind(0);
                for ( size_t src=0; src<srcCHA; src++ )
                {
                    for ( long long ke1=-kE1half; ke1<=kE1half; ke1++ ) 
                    {
                        for ( long long kro=-kROhalf; kro<=kROhalf; kro++ ) 
                        {
                            if ( kro!=oro || ke1!=oe1 )
                            {
                                for ( size_t dst=0; dst<dstCHA; dst++ )
                                {
                                    ker(kro+kROhalf, ke1+kE1half, src, dst, oro+oROhalf, oe1+oE1half) = x(ind, dst);
                                }
                                ind++;
                            }
                            else
                            {
                                for ( size_t dst=0; dst<dstCHA; dst++ )
                                {
                                    ker(kro+kROhalf, ke1+kE1half, src, dst, oro+oROhalf, oe1+oE1half) = 0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRIT<T>::calib(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
imageDomainKernel(const ho6DArray<T>& ker, size_t kRO, size_t kE1, size_t oRO, size_t oE1, size_t ro, size_t e1, hoNDArray<T>& kIm, bool minusI)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(2));
        long long dstCHA = (long long)(ker.get_size(3));

        GADGET_CHECK_RETURN_FALSE(kRO==ker.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kE1==ker.get_size(1));
        GADGET_CHECK_RETURN_FALSE(oRO==ker.get_size(4));
        GADGET_CHECK_RETURN_FALSE(oE1==ker.get_size(5));

        long long kROhalf = kRO/2;
        long long kE1half = kE1/2;
        long long oROhalf = oRO/2;
        long long oE1half = oE1/2;

        // allocate image domain kernel
        kIm.create(ro, e1, srcCHA, dstCHA);

        /// fill the convolution kernels
        long long convKRO = 2*kRO-1;
        long long convKE1 = 2*kE1-1;

        /// fill in convolution kernel
        ho6DArray<T> convKer(convKRO, convKE1, srcCHA, dstCHA, oRO, oE1);
        Gadgetron::clear(&convKer);

        long long oro, oe1, kro, ke1, src, dst;
        for ( oe1=-oE1half; oe1<=oE1half; oe1++ )
        {
            for ( oro=-oROhalf; oro<=oROhalf; oro++ )
            {
                for ( ke1=-kE1half; ke1<=kE1half; ke1++ )
                {
                    for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                    {
                        long long iro = kro - oro + kRO -1;
                        long long ie1 = ke1 - oe1 + kE1 -1;

                        for ( dst=0; dst<dstCHA; dst++ )
                        {
                            for ( src=0; src<srcCHA; src++ )
                            {
                                convKer(iro, ie1, src, dst, oro+oROhalf, oe1+oE1half) = ker(kro+kROhalf, ke1+kE1half, src, dst, oro+oROhalf, oe1+oE1half);
                            }
                        }
                    }
                }
            }
        }

        hoNDArray<T> convKer2;
        hoNDArray<T> conKerMean(convKRO, convKE1, srcCHA, dstCHA, 1, 1);

        //GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer, convKer2));
        //GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer2, conKerMean));

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(convKer, convKer2, 5));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(convKer2, conKerMean, 4));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)(1.0/(oRO*oE1)), conKerMean) );

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
                        convKerFlip( kro, ke1, src, dst) = conKerMean(convKRO-1-kro, convKE1-1-ke1, src, dst, 0, 0);
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

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1)) ), convKerFlip ));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(ro, e1, &convKerFlip, &kIm));
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kIm);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRIT<T>::imageDomainKernel(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, double thres, double overDetermineRatio, 
            size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray<T>& ker)
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

        long long kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib3D(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        long long kE1half = kE1/2;
        if ( 2*kE1half == kE1 )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib3D(...) - 2*kE1half == kE1 " << kE1);
        }
        kE1 = 2*kE1half + 1;

        long long kE2half = kE2/2;
        if ( 2*kE2half == kE2 )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib3D(...) - 2*kE2half == kE2 " << kE2);
        }
        kE2 = 2*kE2half + 1;

        if ( oRO > kRO ) oRO = kRO;
        if ( oE1 > kE1 ) oE1 = kE1;
        if ( oE2 > kE2 ) oE2 = kE2;

        long long oROhalf = oRO/2;
        if ( 2*oROhalf == oRO )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib3D(...) - 2*oROhalf == oRO " << oRO);
        }
        oRO = 2*oROhalf + 1;

        long long oE1half = oE1/2;
        if ( 2*oE1half == oE1 )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib3D(...) - 2*oE1half == oE1 " << oE1);
        }
        oE1 = 2*oE1half + 1;

        long long oE2half = oE2/2;
        if ( 2*oE2half == oE2 )
        {
            GWARN_STREAM("gtPlusSPIRIT<T>::calib3D(...) - 2*oE2half == oE2 " << oE2);
        }
        oE2 = 2*oE2half + 1;

        // allocate kernel
        ker.create(kRO, kE1, kE2, srcCHA, dstCHA, oRO, oE1, oE2);

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
                        GDEBUG_STREAM("gtPlusSPIRIT<T>::calib3D(...) - overDetermineRatio = " << overDetermineRatio << " ; RO data range used : [" << sRO << " " << eRO << "] ...");
                    }
                    catch(...)
                    {
                        GWARN_STREAM("gtPlusSPIRIT<T>::calib3D(...) - overDetermineRatio is ignored ... ");
                        throw;
                    }
                }
                catch (...)
                {
                    GWARN_STREAM("gtPlusSPIRIT<T>::calib3D(...) - overDetermineRatio is ignored ... ");
                    throw;
                }
            }
        }

        size_t rowA = lenRO*lenE1*lenE2;
        size_t colB = dstCHA;

        bool useGPU = (typeid(typename realType<T>::Type)==typeid(float) && calib_use_gpu_);
        if ( useGPU )
        {
            GDEBUG_STREAM("spirit 3D - calling GPU kernel estimation ... ");
        }

        #pragma omp parallel default(none) shared(sRO, eRO, sE1, eE1, sE2, eE2, oRO, oE1, oE2, lenRO, lenE1, lenE2, rowA, colA, colB, kROhalf, kE1half, kE2half, oROhalf, oE1half, oE2half, acsSrc, acsDst, srcCHA, dstCHA, thres, ker, useGPU, std::cout) num_threads( (int)(oRO*oE1*oE2) ) if (oRO*oE1*oE2>=3 && oRO*oE1*oE2<9)
        {
            hoMatrix<T> A(rowA, colA);
            hoMatrix<T> B(rowA, colB);
            hoMatrix<T> x( A.cols(), B.cols() );

            // hoNDArrayMemoryManaged<T> A_mem(colA, rowA, gtPlus_mem_manager_);
            //hoNDArray<T> A_mem(colA, rowA);
            //A.createMatrix( rowA, colA, A_mem.begin() );

            // hoNDArrayMemoryManaged<T> B_mem(colB, rowA, gtPlus_mem_manager_);
            // hoNDArray<T> B_mem(colB, rowA);
            // B.createMatrix( A.rows(), colB, B_mem.begin() );

            T* pA = A.begin();
            T* pB = B.begin();

            long long kInd = 0;
            #pragma omp for
            for ( kInd=0; kInd<(long long)(oRO*oE1*oE2); kInd++ )
            {
                long long oe2 = kInd/(oRO*oE1);
                long long oe1 = kInd - oe2*oRO*oE1;
                oe1 /= oRO;
                long long oro = kInd - oe2*oRO*oE1 - oe1*oRO;

                oe2 -=oE2half;
                oe1 -=oE1half;
                oro -=oROhalf;

                long long dRO, dE1, dE2;

                for ( long long e2=(long long)sE2; e2<=(long long)eE2; e2++ )
                {
                    dE2 = e2 + oe2;

                    for ( long long e1=(long long)sE1; e1<=(long long)eE1; e1++ )
                    {
                        dE1 = e1 + oe1;

                        for ( long long ro=sRO; ro<=(long long)eRO; ro++ )
                        {
                            dRO = ro + oro;

                            long long rInd = (e2-sE2)*lenRO*lenE1 + (e1-sE1)*lenRO + ro-sRO;

                            // fill matrix A
                            size_t col = 0;
                            for ( size_t src=0; src<srcCHA; src++ )
                            {
                                for ( long long ke2=-kE2half; ke2<=kE2half; ke2++ )
                                {
                                    for ( long long ke1=-kE1half; ke1<=kE1half; ke1++ )
                                    {
                                        for ( long long kro=-kROhalf; kro<=kROhalf; kro++ )
                                        {
                                            if ( kro!=oro || ke1!=oe1 || ke2!=oe2 )
                                            {
                                                //A(rInd, col++) = acsSrc(ro+kro, e1+ke1, e2+ke2, src);
                                                pA[rInd+col*rowA] = acsSrc(ro+kro, e1+ke1, e2+ke2, src);
                                                col++;
                                            }
                                        }
                                    }
                                }
                            }

                            // fill matrix B
                            for ( size_t dst=0; dst<dstCHA; dst++ )
                            {
                                //B(rInd, dst) = acsDst(dRO, dE1, dE2, dst);
                                pB[rInd+dst*rowA] = acsDst(dRO, dE1, dE2, dst);
                            }
                        }
                    }
                }

                //GADGET_CHECK_RETURN_FALSE(SolveLinearSystem_Tikhonov(A, B, x, thres));

                //if ( performTiming_ ) { gt_timer3_.start("grappa 3D calibration - solve linear system ... "); }
                //#ifdef USE_CUDA
                //    // go to device
                //    try
                //    {
                //        if ( useGPU )
                //        {
                //            hoNDArray<float_complext> A_tmp(A.get_dimensions(), reinterpret_cast<float_complext*>(A.begin()));
                //            hoNDArray<float_complext> B_tmp(B.get_dimensions(), reinterpret_cast<float_complext*>(B.begin()));

                //            int ret(0);
                //            boost::shared_ptr< hoNDArray<complext<float> > > host_x;
                //            #pragma omp critical(inverse_spirit3D)
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
                //            GWARN_STREAM("GPU inverse_clib_matrix is only available for single-precision, calling the CPU version ... ");
                //            SolveLinearSystem_Tikhonov(A, B, x, thres);
                //        }
                //    }
                //    catch(...)
                //    {
                //        GERROR_STREAM("failed in GPU inverse_clib_matrix for grappa, calling the CPU version ... ");
                //        SolveLinearSystem_Tikhonov(A, B, x, thres);
                //    }
                //#else
                    SolveLinearSystem_Tikhonov(A, B, x, thres);
                //#endif // USE_CUDA
                //if ( performTiming_ ) { gt_timer3_.stop(); }

                // SolveLinearSystem_Tikhonov(A, B, x, thres);

                long long ind(0);

                std::vector<size_t> kerInd(8);
                kerInd[7] = oe2+oE2half;
                kerInd[6] = oe1+oE1half;
                kerInd[5] = oro+oROhalf;

                for ( size_t src=0; src<srcCHA; src++ )
                {
                    kerInd[3] = src;
                    for ( long long ke2=-kE2half; ke2<=kE2half; ke2++ ) 
                    {
                        kerInd[2] = ke2+kE2half;
                        for ( long long ke1=-kE1half; ke1<=kE1half; ke1++ ) 
                        {
                            kerInd[1] = ke1+kE1half;
                            for ( long long kro=-kROhalf; kro<=kROhalf; kro++ ) 
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
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRIT<T>::calib3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
kspaceDomainConvKernel3D(const hoNDArray<T>& ker, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, ho5DArray<T>& convKerFlip, bool minusI, bool ROis3rdDim)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(3));
        long long dstCHA = (long long)(ker.get_size(4));

        GADGET_CHECK_RETURN_FALSE(kRO==ker.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kE1==ker.get_size(1));
        GADGET_CHECK_RETURN_FALSE(kE2==ker.get_size(2));
        GADGET_CHECK_RETURN_FALSE(oRO==ker.get_size(5));
        GADGET_CHECK_RETURN_FALSE(oE1==ker.get_size(6));
        GADGET_CHECK_RETURN_FALSE(oE2==ker.get_size(7));

        long long kROhalf = kRO/2;
        long long kE1half = kE1/2;
        long long kE2half = kE2/2;
        long long oROhalf = oRO/2;
        long long oE1half = oE1/2;
        long long oE2half = oE2/2;

        /// fill the convolution kernels
        long long convKRO = 2*kRO-1;
        long long convKE1 = 2*kE1-1;
        long long convKE2 = 2*kE2-1;

        /// fill in convolution kernel
        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - convert to conv kernel ... "); }

        hoNDArray<T> convKer(convKRO, convKE1, convKE2, srcCHA, dstCHA, oRO, oE1, oE2);
        Gadgetron::clear(&convKer);

        long long oro, oe1, oe2, kro, ke1, ke2, src, dst;
        std::vector<size_t> kerInd(8), convKerInd(8);
        for ( oe2=-oE2half; oe2<=oE2half; oe2++ )
        {
            kerInd[7] = oe2+oE2half;
            convKerInd[7] = oe2+oE2half;

            for ( oe1=-oE1half; oe1<=oE1half; oe1++ )
            {
                kerInd[6] = oe1+oE1half;
                convKerInd[6] = oe1+oE1half;

                for ( oro=-oROhalf; oro<=oROhalf; oro++ )
                {
                    kerInd[5] = oro+oROhalf;
                    convKerInd[5] = oro+oROhalf;

                    for ( ke2=-kE2half; ke2<=kE2half; ke2++ )
                    {
                        kerInd[2] = ke2+kE2half;

                        for ( ke1=-kE1half; ke1<=kE1half; ke1++ )
                        {
                            kerInd[1] = ke1+kE1half;

                            for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                            {
                                long long iro = kro - oro + kRO -1;
                                long long ie1 = ke1 - oe1 + kE1 -1;
                                long long ie2 = ke2 - oe2 + kE2 -1;

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
                }
            }
        }
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - sum over output dimensions ... "); }
        hoNDArray<T> convKer2, convKer3;

        //ho5DArray<T> convKernMean(convKRO, convKE1, convKE2, srcCHA, dstCHA);
        //GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer, convKer2));
        //GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer2, convKer3));
        //GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(convKer3, convKernMean));

        hoNDArray<T> convKernMean(convKRO, convKE1, convKE2, srcCHA, dstCHA, 1, 1, 1);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(convKer, convKer2, 7));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(convKer2, convKer3, 6));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(convKer3, convKernMean, 5));

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)(1.0/(oRO*oE1*oE2)), convKernMean) );
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - flip along dimensions ... "); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(convKernMean, debugFolder_+"convKernMean"); }

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
                                T value = convKernMean(convKRO-1-kro, convKE1-1-ke1, convKE2-1-ke2, src, dst, 0, 0, 0);
                                convKerFlip(ke1, ke2, kro, src, dst) = value;
                            }
                        }
                    }
                }
            }
            if ( performTiming_ ) { gt_timer3_.stop(); }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(convKerFlip, debugFolder_+"convKerFlip"); }

            // minus I
            if ( minusI )
            {
                for ( dst=0; dst<dstCHA; dst++ )
                {
                    T value = convKerFlip(kE1 -1, kE2 -1, kRO -1, dst, dst);
                    convKerFlip(kE1 -1, kE2 -1, kRO -1, dst, dst) = value - T(1.0);
                }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(convKerFlip, debugFolder_+"convKerFlip_minusI"); }
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
            if ( performTiming_ ) { gt_timer3_.stop(); }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(convKerFlip, debugFolder_+"convKerFlip"); }

            // minus I
            if ( minusI )
            {
                for ( dst=0; dst<dstCHA; dst++ )
                {
                    T value = convKerFlip(kRO -1, kE1 -1, kE2 -1, dst, dst);
                    convKerFlip(kRO -1, kE1 -1, kE2 -1, dst, dst) = value - T(1.0);
                }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(convKerFlip, debugFolder_+"convKerFlip_minusI"); }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRIT<T>::kspaceDomainConvKernel3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
imageDomainKernel3D(const hoNDArray<T>& ker, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, size_t ro, size_t e1, size_t e2, hoNDArray<T>& kIm, bool minusI)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(3));
        long long dstCHA = (long long)(ker.get_size(4));

        GADGET_CHECK_RETURN_FALSE(kRO==ker.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kE1==ker.get_size(1));
        GADGET_CHECK_RETURN_FALSE(kE2==ker.get_size(2));
        GADGET_CHECK_RETURN_FALSE(oRO==ker.get_size(5));
        GADGET_CHECK_RETURN_FALSE(oE1==ker.get_size(6));
        GADGET_CHECK_RETURN_FALSE(oE2==ker.get_size(7));

        // allocate image domain kernel
        kIm.create(e1, e2, ro, srcCHA, dstCHA);

        bool ROat3rdDim = true;
        ho5DArray<T> convKerFlip;
        GADGET_CHECK_RETURN_FALSE(this->kspaceDomainConvKernel3D(ker, kRO, kE1,  kE2, oRO, oE1, oE2, convKerFlip, minusI, ROat3rdDim));

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - SNR unit scaling ... "); }
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro*e1*e2)) ), convKerFlip ));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(convKerFlip, debugFolder_+"convKerFlip_scal"); }

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - zero padding ... "); }
        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3D(convKerFlip, e1, e2, ro, kIm));
        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3DNoPresetZeros(convKerFlip, e1, e2, ro, kIm));
        // GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(convKerFlip, e1, e2, ro, kIm, false));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(e1, e2, ro, &convKerFlip, &kIm, false));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kIm, debugFolder_+"convKerFlip_scal_zeropadded"); }

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - conver to image domain ... "); }
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kIm);
        if ( performTiming_ ) { gt_timer3_.stop(); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRIT<T>::imageDomainKernel3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
imageDomainKernelRO3D(const hoNDArray<T>& ker, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, size_t ro, hoNDArray<T>& kImRO, bool minusI)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(3));
        long long dstCHA = (long long)(ker.get_size(4));

        GADGET_CHECK_RETURN_FALSE(kRO==ker.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kE1==ker.get_size(1));
        GADGET_CHECK_RETURN_FALSE(kE2==ker.get_size(2));
        GADGET_CHECK_RETURN_FALSE(oRO==ker.get_size(5));
        GADGET_CHECK_RETURN_FALSE(oE1==ker.get_size(6));
        GADGET_CHECK_RETURN_FALSE(oE2==ker.get_size(7));

        bool ROat3rdDim = false;
        ho5DArray<T> convKerFlip;
        GADGET_CHECK_RETURN_FALSE(this->kspaceDomainConvKernel3D(ker, kRO, kE1,  kE2, oRO, oE1, oE2, convKerFlip, minusI, ROat3rdDim));

        // allocate image domain kernel
        size_t kConvE1 = convKerFlip.get_size(1);
        size_t kConvE2 = convKerFlip.get_size(2);

        kImRO.create(kConvE1, kConvE2, ro, srcCHA, dstCHA);

        hoNDArray<T> kImROTemp(ro, kConvE1, kConvE2, srcCHA, dstCHA);
        Gadgetron::clear(kImROTemp);

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - SNR unit scaling ... "); }
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(ro)) ), convKerFlip ));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(convKerFlip, debugFolder_+"convKerFlip_scal_RO"); }

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - zero padding only for RO ... "); }
        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3DNoPresetZeros(convKerFlip, ro, kConvE1, kConvE2, kImROTemp));
        // GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(convKerFlip, ro, kConvE1, kConvE2, kImROTemp, false));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(ro, kConvE1, kConvE2, &convKerFlip, &kImROTemp, false));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kImROTemp, debugFolder_+"convKerFlip_scal_RO_zeropadded"); }

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - conver to image domain only for RO ... "); }
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(kImROTemp);
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - permute kernel dimensions to be [kE1 kE2 RO ...]  ... "); }

        std::vector<size_t> dim_order(3);
        dim_order[0] = 1;
        dim_order[1] = 2;
        dim_order[2] = 0;

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&kImROTemp, &kImRO, &dim_order));

        if ( performTiming_ ) { gt_timer3_.stop(); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRIT<T>::imageDomainKernelRO3D(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRIT<T>::
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

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - SNR unit scaling for E1 and E2 ... "); }
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal( (typename realType<T>::Type)( std::sqrt((double)(e1*e2)) ), kImROScaled ));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kImROScaled, debugFolder_+"kImROScaledE1E2"); }

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - zero padding for E1 and E2 ... "); }
        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().zeropad3DNoPresetZeros(kImROScaled, e1, e2, dimR[2], kImE1E2RO));
        // GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::zeropad3D(kImROScaled, e1, e2, dimR[2], kImE1E2RO, false));
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(e1, e2, dimR[2], &kImROScaled, &kImE1E2RO, false));
        if ( performTiming_ ) { gt_timer3_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kImE1E2RO, debugFolder_+"kImE1E2RO_zeropadded_E1E2"); }

        if ( performTiming_ ) { gt_timer3_.start("spirit 3D calibration - conver to image domain for E1 and E2 ... "); }
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kImE1E2RO);
        if ( performTiming_ ) { gt_timer3_.stop(); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRIT<T>::imageDomainKernelE1E2RO(...) ... ");
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
        std::vector<size_t> dim, dimAdj, dimOrder;
        kIm.get_dimensions(dim);

        size_t NDim = dim.size();

        dimAdj = dim;
        dimAdj[NDim - 1] = dim[NDim - 2];
        dimAdj[NDim - 2] = dim[NDim - 1];

        if (!adjkIm.dimensions_equal(&dimAdj))
        {
            adjkIm.create(dimAdj);
        }

        dimOrder.resize(NDim);

        for (size_t d = 0; d < NDim; d++)
        {
            dimOrder[d] = d;
        }
        dimOrder[NDim - 2] = NDim - 1;
        dimOrder[NDim - 1] = NDim - 2;

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(const_cast< hoNDArray<T>* >(&kIm), &adjkIm, &dimOrder));

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::conjugate(adjkIm, adjkIm));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRIT<T>::imageDomainAdjointKernel(...) ... ");
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
        #pragma omp parallel default(none) private(d) shared(N, dstCHA, srcCHA, kIm, kImS2D, kImD2S) num_threads( (int)dstCHA ) if (dstCHA > 4)
        {
            hoNDArray<T> ker(N);

            std::vector<size_t> dim(1);
            dim[0] = N;

            hoNDArray<T> dKer, kerS2D, kerD2S;

            #pragma omp for
            for ( d=0; d<dstCHA; d++ )
            {
                for ( long long dprime=0; dprime<dstCHA; dprime++ )
                {
                    dKer.create(&dim, kIm.begin()+d*N+dprime*N*dstCHA);

                    for ( long long s=0; s<srcCHA; s++ )
                    {
                        kerS2D.create(&dim, const_cast<T*>(kImS2D.begin())+s*N+dprime*N*srcCHA);
                        kerD2S.create(&dim, const_cast<T*>(kImD2S.begin())+d*N+s*N*dstCHA);

                        Gadgetron::multiply(kerS2D, kerD2S, ker);
                        Gadgetron::add(dKer, ker, dKer);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRIT<T>::AdjointForwardKernel(...) ... ");
        return false;
    }

    return true;
}

}}
