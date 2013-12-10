/** \file       gtPlusWaveletOperator.h
    \brief      Implement wavelet operator for L1 regularization
    \author     Hui Xue
*/

#pragma once

#include "gtPlusOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusWaveletOperator : public gtPlusOperator<T>
{
public:

    typedef gtPlusOperator<T> BaseClass;

    gtPlusWaveletOperator();
    virtual ~gtPlusWaveletOperator();

    virtual void printInfo(std::ostream& os);

    // forward operator
    virtual bool forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y) = 0;

    // adjoint operator
    virtual bool adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y) = 0;

    // compute L1 norm of wavelet coefficients across CHA
    // waveCoeff: [RO E1 W CHA ...], W is the wavelet coefficient dimension (e.g. for 1 level wavelet decomposition, W=4 for 2D and W=8 for 3D)
    // the W=1 wavelet coefficient is the most low frequent coefficients
    virtual bool L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm);
    virtual bool L1NormTotal(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm, T& L1CoeffNorm);

    // to compute the gradient of wavelet term, divide the wavelet coefficients by its norm
    // if processApproxCoeff = true, the most low frequent coefficients are changed; otherwise, remains unchanged
    virtual bool divideWavCoeffByNorm(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T mu, T p, bool processApproxCoeff=false);

    // soft-threshold or shrink the wavelet coefficients
    // the really applied threshold is mask.*thres
    virtual bool shrinkWavCoeff(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T thres, const hoNDArray<T>& mask, bool processApproxCoeff=false);

    // if the sensitivity S is set, compute gradient of ||wav*F'*S'*(Dc'x+D'y)||1
    // if not, compute gradient of ||wav*F'*(Dc'x+D'y)||1
    // x represents the unacquired kspace points [RO E1 CHA]
    virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g);

    // if the sensitivity S is set, compute cost value of L2 norm ||wav*F'*S'*(Dc'x+D'y)||1
    // if not, compute cost value of L2 norm ||wav*F'*(Dc'x+D'y)||1
    virtual bool obj(const hoNDArray<T>& x, T& obj);

    // number of transformation levels
    size_t numOfWavLevels_;

    // whether to include low frequency approximation coefficients
    bool with_approx_coeff_;

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;
    using BaseClass::gtPlus_mem_manager_;

protected:

    // convert to image domain or back to kspace
    virtual bool convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual bool convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x);

    // compute gradient on the assembled kspace
    virtual bool gradTask(const hoNDArray<T>& x, hoNDArray<T>& g);

    // compute the obj on the assembled kspace
    virtual bool objTask(const hoNDArray<T>& x, T& obj);

    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;
    using BaseClass::coil_senMap_;

    // helper memory
    using BaseClass::kspace_;
    using BaseClass::complexIm_;
    using BaseClass::res_after_apply_kernel_;
    using BaseClass::res_after_apply_kernel_sum_over_;

    hoNDArray<T> wav_coeff_norm_;
    hoNDArray<T> kspace_wav_;
    hoNDArray<T> complexIm_wav_;
    hoNDArray<T> complexIm_norm_;

    using BaseClass::kspace_Managed_;
    using BaseClass::complexIm_Managed_;
    using BaseClass::res_after_apply_kernel_Managed_;
    using BaseClass::res_after_apply_kernel_sum_over_Managed_;
};

template <typename T> 
gtPlusWaveletOperator<T>::gtPlusWaveletOperator() : numOfWavLevels_(1), with_approx_coeff_(false), BaseClass()
{

}

template <typename T> 
gtPlusWaveletOperator<T>::~gtPlusWaveletOperator()
{
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dims = wavCoeff.get_dimensions();

        std::vector<size_t> dimR(*dims);
        dimR[3] = 1;

        if ( !wavCoeffNorm.dimensions_equal(&dimR) )
        {
            wavCoeffNorm.create(&dimR);
        }

        size_t RO = (*dims)[0];
        size_t E1 = (*dims)[1];
        size_t W = (*dims)[2];
        size_t CHA = (*dims)[3];

        // square the coefficients
        GADGET_CHECK_RETURN_FALSE(Gadgetron::multiplyConj(wavCoeff, wavCoeff, complexIm_norm_));
        // sum over CHA
        GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOver4thDimension(complexIm_norm_, wavCoeffNorm));
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletOperator<T>::L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
L1NormTotal(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm, T& L1CoeffNorm)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(this->L1Norm(wavCoeff, wavCoeffNorm));

        GADGET_CHECK_RETURN_FALSE(Gadgetron::sqrt(wavCoeffNorm, kspace_));

        L1CoeffNorm = Gadgetron::asum(&kspace_);
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletOperator<T>::L1NormTotal(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm, T& L1CoeffNorm) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
divideWavCoeffByNorm(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T mu, T p, bool processApproxCoeff)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t E1 = wavCoeff.get_size(1);
        size_t W = wavCoeff.get_size(2);
        size_t CHA = wavCoeff.get_size(3);

        if ( !kspace_.dimensions_equal( &wavCoeffNorm ) )
        {
            kspace_.create( wavCoeffNorm.get_dimensions() );
        }

        long long ii;
        long long N = (long long)wavCoeffNorm.get_number_of_elements();

        const T* pCoeffNorm = wavCoeffNorm.begin();
        T* pBuf = kspace_.begin();

        if ( GT_ABS(std::abs(p) - 1.0) < 0.001 )
        {
            #pragma omp parallel for default(none) private(ii) shared(N, pBuf, pCoeffNorm, mu)
            for ( ii=0; ii<N; ii++ )
            {
                pBuf[ii] = 1.0 / std::sqrt( pCoeffNorm[ii].real() + mu.real() );
            }
        }
        else
        {
            #pragma omp parallel for default(none) private(ii) shared(N, pBuf, pCoeffNorm, mu, p)
            for ( ii=0; ii<N; ii++ )
            {
                pBuf[ii] = std::pow( (double)(pCoeffNorm[ii].real() + mu.real()), (double)(p.real()/2.0-1.0) );
            }
        }

        if ( processApproxCoeff )
        {
            GADGET_CHECK_RETURN_FALSE(Gadgetron::multiplyOver4thDimension(kspace_, wavCoeff, wavCoeff));
        }
        else
        {
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::multiplyOver4thDimensionExcept(kspace_, wavCoeff, 0, wavCoeff, true));
            size_t num = wavCoeff.get_number_of_elements()/(RO*E1*W*CHA);

            #ifdef GCC_OLD_FLAG
                #pragma omp parallel default(none) private(ii) shared(RO, E1, num, W, CHA) if ( num > 1 )
            #else
                #pragma omp parallel default(none) private(ii) shared(RO, E1, num, wavCoeffNorm, wavCoeff, W, CHA) if ( num > 1 )
            #endif
            {

                #pragma omp for
                for ( ii=0; ii<num; ii++ )
                {
                    hoNDArray<T> wavCoeffNormCurr(RO, E1, W-1, kspace_.begin()+ii*RO*E1*W+RO*E1);

                    for ( size_t cha=0; cha<CHA; cha++ )
                    {
                        hoNDArray<T> wavCoeffCurr(RO, E1, W-1, wavCoeff.begin()+ii*RO*E1*W*CHA+cha*RO*E1*W+RO*E1);
                        Gadgetron::multiply(wavCoeffNormCurr, wavCoeffCurr, wavCoeffCurr);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletOperator<T>::divideWavCoeffByNorm(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T mu, T p, bool processApproxCoeff) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
shrinkWavCoeff(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T thres, const hoNDArray<T>& mask, bool processApproxCoeff)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t E1 = wavCoeff.get_size(1);
        size_t W = wavCoeff.get_size(2);
        size_t CHA = wavCoeff.get_size(3);

        if ( !kspace_.dimensions_equal(&wavCoeffNorm) )
        {
            kspace_.create(wavCoeffNorm.get_dimensions());
        }

        if ( !res_after_apply_kernel_.dimensions_equal(&wavCoeffNorm) )
        {
            res_after_apply_kernel_.create(wavCoeffNorm.get_dimensions());
        }

        long long ii;
        long long N = (long long)wavCoeffNorm.get_number_of_elements();
        long long N3D = RO*E1*W;

        long long num = N/N3D;

        const T* pCoeffNorm = wavCoeffNorm.begin();
        T* pMag = kspace_.begin();
        T* pMagInv = res_after_apply_kernel_.begin();

        #pragma omp parallel for default(none) private(ii) shared(N, pMag, pMagInv, pCoeffNorm)
        for ( ii=0; ii<N; ii++ )
        {
            pMag[ii] = std::sqrt( pCoeffNorm[ii].real() );
            pMagInv[ii] = 1.0/(pMag[ii].real()+DBL_EPSILON);
        }

        // phase does not change
        GADGET_CHECK_RETURN_FALSE(Gadgetron::multiplyOver4thDimension(res_after_apply_kernel_, wavCoeff, complexIm_));

        // shrink the magnitude
        if ( mask.dimensions_equal(&wavCoeffNorm) )
        {
            const T* pMask = mask.begin();

            long long n = 0;
            for ( n=0; n<num; n++ )
            {
                long long s=RO*E1; 
                if ( processApproxCoeff )
                {
                    s = 0;
                }

                const T* pMaskCurr = pMask + n*N3D;
                T* pMagCurr = pMag + n*N3D;

                long long nn;

                #pragma omp parallel for private(nn) shared(s, N3D, pMagCurr, pMaskCurr, thres)
                for ( nn=s; nn<N3D; nn++ )
                {
                    if ( std::abs(pMagCurr[nn]) < std::abs(thres*pMaskCurr[nn]) )
                    {
                        pMagCurr[nn] = 0;
                    }
                    else
                    {
                        pMagCurr[nn] -= thres;
                    }
                }
            }
        }
        else
        {
            long long n = 0;
            for ( n=0; n<num; n++ )
            {
                long long s=RO*E1; 
                if ( processApproxCoeff )
                {
                    s = 0;
                }

                T* pMagCurr = pMag + n*N3D;

                long long nn;
                #pragma omp parallel for private(nn) shared(s, N3D, pMagCurr, thres)
                for ( nn=s; nn<N3D; nn++ )
                {
                    if ( std::abs(pMagCurr[nn]) < std::abs(thres) )
                    {
                        pMagCurr[nn] = 0;
                    }
                    else
                    {
                        pMagCurr[nn] -= thres;
                    }
                }
            }
        }

        if ( processApproxCoeff )
        {
            GADGET_CHECK_RETURN_FALSE(Gadgetron::multiplyOver4thDimension(kspace_, complexIm_, wavCoeff));
        }
        else
        {
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::multiplyOver4thDimensionExcept(kspace_, complexIm_, 0, wavCoeff, false));
            num = wavCoeff.get_number_of_elements()/(RO*E1*W*CHA);

            #ifdef GCC_OLD_FLAG
                #pragma omp parallel default(none) private(ii) shared(RO, E1, num, W, CHA) if ( num > 1 )
            #else
                #pragma omp parallel default(none) private(ii) shared(RO, E1, num, wavCoeff, W, CHA) if ( num > 1 )
            #endif
            {
                #pragma omp for
                for ( ii=0; ii<num; ii++ )
                {
                    hoNDArray<T> MagCurr(RO, E1, W-1, kspace_.begin()+ii*RO*E1*W+RO*E1);

                    for ( size_t cha=0; cha<CHA; cha++ )
                    {
                        hoNDArray<T> phaseCurr(RO, E1, W-1, complexIm_.begin()+ii*RO*E1*W*CHA+cha*RO*E1*W+RO*E1);
                        hoNDArray<T> wavCoeffCurr(RO, E1, W-1, wavCoeff.begin()+ii*RO*E1*W*CHA+cha*RO*E1*W+RO*E1);

                        Gadgetron::multiply(MagCurr, phaseCurr, wavCoeffCurr);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletOperator<T>::shrinkWavCoeff(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T thres, const hoNDArray<T>& mask, bool processApproxCoeff) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
grad(const hoNDArray<T>& x, hoNDArray<T>& g)
{
    try
    {
        // D'y+Dc'x
        GADGET_CHECK_RETURN_FALSE(Gadgetron::multiply(unacquired_points_indicator_, x, kspace_));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::add(*acquired_points_, kspace_, kspace_));

        // compute the gradient on assembled kspace
        GADGET_CHECK_RETURN_FALSE(this->gradTask(kspace_, g));

        // only unacquired points are kept
        GADGET_CHECK_RETURN_FALSE(Gadgetron::multiply(unacquired_points_indicator_, g, g));
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletOperator<T>::grad(const hoNDArray<T>& x, hoNDArray<T>& g) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
obj(const hoNDArray<T>& x, T& obj)
{
    try
    {
        // D'y+Dc'x
        GADGET_CHECK_RETURN_FALSE(Gadgetron::multiply(unacquired_points_indicator_, x, kspace_));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::add(*acquired_points_, kspace_, kspace_));

        // compute the objective function on assembled kspace
        GADGET_CHECK_RETURN_FALSE(this->objTask(kspace_, obj));
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletOperator<T>::obj(const hoNDArray<T>& x, T& obj) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
gradTask(const hoNDArray<T>& x, hoNDArray<T>& g)
{
    try
    {
        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t CHA = x.get_size(2);

        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, complexIm_));

        // compute the gradient
        if ( coil_senMap_ && coil_senMap_->get_size(0)==RO && coil_senMap_->get_size(1)==E1 && coil_senMap_->get_size(2)==CHA )
        {
            // perform coil combination
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.coilCombine(complexIm_, *coil_senMap_, res_after_apply_kernel_));

            hoNDArray<T> combined(RO, E1, 1, res_after_apply_kernel_.begin());

            // compute wavelet transform
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(combined, res_after_apply_kernel_sum_over_));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(complexIm_, res_after_apply_kernel_sum_over_));
        }

        // modify coefficients
        GADGET_CHECK_RETURN_FALSE(this->L1Norm(res_after_apply_kernel_sum_over_, wav_coeff_norm_));
        GADGET_CHECK_RETURN_FALSE(this->divideWavCoeffByNorm(res_after_apply_kernel_sum_over_, wav_coeff_norm_, T(1e-15), T(1.0), with_approx_coeff_));

        // go back to image
        GADGET_CHECK_RETURN_FALSE(this->adjointOperator(res_after_apply_kernel_sum_over_, complexIm_wav_));

        if ( coil_senMap_ && coil_senMap_->get_size(0)==RO && coil_senMap_->get_size(1)==E1 && coil_senMap_->get_size(2)==CHA )
        {
            // apply coil sensivity
            GADGET_CHECK_RETURN_FALSE(Gadgetron::multipleMultiply(complexIm_wav_, *coil_senMap_, kspace_wav_));

            // go to kspace
            GADGET_CHECK_RETURN_FALSE(this->convertToKSpace(kspace_wav_, g));
        }
        else
        {
            // go to kspace
            GADGET_CHECK_RETURN_FALSE(this->convertToKSpace(complexIm_wav_, g));
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletOperator<T>::gradTask(const hoNDArray<T>& x, hoNDArray<T>& g) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
objTask(const hoNDArray<T>& x, T& obj)
{
    try
    {
        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t CHA = x.get_size(2);

        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, complexIm_));

        // apply sensitivity
        if ( coil_senMap_ && coil_senMap_->get_size(0)==RO && coil_senMap_->get_size(1)==E1 && coil_senMap_->get_size(2)==CHA )
        {
            // perform coil combination
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.coilCombine(complexIm_, *coil_senMap_, res_after_apply_kernel_));

            hoNDArray<T> combined(RO, E1, 1, res_after_apply_kernel_.begin());

            // compute wavelet transform
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(combined, res_after_apply_kernel_sum_over_));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(complexIm_, res_after_apply_kernel_sum_over_));
        }

        GADGET_CHECK_RETURN_FALSE(this->L1NormTotal(res_after_apply_kernel_sum_over_, wav_coeff_norm_, obj));
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusWaveletOperator<T>::objTask(const hoNDArray<T>& x, T& obj) ... ");
        return false;
    }

    return true;
}

template <typename T> 
inline bool gtPlusWaveletOperator<T>::convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    if ( !complexIm_Managed_.dimensions_equal(&x) )
    {
        complexIm_Managed_.create(x.get_dimensions());
    }

    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(x, im, complexIm_Managed_));

    return true;
}

template <typename T> 
inline bool gtPlusWaveletOperator<T>::convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if ( !kspace_Managed_.dimensions_equal(&im) )
    {
        kspace_Managed_.create(im.get_dimensions());
    }

    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(im, x, kspace_Managed_));

    return true;
}

template <typename T> 
void gtPlusWaveletOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD wavelet operator -----------------------" << endl;
    os << "Wavelet operator for gtPlus ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

}}
