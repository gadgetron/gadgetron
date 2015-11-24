/** \file       gtPlusWaveletOperator.h
    \brief      Implement wavelet operator for L1 regularization
    \author     Hui Xue
*/

#pragma once

#include "gtPlusOperator.h"
#include "hoNDHarrWavelet.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusWaveletOperator : public gtPlusOperator<T>
{
public:

    typedef gtPlusOperator<T> BaseClass;
    typedef typename BaseClass::value_type value_type;

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
    virtual bool shrinkWavCoeff(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, value_type thres, const hoNDArray<T>& mask, bool processApproxCoeff=false);
    virtual bool proximity(hoNDArray<T>& wavCoeff, value_type thres);

    // if the sensitivity S is set, compute gradient of ||wav*S'*F'*(Dc'x+D'y)||1
    // if not, compute gradient of ||wav*F'*(Dc'x+D'y)||1
    // x represents the unacquired kspace points [RO E1 CHA]
    virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g);

    // if the sensitivity S is set, compute cost value of L2 norm ||wav*S'*F'*(Dc'x+D'y)||1
    // if not, compute cost value of L2 norm ||wav*F'*(Dc'x+D'y)||1
    virtual bool obj(const hoNDArray<T>& x, T& obj);

    // number of transformation levels
    size_t numOfWavLevels_;

    // whether to include low frequency approximation coefficients
    bool with_approx_coeff_;

    T scale_factor_first_dimension_;
    T scale_factor_second_dimension_;

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;

public:

    // convert to image domain or back to kspace
    virtual bool convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual bool convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x);

    // compute gradient on the assembled kspace
    virtual bool gradTask(const hoNDArray<T>& x, hoNDArray<T>& g);

    // compute the obj on the assembled kspace
    virtual bool objTask(const hoNDArray<T>& x, T& obj);

    //// utility function
    //void scal(size_t N, float a, float* x);
    //void scal(size_t N, double a, double* x);
    //void scal(size_t N, std::complex<float> a, std::complex<float>* x);
    //void scal(size_t N, std::complex<double> a, std::complex<double>* x);

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
    hoNDArray<T> wav_coeff_norm_approx_;

    hoNDArray<value_type> wav_coeff_norm_mag_;

    hoNDArray<T> kspace_wav_;
    hoNDArray<T> complexIm_wav_;
    hoNDArray<T> complexIm_norm_;

    using BaseClass::kspace_Managed_;
    using BaseClass::complexIm_Managed_;
    using BaseClass::res_after_apply_kernel_Managed_;
    using BaseClass::res_after_apply_kernel_sum_over_Managed_;
};

template <typename T> 
gtPlusWaveletOperator<T>::gtPlusWaveletOperator() : numOfWavLevels_(1), with_approx_coeff_(false), scale_factor_first_dimension_(1.0), scale_factor_second_dimension_(1.0), BaseClass()
{

}

template <typename T> 
gtPlusWaveletOperator<T>::~gtPlusWaveletOperator()
{
}

//template <typename T> 
//void gtPlusWaveletOperator<T>::scal(size_t N, float a, float* x)
//{
//    long long n;
//    #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>64*1024)
//    for (n = 0; n < (long long)N; n++)
//    {
//        x[n] *= a;
//    }
//}
//
//template <typename T> 
//void gtPlusWaveletOperator<T>::scal(size_t N, double a, double* x)
//{
//    long long n;
//    #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>64*1024)
//    for (n = 0; n < (long long)N; n++)
//    {
//        x[n] *= a;
//    }
//}
//
//template <typename T> 
//void gtPlusWaveletOperator<T>::scal(size_t N,  std::complex<float>  a,  std::complex<float> * x)
//{
//    long long n;
//
//    #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>64*1024)
//    for (n = 0; n < (long long)N; n++)
//    {
//        const  std::complex<float> & c = x[n];
//        const float re = c.real();
//        const float im = c.imag();
//
//        const float ar = a.real();
//        const float ai = a.imag();
//
//        reinterpret_cast<float(&)[2]>(x[n])[0] = re*ar-im*ai;
//        reinterpret_cast<float(&)[2]>(x[n])[1] = re*ai+im*ar;
//    }
//}
//
//template <typename T> 
//void gtPlusWaveletOperator<T>::scal(size_t N,  std::complex<double>  a,  std::complex<double> * x)
//{
//    long long n;
//
//    #pragma omp parallel for default(none) private(n) shared(N, x, a) if (N>64*1024)
//    for (n = 0; n < (long long)N; n++)
//    {
//        const  std::complex<double> & c = x[n];
//        const double re = c.real();
//        const double im = c.imag();
//
//        const double ar = a.real();
//        const double ai = a.imag();
//
//        reinterpret_cast<double(&)[2]>(x[n])[0] = re*ar-im*ai;
//        reinterpret_cast<double(&)[2]>(x[n])[1] = re*ai+im*ar;
//    }
//}

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
        Gadgetron::multiplyConj(wavCoeff, wavCoeff, complexIm_norm_);
        // sum over CHA
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(complexIm_norm_, wavCoeffNorm, 3));
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWaveletOperator<T>::L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm) ... ");
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

        Gadgetron::sqrt(wavCoeffNorm, wav_coeff_norm_approx_);

        L1CoeffNorm = Gadgetron::asum(&wav_coeff_norm_approx_);
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWaveletOperator<T>::L1NormTotal(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm, T& L1CoeffNorm) ... ");
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

        if ( !wav_coeff_norm_approx_.dimensions_equal( &wavCoeffNorm ) )
        {
            wav_coeff_norm_approx_.create( wavCoeffNorm.get_dimensions() );
        }

        long long ii;
        long long N = (long long)wavCoeffNorm.get_number_of_elements();

        const T* pCoeffNorm = wavCoeffNorm.begin();
        T* pBuf = wav_coeff_norm_approx_.begin();

        if ( std::abs(std::abs(p) - 1.0) < 0.001 )
        {
            #pragma omp parallel for default(none) private(ii) shared(N, pBuf, pCoeffNorm, mu)
            for ( ii=0; ii<N; ii++ )
            {
                pBuf[ii] = (value_type)(1.0 / std::sqrt( pCoeffNorm[ii].real() + mu.real() ));
            }
        }
        else
        {
            #pragma omp parallel for default(none) private(ii) shared(N, pBuf, pCoeffNorm, mu, p)
            for ( ii=0; ii<N; ii++ )
            {
                pBuf[ii] = (value_type)std::pow( (double)(pCoeffNorm[ii].real() + mu.real()), (double)(p.real()/2.0-1.0) );
            }
        }

        if ( processApproxCoeff )
        {
            size_t num = wavCoeff.get_number_of_elements() / (RO*E1*W*CHA);

#pragma omp parallel default(none) private(ii) shared(RO, E1, num, wavCoeffNorm, wavCoeff, W, CHA) if ( num > 1 )
            {

#pragma omp for
                for (ii = 0; ii<(long long)num; ii++)
                {
                    hoNDArray<T> wavCoeffNormCurr(RO, E1, W, wav_coeff_norm_approx_.begin() + ii*RO*E1*W);

                    for (size_t cha = 0; cha<CHA; cha++)
                    {
                        hoNDArray<T> wavCoeffCurr(RO, E1, W, wavCoeff.begin() + ii*RO*E1*W*CHA + cha*RO*E1*W);
                        Gadgetron::multiply(wavCoeffNormCurr, wavCoeffCurr, wavCoeffCurr);
                    }
                }
            }
        }
        else
        {
            size_t num = wavCoeff.get_number_of_elements()/(RO*E1*W*CHA);

            #pragma omp parallel default(none) private(ii) shared(RO, E1, num, wavCoeffNorm, wavCoeff, W, CHA) if ( num > 1 )
            {

                #pragma omp for
                for ( ii=0; ii<(long long)num; ii++ )
                {
                    hoNDArray<T> wavCoeffNormCurr(RO, E1, W-1, wav_coeff_norm_approx_.begin()+ii*RO*E1*W+RO*E1);

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
        GERROR_STREAM("Errors in gtPlusWaveletOperator<T>::divideWavCoeffByNorm(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T mu, T p, bool processApproxCoeff) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
proximity(hoNDArray<T>& wavCoeff, value_type thres)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(this->L1Norm(wavCoeff, wav_coeff_norm_));
        hoNDArray<T> mask;

        GADGET_CHECK_RETURN_FALSE(this->shrinkWavCoeff(wavCoeff, wav_coeff_norm_, thres, mask, with_approx_coeff_));
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWaveletOperator<T>::proximity(hoNDArray<T>& wavCoeff, value_type thres) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWaveletOperator<T>::
shrinkWavCoeff(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, value_type thres, const hoNDArray<T>& mask, bool processApproxCoeff)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t E1 = wavCoeff.get_size(1);
        size_t W = wavCoeff.get_size(2);
        size_t CHA = wavCoeff.get_size(3);

        if ( !wav_coeff_norm_approx_.dimensions_equal(&wavCoeffNorm) )
        {
            wav_coeff_norm_approx_.create(wavCoeffNorm.get_dimensions());
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
        T* pMag = wav_coeff_norm_approx_.begin();
        T* pMagInv = res_after_apply_kernel_.begin();

        #pragma omp parallel for default(none) private(ii) shared(N, pMag, pMagInv, pCoeffNorm)
        for ( ii=0; ii<N; ii++ )
        {
            pMag[ii] = (value_type)std::sqrt( pCoeffNorm[ii].real() );
            pMagInv[ii] = (value_type)(1.0/(pMag[ii].real()+DBL_EPSILON));
        }

        // phase does not change
#pragma omp parallel default(none) private(ii) shared(RO, E1, num, wavCoeff, W, CHA) if ( num > 1 )
        {
#pragma omp for
            for (ii = 0; ii<num; ii++)
            {
                hoNDArray<T> MagInvCurr(RO, E1, W, res_after_apply_kernel_.begin() + ii*RO*E1*W);

                for (size_t cha = 0; cha<CHA; cha++)
                {
                    hoNDArray<T> wavCoeffCurr(RO, E1, W, wavCoeff.begin() + ii*RO*E1*W*CHA + cha*RO*E1*W);
                    hoNDArray<T> resCurr(RO, E1, W, complexIm_.begin() + ii*RO*E1*W*CHA + cha*RO*E1*W);

                    Gadgetron::multiply(MagInvCurr, wavCoeffCurr, resCurr);
                }
            }
        }

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
                    // if ( std::abs(pMagCurr[nn]) < std::abs(thres*pMaskCurr[nn]) )
                    if ( pMagCurr[nn].real() < thres*pMaskCurr[nn].real() )
                    {
                        pMagCurr[nn] = 0;
                    }
                    else
                    {
                        pMagCurr[nn] -= thres*pMaskCurr[nn];
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
                    // if ( std::abs(pMagCurr[nn]) < std::abs(thres) )
                    if ( pMagCurr[nn].real() < thres )
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

        size_t W_Start = W;
        if ( processApproxCoeff )
        {
            num = wavCoeff.get_number_of_elements() / (RO*E1*W*CHA);

#pragma omp parallel default(none) private(ii) shared(RO, E1, num, wavCoeff, W, CHA) if ( num > 1 )
            {
#pragma omp for
                for (ii = 0; ii<num; ii++)
                {
                    hoNDArray<T> MagCurr(RO, E1, W, wav_coeff_norm_approx_.begin() + ii*RO*E1*W);

                    for (size_t cha = 0; cha<CHA; cha++)
                    {
                        hoNDArray<T> phaseCurr(RO, E1, W, complexIm_.begin() + ii*RO*E1*W*CHA + cha*RO*E1*W);
                        hoNDArray<T> wavCoeffCurr(RO, E1, W, wavCoeff.begin() + ii*RO*E1*W*CHA + cha*RO*E1*W);

                        Gadgetron::multiply(MagCurr, phaseCurr, wavCoeffCurr);
                    }
                }
            }
        }
        else
        {
            num = wavCoeff.get_number_of_elements()/(RO*E1*W*CHA);

            #pragma omp parallel default(none) private(ii) shared(RO, E1, num, wavCoeff, W, CHA) if ( num > 1 )
            {
                #pragma omp for
                for ( ii=0; ii<num; ii++ )
                {
                    hoNDArray<T> MagCurr(RO, E1, W-1, wav_coeff_norm_approx_.begin()+ii*RO*E1*W+RO*E1);

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
        GERROR_STREAM("Errors in gtPlusWaveletOperator<T>::shrinkWavCoeff(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T thres, const hoNDArray<T>& mask, bool processApproxCoeff) ... ");
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
        //gt_timer1_.start("1");
        Gadgetron::multiply(unacquired_points_indicator_, x, kspace_);
        //gt_timer1_.stop();

        //gt_timer1_.start("2");
        Gadgetron::add(*acquired_points_, kspace_, kspace_);
        //gt_timer1_.stop();

        // compute the gradient on assembled kspace
        GADGET_CHECK_RETURN_FALSE(this->gradTask(kspace_, g));

        // only unacquired points are kept
        //gt_timer1_.start("12");
        Gadgetron::multiply(unacquired_points_indicator_, g, g);
        //gt_timer1_.stop();
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWaveletOperator<T>::grad(const hoNDArray<T>& x, hoNDArray<T>& g) ... ");
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
        //gt_timer1_.start("1");
        Gadgetron::multiply(unacquired_points_indicator_, x, kspace_);
        //gt_timer1_.stop();

        //gt_timer1_.start("2");
        Gadgetron::add(*acquired_points_, kspace_, kspace_);
        //gt_timer1_.stop();

        // compute the objective function on assembled kspace
        GADGET_CHECK_RETURN_FALSE(this->objTask(kspace_, obj));
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWaveletOperator<T>::obj(const hoNDArray<T>& x, T& obj) ... ");
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
        //gt_timer2_.start("3");
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, complexIm_));
        //gt_timer2_.stop();

        // compute the gradient
        if ( coil_senMap_ && coil_senMap_->get_size(0)==RO && coil_senMap_->get_size(1)==E1 && coil_senMap_->get_size(2)==CHA )
        {
            // perform coil combination
            //gt_timer2_.start("4");
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.coilCombine(complexIm_, *coil_senMap_, res_after_apply_kernel_));
            //gt_timer2_.stop();

            //gt_timer2_.start("5");
            hoNDArray<T> combined(RO, E1, 1, res_after_apply_kernel_.begin());
            //gt_timer2_.stop();

            // compute wavelet transform
            //gt_timer2_.start("6");
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(combined, res_after_apply_kernel_sum_over_));
            //gt_timer2_.stop();
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(complexIm_, res_after_apply_kernel_sum_over_));
        }

        // modify coefficients
        //gt_timer2_.start("7");
        GADGET_CHECK_RETURN_FALSE(this->L1Norm(res_after_apply_kernel_sum_over_, wav_coeff_norm_));
        //gt_timer2_.stop();

        //gt_timer2_.start("8");
        GADGET_CHECK_RETURN_FALSE(this->divideWavCoeffByNorm(res_after_apply_kernel_sum_over_, wav_coeff_norm_, T( (value_type)1e-15), T( (value_type)1.0 ), with_approx_coeff_));
        //gt_timer2_.stop();

        // go back to image
        //gt_timer2_.start("9");
        GADGET_CHECK_RETURN_FALSE(this->adjointOperator(res_after_apply_kernel_sum_over_, complexIm_wav_));
        //gt_timer2_.stop();

        if ( coil_senMap_ && coil_senMap_->get_size(0)==RO && coil_senMap_->get_size(1)==E1 && coil_senMap_->get_size(2)==CHA )
        {
            // apply coil sensivity
            //gt_timer2_.start("10");
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(*coil_senMap_, complexIm_wav_, kspace_wav_));
            //gt_timer2_.stop();

            // go to kspace
            //gt_timer2_.start("11");
            GADGET_CHECK_RETURN_FALSE(this->convertToKSpace(kspace_wav_, g));
            //gt_timer2_.stop();
        }
        else
        {
            // go to kspace
            GADGET_CHECK_RETURN_FALSE(this->convertToKSpace(complexIm_wav_, g));
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWaveletOperator<T>::gradTask(const hoNDArray<T>& x, hoNDArray<T>& g) ... ");
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
        //gt_timer3_.start("3");
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, complexIm_));
        //gt_timer3_.stop();

        // apply sensitivity
        if ( coil_senMap_ && coil_senMap_->get_size(0)==RO && coil_senMap_->get_size(1)==E1 && coil_senMap_->get_size(2)==CHA )
        {
            // perform coil combination
            //gt_timer3_.start("4");
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.coilCombine(complexIm_, *coil_senMap_, res_after_apply_kernel_));
            //gt_timer3_.stop();

            //gt_timer3_.start("5");
            hoNDArray<T> combined(RO, E1, 1, res_after_apply_kernel_.begin());
            //gt_timer3_.stop();

            // compute wavelet transform
            //gt_timer3_.start("6");
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(combined, res_after_apply_kernel_sum_over_));
            //gt_timer3_.stop();
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(complexIm_, res_after_apply_kernel_sum_over_));
        }

        //gt_timer3_.start("7");
        GADGET_CHECK_RETURN_FALSE(this->L1NormTotal(res_after_apply_kernel_sum_over_, wav_coeff_norm_, obj));
        //gt_timer3_.stop();
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWaveletOperator<T>::objTask(const hoNDArray<T>& x, T& obj) ... ");
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

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(x, im, complexIm_Managed_);

    return true;
}

template <typename T> 
inline bool gtPlusWaveletOperator<T>::convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if ( !kspace_Managed_.dimensions_equal(&im) )
    {
        kspace_Managed_.create(im.get_dimensions());
    }

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(im, x, kspace_Managed_);

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
