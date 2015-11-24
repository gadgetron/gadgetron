/** \file       gtPlusWavelet3DOperator.h
    \brief      Implement 3D wavelet operator for L1 regularization
    \author     Hui Xue

    Redundant haar wavelet transformation is implemented here.
*/

#pragma once

#include "gtPlusWaveletOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusWavelet3DOperator : public gtPlusWaveletOperator<T>
{
public:

    typedef gtPlusWaveletOperator<T> BaseClass;
    typedef typename BaseClass::value_type value_type;

    gtPlusWavelet3DOperator();
    virtual ~gtPlusWavelet3DOperator();

    virtual void printInfo(std::ostream& os);

    // forward operator
    // x : [RO E1 CHA E2 ...]
    // y : [RO E1 E2 W CHA ...]
    virtual bool forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    // adjoint operator
    virtual bool adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    // perform the redundant haar wavelet forward transform
    // in : [RO E1 E2], out : [RO E1 E2 1+7*level]
    // bool dwtRedundantHaar(const hoNDArray<T>& in, hoNDArray<T>& out, size_t level);

    // perform the redundant haar wavelet inverse transform
    // in : [RO E1 E2 1+7*level], out : [RO E1 E2]
    // bool idwtRedundantHaar(const hoNDArray<T>& in, hoNDArray<T>& out, size_t level);

    virtual bool unitary() const { return true; }

    // compute L1 norm of wavelet coefficients across CHA
    // waveCoeff: [RO E1 E2 W CHA ...], W is the wavelet coefficient dimension (e.g. for 1 level wavelet decomposition, W=4 for 2D and W=8 for 3D)
    // the W=1 wavelet coefficient is the most low frequent coefficients
    virtual bool L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm);

    // to compute the gradient of wavelet term, divide the wavelet coefficients by its norm
    // if processApproxCoeff = true, the most low frequent coefficients are changed; otherwise, remains unchanged
    virtual bool divideWavCoeffByNorm(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T mu, T p, bool processApproxCoeff=false);

    // soft-threshold or shrink the wavelet coefficients
    // the really applied threshold is mask.*thres
    virtual bool shrinkWavCoeff(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, value_type thres, const hoNDArray<T>& mask, bool processApproxCoeff=false);
    virtual bool proximity(hoNDArray<T>& wavCoeff, value_type thres);

    // if the sensitivity S is set, compute gradient of ||wav*F'*S'*(Dc'x+D'y)||1
    // if not, compute gradient of ||wav*F'*(Dc'x+D'y)||1
    // x represents the unacquired kspace points [RO E1 CHA E2]
    // virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g);

    // if the sensitivity S is set, compute cost value of L2 norm ||wav*F'*S'*(Dc'x+D'y)||1
    // if not, compute cost value of L2 norm ||wav*F'*(Dc'x+D'y)||1
    // virtual bool obj(const hoNDArray<T>& x, T& obj);

    // scaling along RO
    bool firstDimensionScale(hoNDArray<T>& wavCoeff, T& scaleFactor);
    // scaling along E1
    bool secondDimensionScale(hoNDArray<T>& wavCoeff, T& scaleFactor);
    // scaling along E2
    bool thirdDimensionScale(hoNDArray<T>& wavCoeff, T& scaleFactor);

    // because the spatial resolution of images are often different in through-plane dimension than the other two dimensions
    // sometime it is good to take this into account, so the regularization effects are more isotropic
    // Here only simple scaling factors are used
    // More generally, a weighting matrix can be concatenated with wavelet coefficients to enhance or suppress regularization effects as needed
    // the regularization term can become ||W*wav*F'*(Dc'x+D'y)||1, W is the general weighting matrix
    // in the next version, we shall extend this class with more geneal weighting strategy
    T scale_factor_third_dimension_;

    // in some cases, the boundary high frequency coefficients of the 3rd dimension should not be changed
    bool change_coeffcients_third_dimension_boundary_;

    using BaseClass::scale_factor_first_dimension_;
    using BaseClass::scale_factor_second_dimension_;
    using BaseClass::numOfWavLevels_;
    using BaseClass::with_approx_coeff_;
    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;

public:

    // compute gradient on the assembled kspace
    virtual bool gradTask(const hoNDArray<T>& x, hoNDArray<T>& g);

    // compute the obj on the assembled kspace
    virtual bool objTask(const hoNDArray<T>& x, T& obj);

    // help memory
    hoNDArray<T> mask_;
    hoNDArray<T> forward_buf_;
    hoNDArray<T> adjoint_buf_;

    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;
    using BaseClass::coil_senMap_;

    // helper memory
    using BaseClass::kspace_;
    using BaseClass::complexIm_;
    using BaseClass::complexIm_norm_;
    using BaseClass::res_after_apply_kernel_;
    using BaseClass::res_after_apply_kernel_sum_over_;

    using BaseClass::wav_coeff_norm_;
    using BaseClass::wav_coeff_norm_mag_;
    using BaseClass::wav_coeff_norm_approx_;

    hoNDArray<value_type> wav_coeff_norm_mag_sumCHA_;

    using BaseClass::kspace_wav_;
    using BaseClass::complexIm_wav_;

    using BaseClass::kspace_Managed_;
    using BaseClass::complexIm_Managed_;
    using BaseClass::res_after_apply_kernel_Managed_;
    using BaseClass::res_after_apply_kernel_sum_over_Managed_;
};

template <typename T> 
gtPlusWavelet3DOperator<T>::gtPlusWavelet3DOperator() : 
        scale_factor_third_dimension_(1.0), 
        change_coeffcients_third_dimension_boundary_(true), 
        BaseClass()
{

}

template <typename T> 
gtPlusWavelet3DOperator<T>::~gtPlusWavelet3DOperator()
{
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dims = x.get_dimensions();
        size_t NDim = dims->size();

        size_t RO = (*dims)[0];
        size_t E1 = (*dims)[1];
        size_t CHA = (*dims)[2];
        size_t E2 = (*dims)[3];
        size_t W = 1+7*numOfWavLevels_;

        std::vector<size_t> dimR(NDim+1);
        dimR[0] = RO;
        dimR[1] = E1;
        dimR[2] = E2;
        dimR[3] = W;
        dimR[4] = CHA;

        size_t n;
        for ( n=4; n<NDim; n++ )
        {
            dimR[n+1] = (*dims)[n];
        }

        if ( !y.dimensions_equal(&dimR) )
        {
            y.create(&dimR);
        }

        size_t num = x.get_number_of_elements()/(RO*E1*E2*CHA);

        T* pX = const_cast<T*>(x.begin());
        T* pY = y.begin();

        int t;

        if ( CHA == 1 )
        {
            #pragma omp parallel for default(none) private(t) shared(num, RO, E1, E2, W, pX, pY) if ( num > 1)
            for ( t=0; t<num; t++ )
            {
                hoNDArray<T> in(RO, E1, E2, pX+t*RO*E1*E2);
                hoNDArray<T> out(RO, E1, E2, W, pY+t*RO*E1*E2*W);
                // this->dwtRedundantHaar(in, out, numOfWavLevels_);
                Gadgetron::hoNDHarrWavelet<T> wav;
                wav.transform(in, out, 3, numOfWavLevels_, true);
            }
        }
        else
        {
            // #pragma omp parallel default(none) private(t) shared(num, RO, E1, CHA, E2, W, pX, pY) if ( num > 1 )
            {
                // hoNDArray<T> inPermute(RO, E1, E2, CHA);
                forward_buf_.create(RO, E1, E2, CHA);

                std::vector<size_t> dimOrder(4);
                dimOrder[0] = 0;
                dimOrder[1] = 1;
                dimOrder[2] = 3;
                dimOrder[3] = 2;

                // #pragma omp for
                for ( t=0; t<num; t++ )
                {
                    hoNDArray<T> in(RO, E1, CHA, E2, pX+t*RO*E1*CHA*E2);
                    Gadgetron::permute(&in, &forward_buf_, &dimOrder);

                    long long cha;

                    #pragma omp parallel for default(none) private(cha) shared(num, RO, E1, CHA, E2, W, pY, t) if ( CHA > 4 )
                    for ( cha=0; cha<(long long)CHA; cha++ )
                    {
                        hoNDArray<T> in_dwt(RO, E1, E2, forward_buf_.begin()+cha*RO*E1*E2);
                        hoNDArray<T> out(RO, E1, E2, W, pY+t*RO*E1*E2*W*CHA+cha*RO*E1*E2*W);

                        // this->dwtRedundantHaar(in_dwt, out, numOfWavLevels_);
                        Gadgetron::hoNDHarrWavelet<T> wav;
                        wav.transform(in_dwt, out, 3, numOfWavLevels_, true);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dims = x.get_dimensions();
        size_t NDim = dims->size();

        size_t RO = (*dims)[0];
        size_t E1 = (*dims)[1];
        size_t E2 = (*dims)[2];
        size_t W = (*dims)[3];
        size_t CHA = (*dims)[4];

        std::vector<size_t> dimR(NDim-1);
        dimR[0] = RO;
        dimR[1] = E1;
        dimR[2] = CHA;
        dimR[3] = E2;

        size_t n;
        for ( n=4; n<NDim-1; n++ )
        {
            dimR[n] = (*dims)[n+1];
        }

        if ( !y.dimensions_equal(&dimR) )
        {
            y.create(&dimR);
        }

        size_t num = x.get_number_of_elements()/(RO*E1*E2*W*CHA);

        T* pX = const_cast<T*>(x.begin());
        T* pY = y.begin();

        int t;

        if ( CHA == 1 )
        {
            #pragma omp parallel for default(none) private(t) shared(num, RO, E1, E2, W, pX, pY) if ( num > 1)
            for ( t=0; t<num; t++ )
            {
                hoNDArray<T> in(RO, E1, E2, W, pX+t*RO*E1*E2*W);
                hoNDArray<T> out(RO, E1, E2, pY+t*RO*E1*E2);
                // this->idwtRedundantHaar(in, out, numOfWavLevels_);
                Gadgetron::hoNDHarrWavelet<T> wav;
                wav.transform(in, out, 3, numOfWavLevels_, false);
            }
        }
        else
        {
            // #pragma omp parallel default(none) private(t) shared(num, RO, E1, CHA, E2, W, pX, pY) if ( num > 1 ) num_threads( (int)((num>16) ? 16 : num))
            {
                // hoNDArray<T> outPermute(RO, E1, E2, CHA);
                adjoint_buf_.create(RO, E1, E2, CHA);

                std::vector<size_t> dimOrder(4);
                dimOrder[0] = 0;
                dimOrder[1] = 1;
                dimOrder[2] = 3;
                dimOrder[3] = 2;

                // #pragma omp for
                for ( t=0; t<num; t++ )
                {
                    hoNDArray<T> out(RO, E1, CHA, E2, pY+t*RO*E1*CHA*E2);

                    long long cha;
                    #pragma omp parallel for default(none) private(cha) shared(RO, E1, CHA, E2, W, pX) if ( CHA > 4 )
                    for ( cha=0; cha<(long long)CHA; cha++ )
                    {
                        hoNDArray<T> in(RO, E1, E2, W, pX+cha*RO*E1*E2*W);
                        hoNDArray<T> out_idwt(RO, E1, E2, adjoint_buf_.begin()+cha*RO*E1*E2);

                        // this->idwtRedundantHaar(in, out_idwt, numOfWavLevels_);
                        Gadgetron::hoNDHarrWavelet<T> wav;
                        wav.transform(in, out_idwt, 3, numOfWavLevels_, false);
                    }

                    Gadgetron::permute(&adjoint_buf_, &out, &dimOrder);
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dims = wavCoeff.get_dimensions();

        std::vector<size_t> dimR(*dims);
        dimR[4] = 1;

        if ( !wavCoeffNorm.dimensions_equal(&dimR) )
        {
            wavCoeffNorm.create(&dimR);
        }

        size_t RO = (*dims)[0];
        size_t E1 = (*dims)[1];
        size_t E2 = (*dims)[2];
        size_t W = (*dims)[3];
        size_t CHA = (*dims)[4];

        // square the coefficients
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiplyConj(wavCoeff, wavCoeff, complexIm_norm_));
        // sum over CHA
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(complexIm_norm_, wavCoeffNorm, 4));
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
divideWavCoeffByNorm(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T mu, T p, bool processApproxCoeff)
{
    try
    {
        long long RO = (long long)wavCoeff.get_size(0);
        long long E1 = (long long)wavCoeff.get_size(1);
        long long E2 = (long long)wavCoeff.get_size(2);
        long long W = (long long)wavCoeff.get_size(3);
        long long CHA = (long long)wavCoeff.get_size(4);

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
                pBuf[ii] = (value_type)( 1.0 / std::sqrt( pCoeffNorm[ii].real() + mu.real() ) );
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
            long long num = wavCoeff.get_number_of_elements() / (RO*E1*E2*W*CHA);

            #pragma omp parallel default(none) private(ii) shared(RO, E1, E2, num, wavCoeffNorm, wavCoeff, W, CHA) if ( num > 1 )
            {

                #pragma omp for
                for (ii = 0; ii<num; ii++)
                {
                    hoNDArray<T> wavCoeffNormCurr(RO, E1, E2, W, wav_coeff_norm_approx_.begin() + ii*RO*E1*E2*W);

                    for (long long cha = 0; cha<CHA; cha++)
                    {
                        hoNDArray<T> wavCoeffCurr(RO, E1, E2, W, wavCoeff.begin() + ii*RO*E1*E2*W*CHA + cha*RO*E1*E2*W);
                        Gadgetron::multiply(wavCoeffNormCurr, wavCoeffCurr, wavCoeffCurr);
                    }
                }
            }
        }
        else
        {
            long long num = wavCoeff.get_number_of_elements()/(RO*E1*E2*W*CHA);

            #pragma omp parallel default(none) private(ii) shared(RO, E1, E2, num, wavCoeffNorm, wavCoeff, W, CHA) if ( num > 1 )
            {

                #pragma omp for
                for ( ii=0; ii<num; ii++ )
                {
                    hoNDArray<T> wavCoeffNormCurr(RO, E1, E2, W-1, wav_coeff_norm_approx_.begin()+ii*RO*E1*E2*W+RO*E1*E2);

                    for ( long long cha=0; cha<CHA; cha++ )
                    {
                        hoNDArray<T> wavCoeffCurr(RO, E1, E2, W-1, wavCoeff.begin()+ii*RO*E1*E2*W*CHA+cha*RO*E1*E2*W+RO*E1*E2);
                        Gadgetron::multiply(wavCoeffNormCurr, wavCoeffCurr, wavCoeffCurr);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::divideWavCoeffByNorm(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T mu, T p, bool processApproxCoeff) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
proximity(hoNDArray<T>& wavCoeff, value_type thres)
{
    try
    {
        // GADGET_CHECK_RETURN_FALSE(this->L1Norm(wavCoeff, wav_coeff_norm_));

        // GADGET_CHECK_RETURN_FALSE(Gadgetron::multiplyConj(wavCoeff, wavCoeff, wav_coeff_norm_));
        Gadgetron::abs(wavCoeff, wav_coeff_norm_mag_);

        if ( !mask_.dimensions_equal(&wavCoeff) )
        {
            mask_.create(wavCoeff.get_dimensions());
        }

        Gadgetron::fill(mask_, T(thres) );

        if ( std::abs(std::abs(scale_factor_first_dimension_)-1.0) > 1e-6 )
        {
            GADGET_CHECK_RETURN_FALSE(this->firstDimensionScale(mask_, scale_factor_first_dimension_));
        }

        if ( std::abs(std::abs(scale_factor_second_dimension_)-1.0) > 1e-6 )
        {
            GADGET_CHECK_RETURN_FALSE(this->secondDimensionScale(mask_, scale_factor_second_dimension_));
        }

        if ( std::abs(std::abs(scale_factor_third_dimension_)-1.0) > 1e-6 )
        {
            GADGET_CHECK_RETURN_FALSE(this->thirdDimensionScale(mask_, scale_factor_third_dimension_));
        }

        GADGET_CHECK_RETURN_FALSE(this->shrinkWavCoeff(wavCoeff, wav_coeff_norm_mag_, thres, mask_, this->with_approx_coeff_));
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::proximity(hoNDArray<T>& wavCoeff, T thres) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
shrinkWavCoeff(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, value_type thres, const hoNDArray<T>& mask, bool processApproxCoeff)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dims = wavCoeff.get_dimensions();

        long long RO = (long long)(*dims)[0];
        long long E1 = (long long)(*dims)[1];
        long long E2 = (long long)(*dims)[2];
        long long W = (long long)(*dims)[3];
        long long CHA = (long long)(*dims)[4];

        if ( !wav_coeff_norm_approx_.dimensions_equal(&wavCoeffNorm) )
        {
            wav_coeff_norm_approx_.create(wavCoeffNorm.get_dimensions());
        }

        long long ii;
        long long N = (long long)wavCoeffNorm.get_number_of_elements();
        long long N4D = RO*E1*E2*W;

        long long num = N/N4D;

        value_type* pCoeffNorm = const_cast<value_type*>(wavCoeffNorm.begin());
        T* pMag = wav_coeff_norm_approx_.begin();

        if ( wavCoeffNorm.dimensions_equal(&wavCoeff) )
        {
            #pragma omp parallel for default(none) private(ii) shared(N, pMag, pCoeffNorm)
            for ( ii=0; ii<N; ii++ )
            {
                pMag[ii] = pCoeffNorm[ii];
            }

            Gadgetron::divide(wavCoeff, wav_coeff_norm_approx_, complexIm_);
        }
        else
        {
            if ( !res_after_apply_kernel_.dimensions_equal(&wavCoeffNorm) )
            {
                res_after_apply_kernel_.create(wavCoeffNorm.get_dimensions());
            }

            T* pMagInv = res_after_apply_kernel_.begin();

            #pragma omp parallel for default(none) private(ii) shared(N, pMag, pMagInv, pCoeffNorm)
            for ( ii=0; ii<N; ii++ )
            {
                pMag[ii] = pCoeffNorm[ii];
                pMagInv[ii] = 1/(pCoeffNorm[ii]+FLT_EPSILON);
            }

            // Gadgetron::inv(wav_coeff_norm_approx_, res_after_apply_kernel_);

            // phase does not change
            if ( res_after_apply_kernel_.dimensions_equal(&wavCoeff) )
            {
                Gadgetron::multiply(res_after_apply_kernel_, wavCoeff, complexIm_);
            }
            else
            {
                long long num = wavCoeff.get_number_of_elements() / (RO*E1*E2*W*CHA);

                #pragma omp parallel default(none) private(ii) shared(RO, E1, E2, num, wavCoeff, W, CHA) if ( num > 1 )
                {

                    #pragma omp for
                    for (ii = 0; ii<num; ii++)
                    {
                        hoNDArray<T> magInvCurr(RO, E1, E2, W, res_after_apply_kernel_.begin() + ii*RO*E1*E2*W);

                        for (long long cha = 0; cha<CHA; cha++)
                        {
                            hoNDArray<T> wavCoeffCurr(RO, E1, E2, W, wavCoeff.begin() + ii*RO*E1*E2*W*CHA + cha*RO*E1*E2*W);
                            hoNDArray<T> resCurr(RO, E1, E2, W, complexIm_.begin() + ii*RO*E1*E2*W*CHA + cha*RO*E1*E2*W);

                            Gadgetron::multiply(magInvCurr, wavCoeffCurr, resCurr);
                        }
                    }
                }
            }
        }

        // shrink the magnitude
        if ( mask.dimensions_equal(&wavCoeffNorm) )
        {
            const T* pMask = mask.begin();

            // value_type* pMagCHA = wav_coeff_norm_mag_sumCHA_.begin();

            long long n = 0;
            for ( n=0; n<num; n++ )
            {
                long long s=RO*E1*E2; 
                if ( processApproxCoeff )
                {
                    s = 0;
                }

                const T* pMaskCurr = pMask + n*N4D;
                T* pMagCurr = pMag + n*N4D;

                if ( change_coeffcients_third_dimension_boundary_ )
                {
                    long long nn;
                    #pragma omp parallel for private(nn) shared(s, N4D, pMagCurr, pMaskCurr, thres)
                    for ( nn=s; nn<N4D; nn++ )
                    {
                        // if ( std::abs(pMagCurr[nn]) < std::abs(thres*pMaskCurr[nn]) )
                        if ( pMagCurr[nn].real() < pMaskCurr[nn].real() )
                        // if ( pMagCHA[nn] < pMaskCurr[nn].real() )
                        {
                            pMagCurr[nn] = 0;
                        }
                        else
                        {
                            pMagCurr[nn] -= pMaskCurr[nn];
                        }
                    }
                }
                else
                {
                    // approx coefficents
                    long long nn;
                    #pragma omp parallel for private(nn) shared(s, N4D, pMagCurr, pMaskCurr, thres)
                    for ( nn=s; nn<RO*E1*E2; nn++ )
                    {
                        //if ( std::abs(pMagCurr[nn]) < std::abs(thres*pMaskCurr[nn]) )
                        if ( pMagCurr[nn].real() < pMaskCurr[nn].real() )
                        {
                            pMagCurr[nn] = 0;
                        }
                        else
                        {
                            pMagCurr[nn] -= pMaskCurr[nn];
                        }
                    }

                    size_t level;
                    for ( level=0; level<numOfWavLevels_; level++ )
                    {
                        size_t start = RO*E1*E2 + 7*level;

                        size_t w;
                        for ( w=0; w<7; w++ )
                        {
                            size_t startW = start+w*RO*E1*E2;
                            size_t endW = startW+RO*E1*E2;

                            if ( w >= 3 )
                            {
                                startW += RO*E1;
                                endW -= RO*E1;
                            }

                            long long nn;
                            #pragma omp parallel for private(nn) shared(s, N4D, pMagCurr, pMaskCurr, thres)
                            for ( nn=(long long)startW; nn<(long long)endW; nn++ )
                            {
                                // if ( std::abs(pMagCurr[nn]) < std::abs(thres*pMaskCurr[nn]) )
                                if ( pMagCurr[nn].real() < pMaskCurr[nn].real() )
                                {
                                    pMagCurr[nn] = 0;
                                }
                                else
                                {
                                    pMagCurr[nn] -= pMaskCurr[nn];
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            long long n = 0;
            for ( n=0; n<num; n++ )
            {
                long long s=RO*E1*E2; 
                if ( processApproxCoeff )
                {
                    s = 0;
                }

                T* pMagCurr = pMag + n*N4D;

                if ( change_coeffcients_third_dimension_boundary_ )
                {
                    long long nn;
                    #pragma omp parallel for private(nn) shared(s, N4D, pMagCurr, thres)
                    for ( nn=s; nn<N4D; nn++ )
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
                else
                {
                    // approx coefficents
                    long long nn;
                    #pragma omp parallel for private(nn) shared(s, N4D, pMagCurr, thres)
                    for ( nn=s; nn<RO*E1*E2; nn++ )
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

                    size_t level;
                    for ( level=0; level<numOfWavLevels_; level++ )
                    {
                        size_t start = RO*E1*E2 + 7*level;

                        size_t w;
                        for ( w=0; w<7; w++ )
                        {
                            size_t startW = start+w*RO*E1*E2;
                            size_t endW = startW+RO*E1*E2;

                            if ( w >= 3 )
                            {
                                startW += RO*E1;
                                endW -= RO*E1;
                            }

                            long long nn;
                            #pragma omp parallel for private(nn) shared(s, N4D, pMagCurr, thres)
                            for ( nn=(long long)startW; nn<(long long)endW; nn++ )
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
                }
            }
        }

        if ( processApproxCoeff )
        {
            if ( wav_coeff_norm_approx_.dimensions_equal(&complexIm_) )
            {
                Gadgetron::multiply(wav_coeff_norm_approx_, complexIm_, wavCoeff);
            }
            else
            {
                long long num = wavCoeff.get_number_of_elements() / (RO*E1*E2*W*CHA);

                #pragma omp parallel default(none) private(ii) shared(RO, E1, E2, num, wavCoeffNorm, wavCoeff, W, CHA) if ( num > 1 )
                {

                    #pragma omp for
                    for (ii = 0; ii<num; ii++)
                    {
                        hoNDArray<T> magCurr(RO, E1, E2, W, wav_coeff_norm_approx_.begin() + ii*RO*E1*E2*W);

                        for (long long cha = 0; cha<CHA; cha++)
                        {
                            hoNDArray<T> phaseCurr(RO, E1, E2, W, complexIm_.begin() + ii*RO*E1*E2*W*CHA + cha*RO*E1*E2*W);
                            hoNDArray<T> wavCoeffCurr(RO, E1, E2, W, wavCoeff.begin() + ii*RO*E1*E2*W*CHA + cha*RO*E1*E2*W);

                            Gadgetron::multiply(magCurr, phaseCurr, wavCoeffCurr);
                        }
                    }
                }
            }
        }
        else
        {
            if ( wav_coeff_norm_approx_.dimensions_equal(&wavCoeff) )
            {
                #pragma omp parallel default(none) private(ii) shared(RO, E1, E2, wavCoeffNorm, wavCoeff, W, CHA) if ( CHA > 1 )
                {

                    #pragma omp for
                    for ( ii=0; ii<CHA; ii++ )
                    {
                        hoNDArray<T> magCurr(RO, E1, E2, W-1, wav_coeff_norm_approx_.begin()+ii*RO*E1*E2*W+RO*E1*E2);
                        hoNDArray<T> phaseCurr(RO, E1, E2, W-1, complexIm_.begin()+ii*RO*E1*E2*W+RO*E1*E2);
                        hoNDArray<T> wavCoeffCurr(RO, E1, E2, W-1, wavCoeff.begin()+ii*RO*E1*E2*W+RO*E1*E2);

                        Gadgetron::multiply(magCurr, phaseCurr, wavCoeffCurr);
                    }
                }
            }
            else
            {
                long long num = wavCoeff.get_number_of_elements()/(RO*E1*E2*W*CHA);

                #pragma omp parallel default(none) private(ii) shared(RO, E1, E2, num, wavCoeffNorm, wavCoeff, W, CHA) if ( num > 1 )
                {

                    #pragma omp for
                    for ( ii=0; ii<num; ii++ )
                    {
                        hoNDArray<T> magCurr(RO, E1, E2, W-1, wav_coeff_norm_approx_.begin()+ii*RO*E1*E2*W+RO*E1*E2);

                        for ( long long cha=0; cha<CHA; cha++ )
                        {
                            hoNDArray<T> phaseCurr(RO, E1, E2, W-1, complexIm_.begin()+ii*RO*E1*E2*W*CHA+cha*RO*E1*E2*W+RO*E1*E2);
                            hoNDArray<T> wavCoeffCurr(RO, E1, E2, W-1, wavCoeff.begin()+ii*RO*E1*E2*W*CHA+cha*RO*E1*E2*W+RO*E1*E2);

                            Gadgetron::multiply(magCurr, phaseCurr, wavCoeffCurr);
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::shrinkWavCoeff(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T thres, const hoNDArray<T>& mask, bool processApproxCoeff) ... ");
        return false;
    }
    return true;
}

//template <typename T> 
//bool gtPlusWavelet3DOperator<T>::
//dwtRedundantHaar(const hoNDArray<T>& in, hoNDArray<T>& out, size_t level)
//{
//    try
//    {
//        long long RO = (long long)in.get_size(0);
//        long long E1 = (long long)in.get_size(1);
//        long long E2 = (long long)in.get_size(2);
//
//        T* pOut = out.begin();
//        memcpy(pOut, in.begin(), sizeof(T)*RO*E1*E2);
//
//        long long N2D = RO*E1;
//        long long N3D = RO*E1*E2;
//
//        for (size_t n=0; n<level; n++)
//        {
//            T* lll = pOut;
//            T* llh = lll + n*7*N3D + N3D;
//            T* lhl = llh + N3D;
//            T* lhh = lhl + N3D;
//            T* hll = lhh + N3D;
//            T* hlh = hll + N3D;
//            T* hhl = hlh + N3D;
//            T* hhh = hhl + N3D;
//
//            long long e2;
//            #pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, N2D, lll, llh)
//            for (e2=0; e2<E2; e2++)
//            {
//                long long ind3D = e2 * N2D;
//                for (long long ro=0; ro<RO; ro++)
//                {
//                    T v1 = lll[ro + ind3D];
//
//                    long long ind = ro + ind3D;
//                    for (long long e1=0; e1<E1-1; e1++)
//                    {
//                        llh[ind] = lll[ind] - lll[ind+RO];
//                        lll[ind] += lll[ind+RO];
//                        ind += RO;
//                    }
//
//                    llh[ind] = lll[ind] - v1;
//                    lll[ind] += v1;
//                }
//            }
//
//            this->scal( N3D, T(0.5), lll);
//            this->scal( N3D, T(0.5), llh );
//
//            #pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, N2D, lll, llh, lhh, lhl)
//            for (e2=0; e2<E2; e2++)
//            {
//                long long ind3D = e2*N2D;
//                for (long long e1=0; e1<E1; e1++)
//                {
//                    T v1 = lll[e1*RO + ind3D];
//                    T v2 = llh[e1*RO + ind3D];
//
//                    long long ind = e1*RO + ind3D;
//                    for (long long ro=0; ro<RO-1; ro++)
//                    {
//                        lhh[ind] = llh[ind] - llh[ind + 1];
//                        llh[ind] += llh[ind + 1];
//
//                        lhl[ind] = lll[ind] - lll[ind + 1];
//                        lll[ind] += lll[ind + 1];
//
//                        ind++;
//                    }
//
//                    lhl[ind] = lll[ind] - v1;
//                    lll[ind] += v1;
//
//                    lhh[ind] = llh[ind] - v2;
//                    llh[ind] += v2;
//                }
//            }
//
//            #pragma omp parallel sections
//            {
//                #pragma omp section
//                this->scal( N3D, T(0.5), lll );
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), lhl );
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), llh );
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), lhh );
//            }
//
//            long long e1;
//            #pragma omp parallel for default(none) private(e1) shared(RO, E1, E2, N2D, lll, hll, lhl, hhl, llh, hlh, lhh, hhh)
//            for (e1=0; e1<E1; e1++)
//            {
//                for (long long ro=0; ro<RO; ro++)
//                {
//                    long long ind2D = e1*RO + ro;
//
//                    T v1 = lll[ind2D];
//                    T v2 = lhl[ind2D];
//                    T v3 = llh[ind2D];
//                    T v4 = lhh[ind2D];
//
//                    long long ind = ind2D;
//                    for (long long e2=0; e2<E2-1; e2++)
//                    {
//                        hll[ind] = lll[ind] - lll[ind + N2D];
//                        lll[ind] += lll[ind + N2D];
//
//                        hhl[ind] = lhl[ind] - lhl[ind + N2D];
//                        lhl[ind] += lhl[ind + N2D];
//
//                        hlh[ind] = llh[ind] - llh[ind + N2D];
//                        llh[ind] += llh[ind + N2D];
//
//                        hhh[ind] = lhh[ind] - lhh[ind + N2D];
//                        lhh[ind] += lhh[ind + N2D];
//
//                        ind += N2D;
//                    }
//
//                    if ( E2 > 1 )
//                    {
//                        hll[ind] = lll[ind] - v1;
//                        lll[ind] += v1;
//
//                        hhl[ind] = lhl[ind] - v2;
//                        lhl[ind] += v2;
//
//                        hlh[ind] = llh[ind] - v3;
//                        llh[ind] += v3;
//
//                        hhh[ind] = lhh[ind] - v4;
//                        lhh[ind] += v4;
//                    }
//                }
//            }
//
//            #pragma omp parallel sections
//            {
//                #pragma omp section
//                this->scal( N3D, T(0.5), lll);
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), hll);
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), lhl);
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), hhl);
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), llh);
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), hlh);
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), lhh);
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), hhh);
//            }
//        }
//    }
//    catch (...)
//    {
//        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::dwtRedundantHaar(const hoNDArray<T>& in, hoNDArray<T>& out, size_t level) ... ");
//        return false;
//    }
//    return true;
//}
//
//template <typename T> 
//bool gtPlusWavelet3DOperator<T>::
//idwtRedundantHaar(const hoNDArray<T>& in, hoNDArray<T>& out, size_t level)
//{
//    try
//    {
//        long long RO = (long long)in.get_size(0);
//        long long E1 = (long long)in.get_size(1);
//        long long E2 = (long long)in.get_size(2);
//
//        T* pIn = const_cast<T*>(in.begin());
//        T* pOut = out.begin();
//        memcpy(pOut, in.begin(), sizeof(T)*RO*E1*E2);
//
//        long long N2D = RO*E1;
//        long long N3D = RO*E1*E2;
//
//        hoNDArray<T> LL(N3D);
//        T* pLL = LL.begin();
//
//        hoNDArray<T> HL(N3D);
//        T* pHL = HL.begin();
//
//        hoNDArray<T> LH(N3D);
//        T* pLH = LH.begin();
//
//        hoNDArray<T> HH(N3D);
//        T* pHH = HH.begin();
//
//        long long n;
//        for (n=(long long)level-1; n>=0; n--)
//        {
//            T* lll = pOut;
//            T* llh = pIn + n*7*N3D + N3D;
//            T* lhl = llh + N3D;
//            T* lhh = lhl + N3D;
//            T* hll = lhh + N3D;
//            T* hlh = hll + N3D;
//            T* hhl = hlh + N3D;
//            T* hhh = hhl + N3D;
//
//            long long e1;
//            #pragma omp parallel for default(none) private(e1) shared(RO, E1, E2, N2D, lll, hll, lhl, hhl, llh, hlh, lhh, hhh, pLL, pHL, pLH, pHH) 
//            for (e1=0; e1<E1; e1++)
//            {
//                for (long long ro=0; ro<RO; ro++)
//                {
//                    long long ind2D = e1*RO + ro;
//
//                    long long ind;
//                    for (long long e2=E2-1; e2>0; e2--)
//                    {
//                        ind = ind2D + e2*N2D;
//                        pLL[ind] = (lll[ind]+lll[ind-N2D]) + (hll[ind]-hll[ind-N2D]);
//                        pHL[ind] = (lhl[ind]+lhl[ind-N2D]) + (hhl[ind]-hhl[ind-N2D]);
//                        pLH[ind] = (llh[ind]+llh[ind-N2D]) + (hlh[ind]-hlh[ind-N2D]);
//                        pHH[ind] = (lhh[ind]+lhh[ind-N2D]) + (hhh[ind]-hhh[ind-N2D]);
//                    }
//
//                    if ( E2 > 1 )
//                    {
//                        ind = ind2D + (E2-1)*N2D;
//                        pLL[ind2D] = (lll[ind2D]+lll[ind]) + (hll[ind2D]-hll[ind]);
//                        pHL[ind2D] = (lhl[ind2D]+lhl[ind]) + (hhl[ind2D]-hhl[ind]);
//                        pLH[ind2D] = (llh[ind2D]+llh[ind]) + (hlh[ind2D]-hlh[ind]);
//                        pHH[ind2D] = (lhh[ind2D]+lhh[ind]) + (hhh[ind2D]-hhh[ind]);
//                    }
//                }
//            }
//
//            #pragma omp parallel sections
//            {
//                #pragma omp section
//                this->scal( N3D, T(0.5), pLL );
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), pHL );
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), pLH );
//
//                #pragma omp section
//                this->scal( N3D, T(0.5), pHH );
//            }
//
//            long long e2;
//            #pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, N2D, pLL, pHL, pLH, pHH) 
//            for (e2=0; e2<(long long)E2; e2++)
//            {
//                long long ind3D = e2*N2D;
//                for (long long e1=0; e1<(long long)E1; e1++)
//                {
//                    long long ind = e1*RO + RO-1 + ind3D;
//
//                    T v1 = pLL[ind];
//                    T v2 = pLH[ind];
//
//                    for (long long ro=(long long)RO-1; ro>0; ro--)
//                    {
//                        pLL[ind] = (pLL[ind]+pLL[ind-1]) + (pHL[ind]-pHL[ind-1]);
//                        pLH[ind] = (pLH[ind]+pLH[ind-1]) + (pHH[ind]-pHH[ind-1]);
//                        ind--;
//                    }
//
//                    pLL[ind] = (pLL[ind]+v1) + (pHL[ind]-pHL[ind+RO-1]);
//                    pLH[ind] = (pLH[ind]+v2) + (pHH[ind]-pHH[ind+RO-1]);
//                }
//            }
//
//            this->scal( N3D, T(0.5), pLL );
//            this->scal( N3D, T(0.5), pLH );
//
//            #pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, N2D, pLL,pLH, pOut) 
//            for (e2=0; e2<(long long)E2; e2++)
//            {
//                long long ind3D = e2*N2D;
//                for (long long ro=0; ro<(long long)RO; ro++)
//                {
//                    long long ind = (E1-1)*RO + ro + ind3D;
//                    for (long long e1=(long long)E1-1; e1>0; e1--)
//                    {
//                        pOut[ind] = (pLL[ind]+pLL[ind-RO]) + (pLH[ind]-pLH[ind-RO]);
//                        ind -= RO;
//                    }
//
//                    pOut[ind] = (pLL[ind]+pLL[ind+(E1-1)*RO]) + (pLH[ind]-pLH[ind+(E1-1)*RO]);
//                }
//            }
//
//            this->scal( N3D, T(0.5), pOut );
//        }
//    }
//    catch (...)
//    {
//        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::idwtRedundantHaar(const hoNDArray<T>& in, hoNDArray<T>& out, size_t level) ... ");
//        return false;
//    }
//    return true;
//}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
gradTask(const hoNDArray<T>& x, hoNDArray<T>& g)
{
    try
    {
        // x to image domain
        //gt_timer2_.start("3");
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, complexIm_));
        //gt_timer2_.stop();

        size_t RO = complexIm_.get_size(0);
        size_t E1 = complexIm_.get_size(1);
        size_t CHA = complexIm_.get_size(2);
        size_t E2 = complexIm_.get_size(3);

        // compute the gradient
        if ( coil_senMap_ && coil_senMap_->get_size(0)==RO && coil_senMap_->get_size(1)==E1 && coil_senMap_->get_size(2)==CHA )
        {
            // perform coil combination
            //gt_timer2_.start("4");
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.coilCombine(complexIm_, *coil_senMap_, res_after_apply_kernel_));
            //gt_timer2_.stop();

            //gt_timer2_.start("5");
            hoNDArray<T> combined(RO, E1, 1, E2, res_after_apply_kernel_.begin());
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
        GADGET_CHECK_RETURN_FALSE(this->divideWavCoeffByNorm(res_after_apply_kernel_sum_over_, wav_coeff_norm_, (value_type)(1e-15), (value_type)(1.0), with_approx_coeff_));
        //gt_timer2_.stop();

        // first dimension scaling
        //gt_timer2_.start("9");
        if ( std::abs(std::abs(scale_factor_first_dimension_)-1.0) > 1e-6 )
        {
            GADGET_CHECK_RETURN_FALSE(this->firstDimensionScale(res_after_apply_kernel_sum_over_, scale_factor_first_dimension_));
        }

        // second dimension scaling
        if ( std::abs(std::abs(scale_factor_second_dimension_)-1.0) > 1e-6 )
        {
            GADGET_CHECK_RETURN_FALSE(this->secondDimensionScale(res_after_apply_kernel_sum_over_, scale_factor_second_dimension_));
        }

        // third dimension scaling
        if ( std::abs(std::abs(scale_factor_third_dimension_)-1.0) > 1e-6 )
        {
            GADGET_CHECK_RETURN_FALSE(this->thirdDimensionScale(res_after_apply_kernel_sum_over_, scale_factor_third_dimension_));
        }
        //gt_timer2_.stop();

        // go back to image
        //gt_timer2_.start("10");
        GADGET_CHECK_RETURN_FALSE(this->adjointOperator(res_after_apply_kernel_sum_over_, complexIm_wav_));
        //gt_timer2_.stop();

        if ( coil_senMap_ && coil_senMap_->get_size(0)==RO && coil_senMap_->get_size(1)==E1 && coil_senMap_->get_size(2)==CHA )
        {
            // apply coil sensivity
            //gt_timer2_.start("11");
            if ( !kspace_wav_.dimensions_equal(&complexIm_) )
            {
                kspace_wav_.create(RO, E1, CHA, E2);
            }

            for ( size_t e2=0; e2<E2; e2++ )
            {
                hoNDArray<T> complexImE2(RO, E1, complexIm_wav_.begin()+e2*RO*E1);
                hoNDArray<T> kspace_wavE2(RO, E1, CHA, kspace_wav_.begin()+e2*RO*E1*CHA);

                if ( coil_senMap_->get_size(3) == E2 )
                {
                    hoNDArray<T> coilMapE2(RO, E1, CHA, coil_senMap_->begin()+e2*RO*E1*CHA);
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(coilMapE2, complexImE2, kspace_wavE2));
                }
                else
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(*coil_senMap_, complexImE2, kspace_wavE2));
                }
            }
            //gt_timer2_.stop();

            // go to kspace
            //gt_timer2_.start("12");
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
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::gradTask(const hoNDArray<T>& x, hoNDArray<T>& g) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
objTask(const hoNDArray<T>& x, T& obj)
{
    try
    {
        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, complexIm_));

        size_t RO = complexIm_.get_size(0);
        size_t E1 = complexIm_.get_size(1);
        size_t CHA = complexIm_.get_size(2);
        size_t E2 = complexIm_.get_size(3);

        // apply sensitivity
        if (  coil_senMap_ && coil_senMap_->get_size(0)==RO && coil_senMap_->get_size(1)==E1 && coil_senMap_->get_size(2)==CHA )
        {
            // perform coil combination
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.coilCombine(complexIm_, *coil_senMap_, res_after_apply_kernel_));

            hoNDArray<T> combined(RO, E1, 1, E2, res_after_apply_kernel_.begin());

            // compute wavelet transform
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(combined, res_after_apply_kernel_sum_over_));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(complexIm_, res_after_apply_kernel_sum_over_));
        }

        if ( std::abs(std::abs(scale_factor_third_dimension_)-1.0) > 1e-6 )
        {
            GADGET_CHECK_RETURN_FALSE(this->thirdDimensionScale(res_after_apply_kernel_sum_over_, scale_factor_third_dimension_));
        }

        GADGET_CHECK_RETURN_FALSE(this->L1NormTotal(res_after_apply_kernel_sum_over_, wav_coeff_norm_, obj));
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::objTask(const hoNDArray<T>& x, T& obj) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
firstDimensionScale(hoNDArray<T>& wavCoeff, T& scaleFactor)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t E1 = wavCoeff.get_size(1);
        size_t E2 = wavCoeff.get_size(2);
        size_t W = wavCoeff.get_size(3);

        size_t num = wavCoeff.get_number_of_elements()/(RO*E1*E2*W);

        // coeff 2, 3, 6, 7 are for RO high frequency

        size_t ii;
        for ( ii=0; ii<num; ii++ )
        {
            for ( size_t n=0; n<numOfWavLevels_; n++ )
            {
                // 2, 3
                hoNDArray<T> coeff(RO, E1, E2, 2, wavCoeff.begin()+ii*RO*E1*E2*W+(7*n+2)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff);

                // 6, 7
                hoNDArray<T> coeff2(RO, E1, E2, 2, wavCoeff.begin()+ii*RO*E1*E2*W+(7*n+6)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff2);
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::firstDimensionScale(hoNDArray<T>& wavCoeff, T& scaleFactor) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
secondDimensionScale(hoNDArray<T>& wavCoeff, T& scaleFactor)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t E1 = wavCoeff.get_size(1);
        size_t E2 = wavCoeff.get_size(2);
        size_t W = wavCoeff.get_size(3);

        size_t num = wavCoeff.get_number_of_elements()/(RO*E1*E2*W);

        // coeff 1, 3, 5, 7 are for E1 high frequency

        size_t ii;
        for ( ii=0; ii<num; ii++ )
        {
            for ( size_t n=0; n<numOfWavLevels_; n++ )
            {
                hoNDArray<T> coeff(RO, E1, E2, 1, wavCoeff.begin()+ii*RO*E1*E2*W+(7*n+1)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff);

                hoNDArray<T> coeff1(RO, E1, E2, 1, wavCoeff.begin()+ii*RO*E1*E2*W+(7*n+3)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff1);

                hoNDArray<T> coeff2(RO, E1, E2, 1, wavCoeff.begin()+ii*RO*E1*E2*W+(7*n+5)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff2);

                hoNDArray<T> coeff3(RO, E1, E2, 1, wavCoeff.begin()+ii*RO*E1*E2*W+(7*n+7)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff3);
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::secondDimensionScale(hoNDArray<T>& wavCoeff, T& scaleFactor) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusWavelet3DOperator<T>::
thirdDimensionScale(hoNDArray<T>& wavCoeff, T& scaleFactor)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t E1 = wavCoeff.get_size(1);
        size_t E2 = wavCoeff.get_size(2);
        size_t W = wavCoeff.get_size(3);

        size_t num = wavCoeff.get_number_of_elements()/(RO*E1*E2*W);

        // coeff 4, 5, 6, 7 are for E2 high frequency
        size_t ii;
        for ( ii=0; ii<num; ii++ )
        {
            for ( size_t n=0; n<numOfWavLevels_; n++ )
            {
                hoNDArray<T> coeff(RO, E1, E2, 4, wavCoeff.begin()+ii*RO*E1*E2*W+(7*n+4)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff);
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusWavelet3DOperator<T>::thirdDimensionScale(hoNDArray<T>& wavCoeff, T& scaleFactor) ... ");
        return false;
    }

    return true;
}

template <typename T> 
void gtPlusWavelet3DOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD wavelet 3D operator -----------------------" << endl;
    os << "Wavelet operator for gtPlus ISMRMRD package" << endl;
    os << "-------------------------------------------------------------------------" << endl;
}

}}
