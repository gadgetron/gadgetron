
#include "hoWavelet1DOperator.h"
#include "hoNDFFT.h"

#ifdef min
    #undef min
#endif // min

namespace Gadgetron
{

template <typename T>
hoWavelet1DBaseOperator<T>::hoWavelet1DBaseOperator(const std::vector<size_t>& dims) : BaseClass(dims)
{
}

template <typename T>
hoWavelet1DBaseOperator<T>::~hoWavelet1DBaseOperator()
{
}

template <typename T>
void hoWavelet1DBaseOperator<T>::forward_wav(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        size_t RO = x.get_size(0);
        size_t CHA = x.get_size(1);
        size_t W = 1 + num_of_wav_levels_;

        size_t num = x.get_number_of_elements() / (RO*CHA);

        if (y.get_size(0) != RO || y.get_size(1) != W || y.get_size(2) != CHA)
        {
            std::vector<size_t> dims = x.get_dimensions();
            size_t NDim = dims.size();

            std::vector<size_t> dimR(NDim + 1);
            dimR[0] = RO;
            dimR[1] = W;
            dimR[2] = CHA;

            size_t n;
            for (n = 2; n<NDim; n++)
            {
                dimR[n + 1] = dims[n];
            }

            if (!y.dimensions_equal(dimR))
            {
                y.create(dimR);
            }
        }

        p_active_wav_->transform(x, y, 1, num_of_wav_levels_, true);

        if (wav_coeff_scaling_.get_size(0) == RO && wav_coeff_scaling_.get_size(1) == W)
        {
            size_t n, c;
            for (n = 0; n < num; n++)
            {
                for (c = 0; c < CHA; c++)
                {
                    T* pY = y.begin() + n*RO*W*CHA + c*RO*W;
                    for (size_t ii = 0; ii < RO*W; ii++)
                    {
                        pY[ii] *= wav_coeff_scaling_(ii);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet1DBaseOperator<T>::forward_wav(...) ... ");
    }
}

template <typename T>
void hoWavelet1DBaseOperator<T>::adjoint_wav(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        size_t RO = x.get_size(0);
        size_t W = x.get_size(1);
        size_t CHA = x.get_size(2);

        size_t num = x.get_number_of_elements() / (RO*W*CHA);

        if (!adj_x_.dimensions_equal(x))
        {
            adj_x_ = x;
        }

        memcpy(adj_x_.begin(), x.begin(), x.get_number_of_bytes());

        if (wav_coeff_scaling_.get_size(0) == RO && wav_coeff_scaling_.get_size(1) == W)
        {
            size_t n, c;
            for (n = 0; n < num; n++)
            {
                for (c = 0; c < CHA; c++)
                {
                    T* pX = adj_x_.begin() + n*RO*W*CHA + c*RO*W;
                    for (size_t ii = 0; ii < RO*W; ii++)
                    {
                        pX[ii] /= wav_coeff_scaling_(ii);
                    }
                }
            }
        }

        if (y.get_size(0) != RO || y.get_size(1) != CHA)
        {
            std::vector<size_t> dims = x.get_dimensions();
            size_t NDim = dims.size();

            std::vector<size_t> dimR(NDim - 1);
            dimR[0] = RO;
            dimR[1] = CHA;

            size_t n;
            for (n = 2; n < NDim - 1; n++)
            {
                dimR[n] = dims[n + 1];
            }

            if (!y.dimensions_equal(dimR))
            {
                y.create(dimR);
            }
        }

        p_active_wav_->transform(adj_x_, y, 1, num_of_wav_levels_, false);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet1DBaseOperator<T>::adjoint_wav(...) ... ");
    }
}

template <typename T>
void hoWavelet1DBaseOperator<T>::gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            complexIm_wav_ = *g;
        }

        this->mult_M(x, &wav_coeff_, false);
        this->L1Norm(wav_coeff_, wav_coeff_norm_);
        this->divide_wav_coeff_by_norm(wav_coeff_, wav_coeff_norm_, (value_type)(1e-15), (value_type)(1.0), with_approx_coeff_);
        this->mult_MH(&wav_coeff_, g, false);

        if (accumulate)
        {
            Gadgetron::add(complexIm_wav_, *g, *g);
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in hoWavelet1DBaseOperator<T>::gradient(...) ... ");
    }
}

template <typename T>
typename hoWavelet1DBaseOperator<T>::REAL hoWavelet1DBaseOperator<T>::magnitude(ARRAY_TYPE* x)
{
    try
    {
        this->mult_M(x, &wav_coeff_, false);
        this->L1Norm(wav_coeff_, wav_coeff_norm_);

        typename hoWavelet1DBaseOperator<T>::REAL L1CoeffNorm = Gadgetron::asum(&wav_coeff_norm_);

        return L1CoeffNorm;
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet1DBaseOperator<T>::magnitude(...) ... ");
    }
}

template <typename T>
void hoWavelet1DBaseOperator<T>::proximity(hoNDArray<T>& wavCoeff, value_type thres)
{
    try
    {
        this->L1Norm(wavCoeff, wav_coeff_norm_);
        this->shrink_wav_coeff(wavCoeff, wav_coeff_norm_, thres, with_approx_coeff_);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet1DBaseOperator<T>::proximity(...) ... ");
    }
}

template <typename T>
void hoWavelet1DBaseOperator<T>::L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<value_type>& wavCoeffNorm)
{
    try
    {
        std::vector<size_t> dims = wavCoeff.get_dimensions();

        std::vector<size_t> dimR(dims);
        dimR[2] = 1;

        if (!wavCoeffNorm.dimensions_equal(dimR))
        {
            wavCoeffNorm.create(dimR);
        }

        size_t RO = dims[0];
        size_t W = dims[1];
        size_t CHA = dims[2];

        if (CHA > 1)
        {
            // square the coefficients
            if (!complexIm_norm_.dimensions_equal(wavCoeff))
            {
                complexIm_norm_.create(wavCoeff.dimensions());
            }

            size_t n;
            size_t N = wavCoeff.get_number_of_elements();
            const T* pWavCoeff = wavCoeff.begin();
            for (n = 0; n < N; n++)
            {
                value_type v = std::abs(pWavCoeff[n]);
                complexIm_norm_(n) = v*v;
            }

            // Gadgetron::multiplyConj(wavCoeff, wavCoeff, complexIm_norm_);

            // sum over CHA
            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(complexIm_norm_, wavCoeffNorm, 2));

            GADGET_CATCH_THROW(Gadgetron::sqrt(wavCoeffNorm, wavCoeffNorm));
        }
        else
        {
            GADGET_CATCH_THROW(Gadgetron::abs(wavCoeff, wavCoeffNorm));
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet1DBaseOperator<T>::L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm) ... ");
    }
}

template <typename T>
void hoWavelet1DBaseOperator<T>::shrink_wav_coeff(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, value_type thres, bool processApproxCoeff)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t W = wavCoeff.get_size(1);
        size_t CHA = wavCoeff.get_size(2);

        size_t num = wavCoeff.get_number_of_elements() / (RO*W*CHA);

        size_t ro, w, cha;

        size_t sw = 0;
        if (!processApproxCoeff) sw = 1;

        value_type minV = std::numeric_limits<value_type>::min();

        long long n(0);

#pragma omp parallel for default(none) private(n, w, ro, cha) shared(num, wavCoeffNorm, wavCoeff, RO, W, CHA, minV, thres, sw)
        for (n = 0; n < (long long)num; n++)
        {
            const value_type* pWavCoeffNorm = wavCoeffNorm.begin() + n*RO*W;
            T* pWavCoeff = wavCoeff.begin() + n*RO*W*CHA;

            for (w = sw; w < W; w++)
            {
                for (ro = 0; ro < RO; ro++)
                {
                    if ( std::abs(pWavCoeffNorm[ro+w*RO])< thres)
                    {
                        for (cha = 0; cha < CHA; cha++)
                        {
                            pWavCoeff[ro+w*RO+cha*RO*W] = 0;
                        }
                    }
                    else
                    {
                        for (cha = 0; cha < CHA; cha++)
                        {

                            T v = pWavCoeff[ro + w*RO + cha*RO*W];
                            value_type m = std::abs(v);
                            if (m < minV) m = minV;
                            T v2 = v / m;
                            m -= thres;
                            pWavCoeff[ro + w*RO + cha*RO*W] = m*v2;
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoWavelet1DBaseOperator<T>::shrink_wav_coeff(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T thres, const hoNDArray<T>& mask, bool processApproxCoeff) ... ");
    }
}

template <typename T>
void hoWavelet1DBaseOperator<T>::divide_wav_coeff_by_norm(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, value_type mu, value_type p, bool processApproxCoeff)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t W = wavCoeff.get_size(1);
        size_t CHA = wavCoeff.get_size(2);

        if (!wav_coeff_norm_approx_.dimensions_equal(wavCoeffNorm))
        {
            wav_coeff_norm_approx_.create(wavCoeffNorm.dimensions());
        }

        long long N = (long long)wavCoeffNorm.get_number_of_elements();

        const value_type* pCoeffNorm = wavCoeffNorm.begin();
        value_type* pBuf = wav_coeff_norm_approx_.begin();

        long long ii(0);

        if (std::abs(std::abs(p) - 1.0) < 0.001)
        {
#pragma omp parallel for default(none) private(ii) shared(N, pBuf, pCoeffNorm, mu)
            for (ii = 0; ii<N; ii++)
            {
                pBuf[ii] = (value_type)(1.0 / std::sqrt(pCoeffNorm[ii] * pCoeffNorm[ii] + mu));
            }
        }
        else
        {
#pragma omp parallel for default(none) private(ii) shared(N, pBuf, pCoeffNorm, mu, p)
            for (ii = 0; ii<N; ii++)
            {
                pBuf[ii] = (value_type)std::pow((double)(pCoeffNorm[ii] * pCoeffNorm[ii] + mu), (double)(p / 2.0 - 1.0));
            }
        }

        long long num = wavCoeff.get_number_of_elements() / (RO*W*CHA);

        if (processApproxCoeff)
        {
#pragma omp parallel default(none) private(ii) shared(RO, num, wavCoeffNorm, wavCoeff, W, CHA) if ( num > 1 )
            {

#pragma omp for
                for (ii = 0; ii<num; ii++)
                {
                    hoNDArray<value_type> wavCoeffNormCurr(RO, W, wav_coeff_norm_approx_.begin() + ii*RO*W);

                    for (size_t cha = 0; cha<CHA; cha++)
                    {
                        hoNDArray<T> wavCoeffCurr(RO, W, wavCoeff.begin() + ii*RO*W*CHA + cha*RO*W);
                        Gadgetron::multiply(wavCoeffNormCurr, wavCoeffCurr, wavCoeffCurr);
                    }
                }
            }
        }
        else
        {
#pragma omp parallel default(none) private(ii) shared(RO, num, wavCoeffNorm, wavCoeff, W, CHA) if ( num > 1 )
            {

#pragma omp for
                for (ii = 0; ii<num; ii++)
                {
                    hoNDArray<value_type> wavCoeffNormCurr(RO, W - 1, wav_coeff_norm_approx_.begin() + ii*RO*W + RO);

                    for (size_t cha = 0; cha<CHA; cha++)
                    {
                        hoNDArray<T> wavCoeffCurr(RO, W - 1, wavCoeff.begin() + ii*RO*W*CHA + cha*RO + RO);
                        Gadgetron::multiply(wavCoeffNormCurr, wavCoeffCurr, wavCoeffCurr);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoWavelet1DBaseOperator<T>::divide_wav_coeff_by_norm(...) ... ");
    }
}

// ------------------------------------------------------------

template <typename T>
hoWavelet1DOperator<T>::hoWavelet1DOperator(const std::vector<size_t>& dims) : BaseClass(dims)
{
}

template <typename T>
hoWavelet1DOperator<T>::~hoWavelet1DOperator()
{
}

template <typename T>
void hoWavelet1DOperator<T>::convert_to_image(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    if (!complexIm_fft_.dimensions_equal(x))
    {
        complexIm_fft_.create(x.dimensions());
    }

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(x, im, complexIm_fft_);
}

template <typename T>
void hoWavelet1DOperator<T>::convert_to_kspace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if (!kspace_fft_.dimensions_equal(im))
    {
        kspace_fft_.create(im.dimensions());
    }

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(im, x, kspace_fft_);
}

template <typename T>
void hoWavelet1DOperator<T>::mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            wav_coeff_ = *y;
        }

        if (input_in_kspace_)
        {
            if (!no_null_space_)
            {
                // Dc'x+D'a
                Gadgetron::multiply(unacquired_points_indicator_, *x, kspace_);
                Gadgetron::add(acquired_points_, kspace_, kspace_);

                // F'
                this->convert_to_image(kspace_, complexIm_);
            }
            else
            {
                this->convert_to_image(*x, complexIm_);
            }
        }
        else
        {
            complexIm_ = *x;
        }

        // W
        this->forward_wav(complexIm_, *y);

        if (accumulate)
        {
            Gadgetron::add(wav_coeff_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet1DOperator<T>::mult_M(...) ... ");
    }
}

template <typename T>
void hoWavelet1DOperator<T>::mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            complexIm_ = *y;
        }

        // W'
        this->adjoint_wav(*x, this->complexIm_wav_);

        if (input_in_kspace_)
        {
            // F
            this->convert_to_kspace(this->complexIm_wav_, *y);

            if (!no_null_space_)
            {
                // Dc
                Gadgetron::multiply(unacquired_points_indicator_, *y, *y);
            }
        }
        else
        {
            *y = this->complexIm_wav_;
        }

        if (accumulate)
        {
            Gadgetron::add(complexIm_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet1DOperator<T>::mult_MH(...) ... ");
    }
}

// ------------------------------------------------------------

template <typename T>
hoWavelet1DREALOperator<T>::hoWavelet1DREALOperator(const std::vector<size_t>& dims) : BaseClass(dims)
{
    input_in_kspace_ = false;
}

template <typename T>
hoWavelet1DREALOperator<T>::~hoWavelet1DREALOperator()
{
}

template <typename T>
void hoWavelet1DREALOperator<T>::mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            wav_coeff_ = *y;
        }

        // W
        this->forward_wav(*x, *y);

        if (accumulate)
        {
            Gadgetron::add(wav_coeff_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet1DREALOperator<T>::mult_M(...) ... ");
    }
}

template <typename T>
void hoWavelet1DREALOperator<T>::mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            complexIm_ = *y;
        }

        // W'
        this->adjoint_wav(*x, *y);

        if (accumulate)
        {
            Gadgetron::add(complexIm_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet1DREALOperator<T>::mult_MH(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class hoWavelet1DREALOperator< float >;
template class hoWavelet1DREALOperator< double >;

template class hoWavelet1DOperator< std::complex<float> >;
template class hoWavelet1DOperator< std::complex<double> >;

}
