
#include "hoWavelet2DTOperator.h"
#include "mri_core_coil_map_estimation.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"

#ifdef min
    #undef min
#endif // min

namespace Gadgetron
{

template <typename T>
hoWavelet2DTOperator<T>::hoWavelet2DTOperator(std::vector<size_t> *dims) : BaseClass(dims)
{
    scale_factor_first_dimension_ = 1;
    scale_factor_second_dimension_ = 1;
    scale_factor_third_dimension_ = 1;

    change_coeffcients_third_dimension_boundary_ = true;
}

template <typename T>
hoWavelet2DTOperator<T>::~hoWavelet2DTOperator()
{
}

template <typename T>
void hoWavelet2DTOperator<T>::convert_to_image(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    if (!complexIm_fft_.dimensions_equal(&x))
    {
        complexIm_fft_.create(x.dimensions());
    }

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(x, im, complexIm_fft_);
}

template <typename T>
void hoWavelet2DTOperator<T>::convert_to_kspace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if (!kspace_fft_.dimensions_equal(&im))
    {
        kspace_fft_.create(im.dimensions());
    }

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(im, x, kspace_fft_);
}

template <typename T>
void hoWavelet2DTOperator<T>::forward_wav(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        auto dims = x.dimensions();
        size_t NDim = dims.size();

        size_t RO = dims[0];
        size_t E1 = dims[1];
        size_t CHA = dims[2];
        size_t E2 = dims[3];
        size_t W = 1 + 7 * num_of_wav_levels_;

        std::vector<size_t> dimR(NDim + 1);
        dimR[0] = RO;
        dimR[1] = E1;
        dimR[2] = E2;
        dimR[3] = W;
        dimR[4] = CHA;

        size_t n;
        for (n = 4; n<NDim; n++)
        {
            dimR[n + 1] = dims[n];
        }

        if (!y.dimensions_equal(&dimR))
        {
            y.create(dimR);
        }

        size_t num = x.get_number_of_elements() / (RO*E1*E2*CHA);

        T* pX = const_cast<T*>(x.begin());
        T* pY = y.begin();

        int t;

        if (CHA == 1)
        {
#pragma omp parallel for default(none) private(t) shared(num, RO, E1, E2, W, pX, pY) if ( num > 1)
            for (t = 0; t<num; t++)
            {
                hoNDArray<T> in(RO, E1, E2, pX + t*RO*E1*E2);
                hoNDArray<T> out(RO, E1, E2, W, pY + t*RO*E1*E2*W);
                p_active_wav_->transform(in, out, 3, num_of_wav_levels_, true);
            }
        }
        else
        {
            {
                forward_buf_.create(RO, E1, E2, CHA);

                std::vector<size_t> dimOrder(4);
                dimOrder[0] = 0;
                dimOrder[1] = 1;
                dimOrder[2] = 3;
                dimOrder[3] = 2;

                for (t = 0; t<num; t++)
                {
                    hoNDArray<T> in(RO, E1, CHA, E2, pX + t*RO*E1*CHA*E2);
                    Gadgetron::permute(in, forward_buf_, dimOrder);

                    long long cha;

#pragma omp parallel for default(none) private(cha) shared(num, RO, E1, CHA, E2, W, pY, t) if ( CHA > 4 )
                    for (cha = 0; cha<CHA; cha++)
                    {
                        hoNDArray<T> in_dwt(RO, E1, E2, forward_buf_.begin() + cha*RO*E1*E2);
                        hoNDArray<T> out(RO, E1, E2, W, pY + t*RO*E1*E2*W*CHA + cha*RO*E1*E2*W);
                        p_active_wav_->transform(in_dwt, out, 3, num_of_wav_levels_, true);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::forward_wav(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::adjoint_wav(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        auto dims = x.dimensions();
        size_t NDim = dims.size();

        size_t RO = dims[0];
        size_t E1 = dims[1];
        size_t E2 = dims[2];
        size_t W = dims[3];
        size_t CHA = dims[4];

        std::vector<size_t> dimR(NDim - 1);
        dimR[0] = RO;
        dimR[1] = E1;
        dimR[2] = CHA;
        dimR[3] = E2;

        size_t n;
        for (n = 4; n<NDim - 1; n++)
        {
            dimR[n] = dims[n + 1];
        }

        if (!y.dimensions_equal(&dimR))
        {
            y.create(dimR);
        }

        size_t num = x.get_number_of_elements() / (RO*E1*E2*W*CHA);

        T* pX = const_cast<T*>(x.begin());
        T* pY = y.begin();

        int t;

        if (CHA == 1)
        {
#pragma omp parallel for default(none) private(t) shared(num, RO, E1, E2, W, pX, pY) if ( num > 1)
            for (t = 0; t<num; t++)
            {
                hoNDArray<T> in(RO, E1, E2, W, pX + t*RO*E1*E2*W);
                hoNDArray<T> out(RO, E1, E2, pY + t*RO*E1*E2);
                p_active_wav_->transform(in, out, 3, num_of_wav_levels_, false);
            }
        }
        else
        {
            {
                adjoint_buf_.create(RO, E1, E2, CHA);

                std::vector<size_t> dimOrder(4);
                dimOrder[0] = 0;
                dimOrder[1] = 1;
                dimOrder[2] = 3;
                dimOrder[3] = 2;

                for (t = 0; t<num; t++)
                {
                    hoNDArray<T> out(RO, E1, CHA, E2, pY + t*RO*E1*CHA*E2);

                    long long cha;
#pragma omp parallel for default(none) private(cha) shared(RO, E1, CHA, E2, W, pX) if ( CHA > 4 )
                    for (cha = 0; cha<CHA; cha++)
                    {
                        hoNDArray<T> in(RO, E1, E2, W, pX + cha*RO*E1*E2*W);
                        hoNDArray<T> out_idwt(RO, E1, E2, adjoint_buf_.begin() + cha*RO*E1*E2);
                        p_active_wav_->transform(in, out_idwt, 3, num_of_wav_levels_, false);
                    }

                    Gadgetron::permute(adjoint_buf_, out, dimOrder);
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::apply_scale_first_dimension(hoNDArray<T>& wavCoeff, value_type& scaleFactor)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t E1 = wavCoeff.get_size(1);
        size_t E2 = wavCoeff.get_size(2);
        size_t W = wavCoeff.get_size(3);

        size_t num = wavCoeff.get_number_of_elements() / (RO*E1*E2*W);

        // coeff 2, 3, 6, 7 are for RO high frequency

        size_t ii;
        for (ii = 0; ii<num; ii++)
        {
            for (size_t n = 0; n<num_of_wav_levels_; n++)
            {
                // 2, 3
                hoNDArray<T> coeff(RO, E1, E2, 2, wavCoeff.begin() + ii*RO*E1*E2*W + (7 * n + 2)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff);

                // 6, 7
                hoNDArray<T> coeff2(RO, E1, E2, 2, wavCoeff.begin() + ii*RO*E1*E2*W + (7 * n + 6)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff2);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::apply_scale_first_dimension(...) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::apply_scale_second_dimension(hoNDArray<T>& wavCoeff, value_type& scaleFactor)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t E1 = wavCoeff.get_size(1);
        size_t E2 = wavCoeff.get_size(2);
        size_t W = wavCoeff.get_size(3);

        size_t num = wavCoeff.get_number_of_elements() / (RO*E1*E2*W);

        // coeff 1, 3, 5, 7 are for E1 high frequency

        size_t ii;
        for (ii = 0; ii<num; ii++)
        {
            for (size_t n = 0; n<num_of_wav_levels_; n++)
            {
                hoNDArray<T> coeff(RO, E1, E2, 1, wavCoeff.begin() + ii*RO*E1*E2*W + (7 * n + 1)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff);

                hoNDArray<T> coeff1(RO, E1, E2, 1, wavCoeff.begin() + ii*RO*E1*E2*W + (7 * n + 3)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff1);

                hoNDArray<T> coeff2(RO, E1, E2, 1, wavCoeff.begin() + ii*RO*E1*E2*W + (7 * n + 5)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff2);

                hoNDArray<T> coeff3(RO, E1, E2, 1, wavCoeff.begin() + ii*RO*E1*E2*W + (7 * n + 7)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff3);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::apply_scale_second_dimension(...) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::apply_scale_third_dimension(hoNDArray<T>& wavCoeff, value_type& scaleFactor)
{
    try
    {
        size_t RO = wavCoeff.get_size(0);
        size_t E1 = wavCoeff.get_size(1);
        size_t E2 = wavCoeff.get_size(2);
        size_t W = wavCoeff.get_size(3);

        size_t num = wavCoeff.get_number_of_elements() / (RO*E1*E2*W);

        // coeff 4, 5, 6, 7 are for E2 high frequency
        size_t ii;
        for (ii = 0; ii<num; ii++)
        {
            for (size_t n = 0; n<num_of_wav_levels_; n++)
            {
                hoNDArray<T> coeff(RO, E1, E2, 4, wavCoeff.begin() + ii*RO*E1*E2*W + (7 * n + 4)*RO*E1*E2);
                Gadgetron::scal(scaleFactor, coeff);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::apply_scale_third_dimension(...) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
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

        //GDEBUG_STREAM("complexIm_ = " << Gadgetron::norm2(complexIm_));

        size_t RO = x->get_size(0);
        size_t E1 = x->get_size(1);
        size_t CHA = x->get_size(2);
        size_t N = x->get_size(3);

        std::vector<size_t> dimMOCO(3);
        dimMOCO[0] = RO;
        dimMOCO[1] = E1;
        dimMOCO[2] = N;

        bool hasMOCO = false;
        if (mocoer_.dx_.dimensions_equal(&dimMOCO) && mocoer_.dy_.dimensions_equal(&dimMOCO))
        {
            hasMOCO = true;
        }

        if(coil_map_.get_size(0) == RO && coil_map_.get_size(1) == E1 && coil_map_.get_size(2) == CHA)
        {
            // S'
            Gadgetron::coil_combine(complexIm_, coil_map_, 2, complexIm_after_apply_coil_map_);

            hoNDArray<T> combined(RO, E1, 1, N, complexIm_after_apply_coil_map_.begin());

            //GDEBUG_STREAM("combined = " << Gadgetron::norm2(combined));

            if(hasMOCO)
            {
                // M
                mocoer_.mult_M(&combined, &moco_image_, false);

                //GDEBUG_STREAM("moco_image_ = " << Gadgetron::norm2(moco_image_));

                // W
                this->forward_wav(moco_image_, *y);
            }
            else
            {
                this->forward_wav(combined, *y);
            }
        }
        else
        {
            if (hasMOCO)
            {
                mocoer_.mult_M(&complexIm_, &moco_image_, false);
                this->forward_wav(moco_image_, *y);
            }
            else
            {
                this->forward_wav(complexIm_, *y);

                //GDEBUG_STREAM("y = " << Gadgetron::norm2(*y));
            }
        }

        if (accumulate)
        {
            Gadgetron::add(wav_coeff_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::mult_M(...) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            complexIm_ = *y;
        }

        size_t RO = x->get_size(0);
        size_t E1 = x->get_size(1);
        size_t N = x->get_size(2);
        size_t W = x->get_size(3);
        size_t CHA = x->get_size(4);

        std::vector<size_t> dimR(4);
        dimR[0] = RO;
        dimR[1] = E1;
        dimR[2] = CHA;
        dimR[3] = N;

        bool hasCoilMap = false;
        if (coil_map_.get_size(0) == RO && coil_map_.get_size(1) == E1)
        {
            hasCoilMap = true;
        }

        if (hasCoilMap)
        {
            dimR[2] = coil_map_.get_size(2);
        }

        std::vector<size_t> dimMOCO(3);
        dimMOCO[0] = RO;
        dimMOCO[1] = E1;
        dimMOCO[2] = N;

        bool hasMOCO = false;
        if (mocoer_.adj_dx_.dimensions_equal(&dimMOCO) && mocoer_.adj_dy_.dimensions_equal(&dimMOCO))
        {
            hasMOCO = true;
        }

        if (!y->dimensions_equal(&dimR))
        {
            y->create(&dimR);
        }

        // W'
        this->adjoint_wav(*x, this->complexIm_wav_);

        // M'
        if(hasMOCO)
        {
            mocoer_.mult_MH(&complexIm_wav_, &adj_moco_image_, false);
            this->complexIm_wav_ = this->adj_moco_image_;
        }

        if (hasCoilMap)
        {
            // S
            if (kspace_wav_.get_number_of_elements() != RO*E1*dimR[2] * N)
            {
                kspace_wav_.create(RO, E1, dimR[2], N);
            }

            size_t CHA = dimR[2];

            size_t n;
            for (size_t n = 0; n<N; n++)
            {
                hoNDArray<T> complexImN(RO, E1, complexIm_wav_.begin() + n*RO*E1);
                hoNDArray<T> kspace_wavN(RO, E1, CHA, kspace_wav_.begin() + n*RO*E1*CHA);

                if (coil_map_.get_size(3) == N)
                {
                    hoNDArray<T> coilMapN(RO, E1, CHA, coil_map_.begin() + n*RO*E1*CHA);
                    Gadgetron::multiply(coilMapN, complexImN, kspace_wavN);
                }
                else
                {
                    Gadgetron::multiply(coil_map_, complexImN, kspace_wavN);
                }
            }

            if (input_in_kspace_)
            {
                // F
                this->convert_to_kspace(this->kspace_wav_, *y);

                if (!no_null_space_)
                {
                    // Dc
                    Gadgetron::multiply(unacquired_points_indicator_, *y, *y);
                }
            }
            else
            {
                *y = this->kspace_wav_;
            }
        }
        else
        {
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
        }

        if (accumulate)
        {
            Gadgetron::add(complexIm_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::mult_MH(...) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate)
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

        // apply the scaling factor
        if (std::abs(scale_factor_first_dimension_ - 1.0) > 1e-6)
        {
            this->apply_scale_first_dimension(wav_coeff_, scale_factor_first_dimension_);
        }

        // second dimension scaling
        if (std::abs(scale_factor_second_dimension_ - 1.0) > 1e-6)
        {
            this->apply_scale_second_dimension(wav_coeff_, scale_factor_second_dimension_);
        }

        // third dimension scaling
        if (std::abs(scale_factor_third_dimension_ - 1.0) > 1e-6)
        {
            this->apply_scale_third_dimension(wav_coeff_, scale_factor_third_dimension_);
        }

        this->mult_MH(&wav_coeff_, g, false);

        if (accumulate)
        {
            Gadgetron::add(complexIm_wav_, *g, *g);
        }
    }
    catch(...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::gradient(...) ... ");
    }
}

template <typename T>
typename hoWavelet2DTOperator<T>::REAL hoWavelet2DTOperator<T>::magnitude(ARRAY_TYPE* x)
{
    try
    {
        this->mult_M(x, &wav_coeff_, false);

        if (std::abs(scale_factor_first_dimension_ - 1.0) > 1e-6)
        {
            this->apply_scale_first_dimension(wav_coeff_, scale_factor_first_dimension_);
        }

        if (std::abs(scale_factor_second_dimension_ - 1.0) > 1e-6)
        {
            this->apply_scale_second_dimension(wav_coeff_, scale_factor_second_dimension_);
        }

        if (std::abs(scale_factor_third_dimension_ - 1.0) > 1e-6)
        {
            this->apply_scale_third_dimension(wav_coeff_, scale_factor_third_dimension_);
        }

        this->L1Norm(wav_coeff_, wav_coeff_norm_);

        typename hoWavelet2DTOperator<T>::REAL L1CoeffNorm = Gadgetron::asum(&wav_coeff_norm_);

        return L1CoeffNorm;
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::magnitude(...) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::proximity(hoNDArray<T>& wavCoeff, value_type thres)
{
    try
    {
        if (this->proximity_across_cha_)
        {
            this->L1Norm(wavCoeff, wav_coeff_norm_);
            Gadgetron::abs(wav_coeff_norm_, wav_coeff_norm_mag_);
        }
        else
        {
            Gadgetron::abs(wavCoeff, wav_coeff_norm_mag_);
        }

        if (!mask_.dimensions_equal(&wavCoeff))
        {
            mask_.create(wavCoeff.dimensions());
        }

        Gadgetron::fill(mask_, T(thres));

        if (std::abs(scale_factor_first_dimension_ - 1.0) > 1e-6)
        {
            this->apply_scale_first_dimension(mask_, scale_factor_first_dimension_);
        }

        if (std::abs(scale_factor_second_dimension_ - 1.0) > 1e-6)
        {
            this->apply_scale_second_dimension(mask_, scale_factor_second_dimension_);
        }

        if (std::abs(scale_factor_third_dimension_ - 1.0) > 1e-6)
        {
            this->apply_scale_third_dimension(mask_, scale_factor_third_dimension_);
        }

        this->shrink_wav_coeff(wavCoeff, wav_coeff_norm_mag_, mask_, with_approx_coeff_);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::proximity(...) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<value_type>& wavCoeffNorm)
{
    try
    {
        auto dims = wavCoeff.dimensions();

        auto dimR = dims;
        dimR[4] = 1;

        if (!wavCoeffNorm.dimensions_equal(&dimR))
        {
            wavCoeffNorm.create(dimR);
        }

        size_t RO = dims[0];
        size_t E1 = dims[1];
        size_t E2 = dims[2];
        size_t W = dims[3];
        size_t CHA = dims[4];

        if (CHA > 1)
        {
            if (!complexIm_norm_.dimensions_equal(&wavCoeff))
            {
                complexIm_norm_.create(wavCoeff.dimensions());
            }

            size_t n;
            size_t N = wavCoeff.get_number_of_elements();
            const T* pWavCoeff = wavCoeff.begin();
            for (n = 0; n < N; n++)
            {
                /*T v = pWavCoeff[n] * std::conj(pWavCoeff[n]);
                complexIm_norm_(n) = std::abs(v);*/

                value_type v = std::abs(pWavCoeff[n]);
                complexIm_norm_(n) = v*v;
            }

            Gadgetron::sum_over_dimension(complexIm_norm_, wavCoeffNorm, 4);
            GADGET_CATCH_THROW(Gadgetron::sqrt(wavCoeffNorm, wavCoeffNorm));
        }
        else
        {
            GADGET_CATCH_THROW(Gadgetron::abs(wavCoeff, wavCoeffNorm));
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet2DTOperator<T>::L1Norm(const hoNDArray<T>& wavCoeff, hoNDArray<T>& wavCoeffNorm) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::shrink_wav_coeff(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, const hoNDArray<T>& mask, bool processApproxCoeff)
{
    try
    {
        auto dims = wavCoeff.dimensions();

        long long RO = (long long)dims[0];
        long long E1 = (long long)dims[1];
        long long E2 = (long long)dims[2];
        long long W = (long long)dims[3];
        long long CHA = (long long)dims[4];

        T* pCoeff = wavCoeff.begin();

        long long CHA_NORM = wavCoeffNorm.get_size(4);
        value_type* pCoeffNorm = const_cast<value_type*>(wavCoeffNorm.begin());

        const T* pMask = mask.begin();

        long long startW = 1;
        if (processApproxCoeff)
        {
            startW = 0;
        }

        size_t num = wavCoeffNorm.get_number_of_elements();

        long long cha;

#pragma omp parallel for default(none) private(cha) shared(CHA, CHA_NORM, startW, W, RO, E1, E2, pCoeffNorm, pMask, pCoeff) if(CHA>1)
        for (cha = 0; cha < CHA; cha++)
        {
            size_t cha_norm = cha;
            if (cha_norm >= CHA_NORM) cha_norm = CHA_NORM - 1;

            long long w, e2, e1, ro;
            for (w = startW; w < W; w++)
            {
                for (e2 = 0; e2 < E2; e2++)
                {
                    for (e1 = 0; e1 < E1; e1++)
                    {
                        for (ro = 0; ro < RO; ro++)
                        {
                            size_t offset = ro + e1*RO + e2*RO*E1 + w*RO*E1*E2;
                            size_t ind = offset + cha*RO*E1*E2*W;
                            size_t indNorm = offset + cha_norm*RO*E1*E2*W;

                            if (pCoeffNorm[indNorm] < pMask[indNorm].real())
                            {
                                pCoeff[ind] = 0;
                            }
                            else
                            {
                                T v = pCoeff[ind];
                                value_type m = std::abs(v);
                                if (m > FLT_EPSILON)
                                {
                                    T v2 = v / m;
                                    m -= pMask[indNorm].real();
                                    pCoeff[ind] = m*v2;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoWavelet2DTOperator<T>::shrink_wav_coeff(hoNDArray<T>& wavCoeff, const hoNDArray<T>& wavCoeffNorm, T thres, const hoNDArray<T>& mask, bool processApproxCoeff) ... ");
    }
}

template <typename T>
void hoWavelet2DTOperator<T>::divide_wav_coeff_by_norm(hoNDArray<T>& wavCoeff, const hoNDArray<value_type>& wavCoeffNorm, value_type mu, value_type p, bool processApproxCoeff)
{
    try
    {
        long long RO = (long long)wavCoeff.get_size(0);
        long long E1 = (long long)wavCoeff.get_size(1);
        long long E2 = (long long)wavCoeff.get_size(2);
        long long W = (long long)wavCoeff.get_size(3);
        long long CHA = (long long)wavCoeff.get_size(4);

        long long CHA_NORM = (long long)wavCoeffNorm.get_size(4);

        if (!wav_coeff_norm_approx_.dimensions_equal(&wavCoeffNorm))
        {
            wav_coeff_norm_approx_.create(wavCoeffNorm.dimensions());
        }

        long long ii;
        long long N = (long long)wavCoeffNorm.get_number_of_elements();

        const value_type* pCoeffNorm = wavCoeffNorm.begin();
        value_type* pBuf = wav_coeff_norm_approx_.begin();

        if (std::abs(p - 1.0) < 0.001)
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

        if (processApproxCoeff)
        {
            Gadgetron::multiply(wavCoeff, wav_coeff_norm_approx_, wavCoeff);
        }
        else
        {
            for (long long cha = 0; cha<CHA; cha++)
            {
                size_t ind = (cha<CHA_NORM) ? cha : CHA_NORM-1;
                hoNDArray<value_type> wavCoeffNormCurr(RO, E1, E2, W - 1, wav_coeff_norm_approx_.begin() + ind*RO*E1*E2*W + RO*E1*E2);

                hoNDArray<T> wavCoeffCurr(RO, E1, E2, W - 1, wavCoeff.begin() + cha*RO*E1*E2*W + RO*E1*E2);

                Gadgetron::multiply(wavCoeffNormCurr, wavCoeffCurr, wavCoeffCurr);
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in hoWavelet2DTOperator<T>::divide_wav_coeff_by_norm(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCPUOPERATOR hoWavelet2DTOperator< std::complex<float> >;
template class EXPORTCPUOPERATOR hoWavelet2DTOperator< std::complex<double> >;

}
