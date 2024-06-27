
#include "hoWavelet3DOperator.h"
#include "mri_core_coil_map_estimation.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"

#ifdef min
    #undef min
#endif // min

namespace Gadgetron
{

template <typename T>
hoWavelet3DOperator<T>::hoWavelet3DOperator(std::vector<size_t> *dims) : BaseClass(dims)
{
}

template <typename T>
hoWavelet3DOperator<T>::~hoWavelet3DOperator()
{
}

template <typename T>
void hoWavelet3DOperator<T>::convert_to_image(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    if (!complexIm_fft_.dimensions_equal(x))
    {
        complexIm_fft_.create(x.dimensions());
    }

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(x, im, complexIm_fft_);
}

template <typename T>
void hoWavelet3DOperator<T>::convert_to_kspace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if (!kspace_fft_.dimensions_equal(im))
    {
        kspace_fft_.create(im.dimensions());
    }

    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(im, x, kspace_fft_);
}

template <typename T>
void hoWavelet3DOperator<T>::forward_wav(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        std::vector<size_t> dims = x.dimensions();
        size_t NDim = dims.size();

        size_t RO = dims[0];
        size_t E1 = dims[1];
        size_t E2 = dims[2];
        size_t CHA = dims[3];
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

        if (!y.dimensions_equal(dimR))
        {
            y.create(dimR);
        }

        p_active_wav_->transform(x, y, 3, num_of_wav_levels_, true);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet3DOperator<T>::forward_wav(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
    }
}

template <typename T>
void hoWavelet3DOperator<T>::adjoint_wav(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        std::vector<size_t> dims = x.dimensions();
        size_t NDim = dims.size();

        size_t RO = dims[0];
        size_t E1 = dims[1];
        size_t E2 = dims[2];
        size_t W = dims[3];
        size_t CHA = dims[4];

        std::vector<size_t> dimR(NDim - 1);
        dimR[0] = RO;
        dimR[1] = E1;
        dimR[2] = E2;
        dimR[3] = CHA;

        size_t n;
        for (n = 4; n<NDim - 1; n++)
        {
            dimR[n] = dims[n + 1];
        }

        if (!y.dimensions_equal(dimR))
        {
            y.create(dimR);
        }

        p_active_wav_->transform(x, y, 3, num_of_wav_levels_, false);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet3DOperator<T>::adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y) ... ");
    }
}

template <typename T>
void hoWavelet3DOperator<T>::mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
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

        size_t RO = x->get_size(0);
        size_t E1 = x->get_size(1);
        size_t E2 = x->get_size(2);
        size_t CHA = x->get_size(3);

        if(coil_map_.get_size(0) == RO && coil_map_.get_size(1) == E1 && coil_map_.get_size(2) == E2 && coil_map_.get_size(2) == CHA)
        {
            Gadgetron::coil_combine(complexIm_, coil_map_, 3, complexIm_after_apply_coil_map_);

            hoNDArray<T> combined(RO, E1, E2, 1, complexIm_after_apply_coil_map_.begin());
            this->forward_wav(combined, *y);
        }
        else
        {
            this->forward_wav(complexIm_, *y);
        }

        if (accumulate)
        {
            Gadgetron::add(wav_coeff_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoWavelet3DOperator<T>::mult_M(...) ... ");
    }
}

template <typename T>
void hoWavelet3DOperator<T>::mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            complexIm_ = *y;
        }

        size_t RO = x->get_size(0);
        size_t E1 = x->get_size(1);
        size_t E2 = x->get_size(2);
        size_t W = x->get_size(3);
        size_t CHA = x->get_size(4);

        std::vector<size_t> dimR(4);
        dimR[0] = RO;
        dimR[1] = E1;
        dimR[2] = E2;
        dimR[3] = CHA;

        bool hasCoilMap = false;
        if (coil_map_.get_size(0) == RO && coil_map_.get_size(1) == E1 && coil_map_.get_size(2) == E2)
        {
            hasCoilMap = true;
        }

        if (hasCoilMap)
        {
            dimR[3] = coil_map_.get_size(3);
            GADGET_CHECK_THROW(CHA==1); // already coil combined
        }

        if (!y->dimensions_equal(dimR))
        {
            y->create(dimR);
        }

        this->adjoint_wav(*x, this->complexIm_wav_);

        if (hasCoilMap)
        {
            if (kspace_wav_.get_number_of_elements() != RO*E1*E2*dimR[3])
            {
                kspace_wav_.create(RO, E1, E2, dimR[3]);
            }

            hoNDArray<T> combined;
            combined.create(RO, E1, E2, this->complexIm_wav_.begin());
            Gadgetron::multiply(coil_map_, combined, kspace_wav_);

            if (input_in_kspace_)
            {
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
        GADGET_THROW("Errors in hoWavelet3DOperator<T>::mult_MH(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCPUOPERATOR hoWavelet3DOperator< std::complex<float> >;
template class EXPORTCPUOPERATOR hoWavelet3DOperator< std::complex<double> >;

}
