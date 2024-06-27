
#include "hoSPIRITOperator.h"
#include "mri_core_spirit.h"

namespace Gadgetron 
{

template <typename T> 
hoSPIRITOperator<T>::hoSPIRITOperator(std::vector<size_t> *dims) : use_non_centered_fft_(false), no_null_space_(false), BaseClass(dims)
{
}

template <typename T> 
hoSPIRITOperator<T>::~hoSPIRITOperator()
{
}

template <typename T> 
void hoSPIRITOperator<T>::restore_acquired_kspace(const ARRAY_TYPE& acquired, ARRAY_TYPE& y)
{
    try
    {
        GADGET_CHECK_THROW(acquired.get_number_of_elements() == y.get_number_of_elements());

        size_t N = acquired.get_number_of_elements();

        const T* pA = acquired.get_data_ptr();
        T* pY = y.get_data_ptr();

        int n(0);
#pragma omp parallel for default(none) private(n) shared(N, pA, pY)
        for (n = 0; n<(int)N; n++)
        {
            if (std::abs(pA[n]) > 0)
            {
                pY[n] = pA[n];
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in hoSPIRITOperator<T>::restore_acquired_kspace(...) ... ");
    }
}

template <typename T> 
void hoSPIRITOperator<T>::restore_acquired_kspace(ARRAY_TYPE& y)
{
    this->restore_acquired_kspace(acquired_points_, y);
}

template <typename T> 
void hoSPIRITOperator<T>::set_acquired_points(ARRAY_TYPE& kspace)
{
    try
    {
        std::vector<size_t> dim = kspace.dimensions();
        acquired_points_.create(dim, kspace.begin());

        acquired_points_indicator_.create(kspace.dimensions());
        Gadgetron::clear(acquired_points_indicator_);

        unacquired_points_indicator_.create(kspace.dimensions());
        Gadgetron::clear(unacquired_points_indicator_);

        size_t N = kspace.get_number_of_elements();

        long long ii(0);

#pragma omp parallel for default(shared) private(ii) shared(N, kspace)
        for (ii = 0; ii<(long long)N; ii++)
        {
            if (std::abs(kspace(ii)) < DBL_EPSILON)
            {
                unacquired_points_indicator_(ii) = T(1.0);
            }
            else
            {
                acquired_points_indicator_(ii) = T(1.0);
            }
        }

        // allocate the helper memory
        kspace_.create(kspace.dimensions());
        complexIm_.create(kspace.dimensions());
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRITOperator<T>::set_acquired_points(...) ... ");
    }
}

template <typename T> 
void hoSPIRITOperator<T>::set_coil_sen_map(ARRAY_TYPE& senMap)
{
    coil_senMap_ = senMap;
}

template <typename T> 
void hoSPIRITOperator<T>::set_forward_kernel(ARRAY_TYPE& forward_kernel, bool compute_adjoint_forward_kernel)
{
    try
    {
        std::vector<size_t> dim;
        forward_kernel.get_dimensions(dim);
        forward_kernel_.create(dim, forward_kernel.begin());

        GADGET_CATCH_THROW(Gadgetron::spirit_image_domain_adjoint_kernel(forward_kernel_, adjoint_kernel_));

        if (compute_adjoint_forward_kernel)
        {
            GADGET_CATCH_THROW(Gadgetron::spirit_adjoint_forward_kernel(adjoint_kernel_, forward_kernel_, adjoint_forward_kernel_));
        }

        // allocate the helper memory
        std::vector<size_t> dims;
        forward_kernel.get_dimensions(dims);
        size_t NDim = dims.size();

        std::vector<size_t> dimSrc(NDim - 1), dimDst(NDim - 1);
        size_t ii;
        for (ii = 0; ii < NDim - 2; ii++)
        {
            dimSrc[ii] = dims[ii];
            dimDst[ii] = dims[ii];
        }

        dimSrc[NDim - 2] = dims[NDim - 2];
        dimDst[NDim - 2] = dims[NDim - 1];

        res_after_apply_kernel_.create(dims);
        res_after_apply_kernel_sum_over_.create(dimDst);
        kspace_dst_.create(dimDst);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRITOperator<T>::set_forward_kernel(...) ... ");
    }
}

template<typename T>
void hoSPIRITOperator<T>::sum_over_src_channel(const ARRAY_TYPE& x, ARRAY_TYPE& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();
        size_t NDim = dim->size();

        if (NDim < 2) return;

        std::vector<size_t> dimR(NDim - 1);
        std::vector<size_t> dimRInternal = *dim;
        dimRInternal[NDim - 2] = 1;

        size_t d;
        for (d = 0; d<NDim - 2; d++)
        {
            dimR[d] = (*dim)[d];
        }
        dimR[NDim - 2] = (*dim)[NDim - 1];

        if (!r.dimensions_equal(dimR))
        {
            r.create(dimR);
        }

        if (x.get_size(NDim - 2) <= 1)
        {
            memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
            return;
        }

        hoNDArray<T> rSum(dimRInternal, r.begin());

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(x, rSum, NDim - 2));
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRITOperator<T>::sum_over_src_channel(const ARRAY_TYPE& x, ARRAY_TYPE& r) ... ");
    }
}

template <typename T>
void hoSPIRITOperator<T>::mult_M(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        if (accumulate) 
        {
            kspace_dst_ = *y;
        }

        if(no_null_space_)
        {
            // (G-I)x
            this->convert_to_image(*x, complexIm_);
        }
        else
        {
            // (G-I)Dc'x
            Gadgetron::multiply(unacquired_points_indicator_, *x, *y);

            // x to image domain
            this->convert_to_image(*y, complexIm_);
        }

        // apply kernel and sum
        Gadgetron::multiply(forward_kernel_, complexIm_, res_after_apply_kernel_);

        this->sum_over_src_channel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_);

        // go back to kspace 
        this->convert_to_kspace(res_after_apply_kernel_sum_over_, *y);

        if(accumulate)
        {
            Gadgetron::add(kspace_dst_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRITOperator<T>::mult_M(...) ... ");
    }
}

template <typename T>
void hoSPIRITOperator<T>::mult_MH(ARRAY_TYPE* x, ARRAY_TYPE* y, bool accumulate)
{
    try
    {
        // Dc(G-I)'x or if no_null_space_ == true, (G-I)'x

        if(accumulate)
        {
            kspace_ = *y;
        }

        // x to image domain
        this->convert_to_image(*x, complexIm_);

        // apply kernel and sum
        Gadgetron::multiply(adjoint_kernel_, complexIm_, res_after_apply_kernel_);
        this->sum_over_src_channel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_);

        // go back to kspace 
        this->convert_to_kspace(res_after_apply_kernel_sum_over_, *y);

        if (!no_null_space_)
        {
            // apply Dc
            Gadgetron::multiply(unacquired_points_indicator_, *y, *y);
        }

        if (accumulate)
        {
            Gadgetron::add(kspace_, *y, *y);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRITOperator<T>::mult_MH(...) ... ");
    }
}

template <typename T>
void hoSPIRITOperator<T>::compute_righ_hand_side(const ARRAY_TYPE& x, ARRAY_TYPE& b)
{
    try
    {
        if (no_null_space_)
        {
            b.create(x.dimensions());
            Gadgetron::clear(b);
        }
        else
        {
            // non-symmetric rhs: -(G-I)D'x

            // need to be done for D'x, acquired points are already in place

            // x to image domain
            this->convert_to_image(x, complexIm_);

            // apply kernel and sum
            GADGET_CATCH_THROW(Gadgetron::multiply(forward_kernel_, complexIm_, res_after_apply_kernel_));

            GADGET_CATCH_THROW(this->sum_over_src_channel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_));

            // go back to kspace 
            this->convert_to_kspace(res_after_apply_kernel_sum_over_, b);

            // multiply by -1
            Gadgetron::scal((typename realType<T>::Type)(-1.0), b);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRITOperator<T>::compute_righ_hand_side(...) ... ");
    }
}

template <typename T>
void hoSPIRITOperator<T>::gradient(ARRAY_TYPE* x, ARRAY_TYPE* g, bool accumulate)
{
    try
    {
        if (accumulate)
        {
            kspace_ = *g;
        }

        if (no_null_space_)
        {
            // gradient of L2 norm is 2*Dc*(G-I)'(G-I)x
            this->convert_to_image(*x, complexIm_);
        }
        else
        {
            // gradient of L2 norm is 2*Dc*(G-I)'(G-I)(D'y+Dc'x)
            Gadgetron::multiply(unacquired_points_indicator_, *x, kspace_);
            Gadgetron::add(acquired_points_, kspace_, kspace_);

            // x to image domain
            this->convert_to_image(kspace_, complexIm_);
        }

        // apply kernel and sum
        Gadgetron::multiply(adjoint_forward_kernel_, complexIm_, res_after_apply_kernel_);
        this->sum_over_src_channel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_);

        // go back to kspace 
        this->convert_to_kspace(res_after_apply_kernel_sum_over_, *g);

        // apply Dc
        Gadgetron::multiply(unacquired_points_indicator_, *g, *g);

        // multiply by 2
        Gadgetron::scal((typename realType<T>::Type)(2.0), *g);

        if (accumulate)
        {
            Gadgetron:add(kspace_, *g, *g);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRITOperator<T>::gradient(...) ... ");
    }
}

template <typename T>
typename hoSPIRITOperator<T>::REAL hoSPIRITOperator<T>::magnitude(ARRAY_TYPE* x)
{
    try
    {
        if (no_null_space_)
        {
            // L2 norm of ||(G-I)x||2
            this->convert_to_image(*x, complexIm_);
        }
        else
        {
            // L2 norm of ||(G-I)(D'y+Dc'x)||2
            // D'y+Dc'x
            Gadgetron::multiply(unacquired_points_indicator_, *x, kspace_);
            Gadgetron::add(acquired_points_, kspace_, kspace_);

            // x to image domain
            this->convert_to_image(kspace_, complexIm_);
        }

        // apply kernel and sum
        Gadgetron::multiply(forward_kernel_, complexIm_, res_after_apply_kernel_);
        this->sum_over_src_channel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_);

        // L2 norm
        T obj = Gadgetron::dot(res_after_apply_kernel_sum_over_, res_after_apply_kernel_sum_over_, true);

        return std::abs(obj);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoSPIRITOperator<T>::magnitude(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCPUOPERATOR hoSPIRITOperator< std::complex<float> >;
template class EXPORTCPUOPERATOR hoSPIRITOperator< std::complex<double> >;

}
