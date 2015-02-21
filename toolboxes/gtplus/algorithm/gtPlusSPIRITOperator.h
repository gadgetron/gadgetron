/** \file       gtPlusSPIRITOperator.h
    \brief      Implement SPIRIT operator functinalities
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSPIRIT.h"
#include "gtPlusOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusSPIRITOperator : public gtPlusSPIRIT<T>, public gtPlusOperator<T>
{
public:

    typedef gtPlusOperator<T> BaseClass;

    gtPlusSPIRITOperator() : use_symmetric_spirit_(false), use_non_centered_fft_(false), BaseClass() {}
    virtual ~gtPlusSPIRITOperator() {}

    virtual void printInfo(std::ostream& os);

    // set forward kernel, compute the adjoint and adjoint_forward kernel
    bool setForwardKernel(boost::shared_ptr< hoNDArray<T> >& forward_kernel, bool computeAdjForwardKernel=true);
    bool setAdjointForwardKernel(boost::shared_ptr< hoNDArray<T> >& adjoint_forward_kernel);

    hoNDArray<T>* getAdjointKernel();
    hoNDArray<T>* getAdjointForwardKernel();

    // apply Dc(G-I)'(G-I)Dc' to x
    virtual bool forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    // adjoint operator
    virtual bool adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    // compute right hand side
    // b = -Dc(G-I)'(G-I)D'x
    virtual bool computeRighHandSide(const hoNDArray<T>& x, hoNDArray<T>& b);

    // compute gradient of ||(G-I)(Dc'x+D'y)||2
    virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g);

    // compute cost value of L2 norm ||(G-I)(Dc'x+D'y)||2
    virtual bool obj(const hoNDArray<T>& x, T& obj);

    // indicate the operator is unitary or not
    // unitary operator, AA' = I
    virtual bool unitary() const { return false; }

    // convert to image domain or back to kspace
    virtual bool convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im) = 0;
    virtual bool convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x) = 0;

    // whether to use symmetric spirit equation
    // symmetric equation: A = Dc(G-I)'(G-I)Dc'
    // non-symmetric equation: A = (G-I)Dc'
    bool use_symmetric_spirit_;

    // if true, use the fft. not fftc
    bool use_non_centered_fft_;

    using gtPlusSPIRIT<T>::calib_use_gpu_;
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

    // G-I, [... srcCHA dstCHA]
    boost::shared_ptr< hoNDArray<T> > forward_kernel_;
    // (G-I)', [... dstCHA srcCHA]
    boost::shared_ptr< hoNDArray<T> > adjoint_kernel_;
    // (G-I)'(G-I), [... dstCHA dstCHA]
    boost::shared_ptr< hoNDArray<T> > adjoint_forward_kernel_;

    using BaseClass::acquired_points_;
    using BaseClass::acquired_points_indicator_;
    using BaseClass::unacquired_points_indicator_;

    // helper memory
    using BaseClass::kspace_;
    using BaseClass::complexIm_;
    using BaseClass::res_after_apply_kernel_;
    using BaseClass::res_after_apply_kernel_sum_over_;

    using BaseClass::kspace_Managed_;
    using BaseClass::complexIm_Managed_;
    using BaseClass::res_after_apply_kernel_Managed_;
    using BaseClass::res_after_apply_kernel_sum_over_Managed_;
};

template <typename T> 
void gtPlusSPIRITOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT operator ------------------" << endl;
    os << "Implementation of SPIRIT operator for ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
inline hoNDArray<T>* gtPlusSPIRITOperator<T>::getAdjointKernel()
{
    return adjoint_kernel_.get();
}

template <typename T> 
inline hoNDArray<T>* gtPlusSPIRITOperator<T>::getAdjointForwardKernel()
{
    return adjoint_forward_kernel_.get();
}

template <typename T> 
bool gtPlusSPIRITOperator<T>::
setForwardKernel(boost::shared_ptr< hoNDArray<T> >& forward_kernel, bool computeAdjForwardKernel)
{
    try
    {
        forward_kernel_ = forward_kernel;

        adjoint_kernel_ = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>());
        GADGET_CHECK_RETURN_FALSE(this->imageDomainAdjointKernel(*forward_kernel_, *adjoint_kernel_));

        if ( computeAdjForwardKernel || use_symmetric_spirit_ )
        {
            adjoint_forward_kernel_ = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>());
            GADGET_CHECK_RETURN_FALSE(this->AdjointForwardKernel(*adjoint_kernel_, *forward_kernel_, *adjoint_forward_kernel_));
        }

        // allocate the helper memory
        boost::shared_ptr< std::vector<size_t> > dims = forward_kernel->get_dimensions();
        size_t NDim = dims->size();

        std::vector<size_t> dimSrc(NDim-1), dimDst(NDim-1);
        size_t ii;
        for ( ii=0; ii<NDim-2; ii++ )
        {
            dimSrc[ii] = (*dims)[ii];
            dimDst[ii] = (*dims)[ii];
        }

        dimSrc[NDim-2] = (*dims)[NDim-2];
        dimDst[NDim-2] = (*dims)[NDim-1];

        kspace_.create(dimSrc);
        complexIm_.create(dimSrc);
        res_after_apply_kernel_.create(dims);
        res_after_apply_kernel_sum_over_.create(dimDst);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRITOperator<T>::setForwardKernel(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITOperator<T>::
setAdjointForwardKernel(boost::shared_ptr< hoNDArray<T> >& adjoint_forward_kernel)
{
    try
    {
        adjoint_forward_kernel_ = adjoint_forward_kernel;

        // allocate the helper memory
        boost::shared_ptr< std::vector<size_t> > dims = adjoint_forward_kernel_->get_dimensions();
        size_t NDim = dims->size();

        std::vector<size_t> dimSrc(NDim-1), dimDst(NDim-1);
        size_t ii;
        for ( ii=0; ii<NDim-2; ii++ )
        {
            dimSrc[ii] = (*dims)[ii];
            dimDst[ii] = (*dims)[ii];
        }

        dimSrc[NDim-2] = (*dims)[NDim-2];
        dimDst[NDim-2] = (*dims)[NDim-1];

        kspace_.create(dimSrc);
        complexIm_.create(dimSrc);
        res_after_apply_kernel_.create(dims);
        res_after_apply_kernel_sum_over_.create(dimDst);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRITOperator<T>::setAdjointForwardKernel(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITOperator<T>::forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        // Dc(G-I)'(G-I)Dc'x

        Gadgetron::multiply(unacquired_points_indicator_, x, y);

        // x to image domain
        this->convertToImage(y, complexIm_);

        // apply kernel and sum
        if ( use_symmetric_spirit_ )
        {
            Gadgetron::multiply(*adjoint_forward_kernel_, complexIm_, res_after_apply_kernel_);
        }
        else
        {
            Gadgetron::multiply(*forward_kernel_, complexIm_, res_after_apply_kernel_);
        }

        this->performSumOverSrcChannel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_);

        // go back to kspace 
        this->convertToKSpace(res_after_apply_kernel_sum_over_, y);

        // apply Dc
        if ( use_symmetric_spirit_ )
        {
            Gadgetron::multiply(unacquired_points_indicator_, y, y);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRITOperator<T>::forwardOperator(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITOperator<T>::adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        if ( use_symmetric_spirit_ )
        {
            // Dc(G-I)'(G-I)Dc' is symmetric
            GADGET_CHECK_RETURN_FALSE(this->forwardOperator(x, y));
        }
        else
        {
            // Dc(G-I)'x

            // x to image domain
            this->convertToImage(x, complexIm_);

            // apply kernel and sum
            Gadgetron::multiply(*adjoint_kernel_, complexIm_, res_after_apply_kernel_);
            this->performSumOverSrcChannel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_);

            // go back to kspace 
            this->convertToKSpace(res_after_apply_kernel_sum_over_, y);

            // apply Dc
            Gadgetron::multiply(unacquired_points_indicator_, y, y);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRITOperator<T>::adjointOperator(...) ... ");
        return false;
    }
    return true;
}

template <typename T> 
bool gtPlusSPIRITOperator<T>::computeRighHandSide(const hoNDArray<T>& x, hoNDArray<T>& b)
{
    try
    {
        // symmetric: -Dc(G-I)'(G-I)D'x
        // non-symmetric: -(G-I)D'x

        // D'x, need to do nothing, acquired points are already in place

        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, complexIm_));

        // apply kernel and sum
        if ( use_symmetric_spirit_ )
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(*adjoint_forward_kernel_, complexIm_, res_after_apply_kernel_));
        }
        else
        {
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(*forward_kernel_, complexIm_, res_after_apply_kernel_));
        }

        GADGET_CHECK_RETURN_FALSE(this->performSumOverSrcChannel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_));

        // go back to kspace 
        GADGET_CHECK_RETURN_FALSE(this->convertToKSpace(res_after_apply_kernel_sum_over_, b));

        // apply Dc
        if ( use_symmetric_spirit_ )
        {
            Gadgetron::multiply(unacquired_points_indicator_, b, b);
        }

        // multiply by -1
        Gadgetron::scal( (typename realType<T>::Type)(-1.0), b);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRITOperator<T>::computeRighHandSide(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITOperator<T>::grad(const hoNDArray<T>& x, hoNDArray<T>& g)
{
    try
    {
        // gradient of L2 norm is
        // 2*Dc*(G-I)'(G-I)(D'y+Dc'x)

        // D'y+Dc'x
        Gadgetron::multiply(unacquired_points_indicator_, x, kspace_);
        Gadgetron::add(*acquired_points_, kspace_, kspace_);

        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(kspace_, complexIm_));

        // apply kernel and sum
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(*adjoint_forward_kernel_, complexIm_, res_after_apply_kernel_));
        GADGET_CHECK_RETURN_FALSE(this->performSumOverSrcChannel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_));

        // go back to kspace 
        GADGET_CHECK_RETURN_FALSE(this->convertToKSpace(res_after_apply_kernel_sum_over_, g));

        // apply Dc
        Gadgetron::multiply(unacquired_points_indicator_, g, g);

        // multiply by 2
        Gadgetron::scal( (typename realType<T>::Type)(2.0), g);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRITOperator<T>::grad(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITOperator<T>::obj(const hoNDArray<T>& x, T& obj)
{
    try
    {
        // L2 norm
        // ||(G-I)(D'y+Dc'x)||2

        // D'y+Dc'x
        Gadgetron::multiply(unacquired_points_indicator_, x, kspace_);
        Gadgetron::add(*acquired_points_, kspace_, kspace_);

        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(kspace_, complexIm_));

        // apply kernel and sum
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(*forward_kernel_, complexIm_, res_after_apply_kernel_));
        GADGET_CHECK_RETURN_FALSE(this->performSumOverSrcChannel(res_after_apply_kernel_, res_after_apply_kernel_sum_over_));

        // L2 norm
        Gadgetron::dotc(res_after_apply_kernel_sum_over_, res_after_apply_kernel_sum_over_, obj);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusSPIRITOperator<T>::grad(...) ... ");
        return false;
    }

    return true;
}

}}
