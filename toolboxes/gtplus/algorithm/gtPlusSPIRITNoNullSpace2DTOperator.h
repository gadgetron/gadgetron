/** \file       gtPlusSPIRITNoNullSpace2DTOperator.h
    \brief      Implement SPIRIT 2DT operator without Null space
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSPIRITNoNullSpace2DOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusSPIRITNoNullSpace2DTOperator : public gtPlusSPIRITNoNullSpace2DOperator<T>
{
public:

    typedef gtPlusSPIRITNoNullSpace2DOperator<T> BaseClass;

    gtPlusSPIRITNoNullSpace2DTOperator() : BaseClass() {}
    virtual ~gtPlusSPIRITNoNullSpace2DTOperator() {}

    virtual void printInfo(std::ostream& os);

    // set forward kernel, compute the adjoint and adjoint_forward kernel
    bool setForwardKernel(boost::shared_ptr< hoNDArray<T> >& forward_kernel, bool computeAdjForwardKernel=true);
    bool setAdjointForwardKernel(boost::shared_ptr< hoNDArray<T> >& adjoint_forward_kernel);
    // set the acquired kspace, unacquired points are set to be zero
    bool setAcquiredPoints(boost::shared_ptr< hoNDArray<T> >& kspace);

    // compute gradient of ||(G-I)x||2
    virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g);

    // compute cost value of L2 norm ||(G-I)x||2
    virtual bool obj(const hoNDArray<T>& x, T& obj);

protected:

};

template <typename T> 
void gtPlusSPIRITNoNullSpace2DTOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT 2DT operator without null space constraint ------------------" << endl;
    os << "Implementation of SPIRIT 2DT operator for ISMRMRD package" << endl;
    os << "--------------------------------------------------------------------------------------------------" << endl;
}

template <typename T> 
bool gtPlusSPIRITNoNullSpace2DTOperator<T>::
setForwardKernel(boost::shared_ptr< hoNDArray<T> >& forward_kernel, bool computeAdjForwardKernel)
{
    try
    {
        this->forward_kernel_ = forward_kernel;

        size_t RO = this->forward_kernel_->get_size(0);
        size_t E1 = this->forward_kernel_->get_size(1);
        size_t srcCHA = this->forward_kernel_->get_size(2);
        size_t dstCHA = this->forward_kernel_->get_size(3);
        size_t N = this->forward_kernel_->get_size(4);

        this->adjoint_kernel_ = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>(RO, E1, dstCHA, srcCHA, N));

        bool computeAdjointForwardKernel = (computeAdjForwardKernel || this->use_symmetric_spirit_);

        if ( computeAdjointForwardKernel )
        {
            this->adjoint_forward_kernel_ = boost::shared_ptr< hoNDArray<T> >(new hoNDArray<T>(RO, E1, dstCHA, dstCHA, N));
        }

        size_t n;
        for ( n=0; n<N; n++ )
        {
            hoNDArray<T> kerCurr(RO, E1, srcCHA, dstCHA, this->forward_kernel_->begin()+n*RO*E1*srcCHA*dstCHA);
            hoNDArray<T> adjKerCurr(RO, E1, dstCHA, srcCHA, this->adjoint_kernel_->begin()+n*RO*E1*dstCHA*srcCHA);

            GADGET_CHECK_RETURN_FALSE(this->imageDomainAdjointKernel(kerCurr, adjKerCurr));

            if ( computeAdjointForwardKernel )
            {
                hoNDArray<T> adjForwardKerCurr(RO, E1, dstCHA, dstCHA, this->adjoint_forward_kernel_->begin()+n*RO*E1*dstCHA*dstCHA);
                GADGET_CHECK_RETURN_FALSE(this->AdjointForwardKernel(adjKerCurr, kerCurr, adjForwardKerCurr));
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRITNoNullSpace2DTOperator<T>::setForwardKernel(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITNoNullSpace2DTOperator<T>::
setAdjointForwardKernel(boost::shared_ptr< hoNDArray<T> >& adjoint_forward_kernel)
{
    try
    {
        this->adjoint_forward_kernel_ = adjoint_forward_kernel;
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRITNoNullSpace2DTOperator<T>::setAdjointForwardKernel(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITNoNullSpace2DTOperator<T>::
setAcquiredPoints(boost::shared_ptr< hoNDArray<T> >& kspace)
{
    try
    {
        this->acquired_points_ = kspace;

        size_t RO = this->acquired_points_->get_size(0);
        size_t E1 = this->acquired_points_->get_size(1);
        size_t srcCHA = this->acquired_points_->get_size(2);
        size_t E2 = this->acquired_points_->get_size(3);

        this->acquired_points_indicator_.create(kspace->get_dimensions());
        Gadgetron::clear(this->acquired_points_indicator_);

        this->unacquired_points_indicator_.create(kspace->get_dimensions());
        Gadgetron::clear(this->unacquired_points_indicator_);

        size_t N = kspace->get_number_of_elements();

        long long ii;

        #ifdef GCC_OLD_FLAG
            #pragma omp parallel for default(none) private(ii) shared(N)
        #else
            #pragma omp parallel for default(none) private(ii) shared(N, kspace)
        #endif
        for ( ii=0; ii<(long long)N; ii++ )
        {
            if ( std::abs( (*kspace)(ii) ) < DBL_EPSILON )
            {
                this->unacquired_points_indicator_(ii) = 1.0;
            }
            else
            {
                this->acquired_points_indicator_(ii) = 1.0;
            }
        }

        // allocate the helper memory
        this->kspace_.create(RO, E1, srcCHA, E2);
        this->complexIm_.create(RO, E1, srcCHA, E2);

        if ( this->forward_kernel_ )
        {
            size_t dstCHA = this->forward_kernel_->get_size(3);
            this->res_after_apply_kernel_.create(RO, E1, srcCHA, dstCHA);
            this->res_after_apply_kernel_sum_over_.create(RO, E1, dstCHA, E2);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRITNoNullSpace2DTOperator<T>::setAcquiredPoints(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITNoNullSpace2DTOperator<T>::grad(const hoNDArray<T>& x, hoNDArray<T>& g)
{
    try
    {
        // gradient of L2 norm is
        // 2*(G-I)'(G-I)x

        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, this->complexIm_));

        // apply kernel and sum
        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t CHA = x.get_size(2);
        size_t N = x.get_size(3);

        size_t dstCHA = this->adjoint_forward_kernel_->get_size(3);
        size_t kernelN = this->adjoint_forward_kernel_->get_size(4);

        this->res_after_apply_kernel_sum_over_.create(RO, E1, dstCHA, N);

        size_t n;
        for ( n=0; n<N; n++)
        {
            hoNDArray<T> currComplexIm(RO, E1, CHA, this->complexIm_.begin()+n*RO*E1*CHA);

            hoNDArray<T> curr_adjoint_forward_kernel;

            if ( n < kernelN )
            {
                curr_adjoint_forward_kernel.create(RO, E1, CHA, dstCHA, this->adjoint_forward_kernel_->begin()+n*RO*E1*CHA*dstCHA);
            }
            else
            {
                curr_adjoint_forward_kernel.create(RO, E1, CHA, dstCHA, this->adjoint_forward_kernel_->begin()+(kernelN-1)*RO*E1*CHA*dstCHA);
            }

            GADGET_CHECK_RETURN_FALSE(Gadgetron::multipleMultiply(currComplexIm, curr_adjoint_forward_kernel, this->res_after_apply_kernel_));

            hoNDArray<T> sumResCurr(RO, E1, dstCHA, this->res_after_apply_kernel_sum_over_.begin()+n*RO*E1*dstCHA);
            GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverSecondLastDimension(this->res_after_apply_kernel_, sumResCurr));
        }

        // go back to kspace 
        GADGET_CHECK_RETURN_FALSE(this->convertToKSpace(this->res_after_apply_kernel_sum_over_, g));

        // multiply by 2
        GADGET_CHECK_RETURN_FALSE(Gadgetron::scal(T(2.0), g));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRITNoNullSpace2DTOperator<T>::grad(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusSPIRITNoNullSpace2DTOperator<T>::obj(const hoNDArray<T>& x, T& obj)
{
    try
    {
        // L2 norm
        // ||(G-I)x||2

        // x to image domain
        GADGET_CHECK_RETURN_FALSE(this->convertToImage(x, this->complexIm_));

        // apply kernel and sum
        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);
        size_t CHA = x.get_size(2);
        size_t N = x.get_size(3);

        size_t dstCHA = this->forward_kernel_->get_size(3);
        size_t kernelN = this->forward_kernel_->get_size(4);

        this->res_after_apply_kernel_sum_over_.create(RO, E1, dstCHA, N);

        size_t n;
        for ( n=0; n<N; n++)
        {
            hoNDArray<T> currComplexIm(RO, E1, CHA, this->complexIm_.begin()+n*RO*E1*CHA);

            hoNDArray<T> curr_forward_kernel;

            if ( n < kernelN )
            {
                curr_forward_kernel.create(RO, E1, CHA, dstCHA, this->forward_kernel_->begin()+n*RO*E1*CHA*dstCHA);
            }
            else
            {
                curr_forward_kernel.create(RO, E1, CHA, dstCHA, this->forward_kernel_->begin()+(kernelN-1)*RO*E1*CHA*dstCHA);
            }

            GADGET_CHECK_RETURN_FALSE(Gadgetron::multipleMultiply(currComplexIm, curr_forward_kernel, this->res_after_apply_kernel_));

            hoNDArray<T> sumResCurr(RO, E1, dstCHA, this->res_after_apply_kernel_sum_over_.begin()+n*RO*E1*dstCHA);
            GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverSecondLastDimension(this->res_after_apply_kernel_, sumResCurr));
        }

        // L2 norm
        GADGET_CHECK_RETURN_FALSE(Gadgetron::dotc(this->res_after_apply_kernel_sum_over_, this->res_after_apply_kernel_sum_over_, obj));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRITNoNullSpace2DTOperator<T>::grad(...) ... ");
        return false;
    }

    return true;
}

}}
