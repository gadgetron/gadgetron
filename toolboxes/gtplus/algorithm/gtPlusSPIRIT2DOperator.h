/** \file       gtPlusSPIRIT2DOperator.h
    \brief      Base class for gtPlus 2D operators
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSPIRITOperator.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusSPIRIT2DOperator : public gtPlusSPIRITOperator<T>
{
public:

    typedef gtPlusSPIRITOperator<T> BaseClass;

    gtPlusSPIRIT2DOperator() : BaseClass() {}
    virtual ~gtPlusSPIRIT2DOperator() {}

    virtual void printInfo(std::ostream& os);

    // convert to image domain or back to kspace
    virtual bool convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im);
    virtual bool convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x);

    // forward
    virtual bool forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    // adjoint operator
    virtual bool adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    using BaseClass::use_symmetric_spirit_;
    using BaseClass::use_non_centered_fft_;
    using BaseClass::calib_use_gpu_;

public:

    // [RO E1 srcCHA dstCHA]
    using BaseClass::forward_kernel_;
    using BaseClass::adjoint_kernel_;
    using BaseClass::adjoint_forward_kernel_;
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
void gtPlusSPIRIT2DOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD SPIRIT 2D operator ------------------" << endl;
    os << "Implementation of SPIRIT 2D operator for ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
inline bool gtPlusSPIRIT2DOperator<T>::convertToImage(const hoNDArray<T>& x, hoNDArray<T>& im)
{
    if ( this->use_non_centered_fft_ )
    {
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2(x, im));
    }
    else
    {
        if ( !complexIm_Managed_.dimensions_equal(&x) )
        {
            complexIm_Managed_.create(x.get_dimensions());
        }

        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(x, im, complexIm_Managed_));
    }

    return true;
}

template <typename T> 
inline bool gtPlusSPIRIT2DOperator<T>::convertToKSpace(const hoNDArray<T>& im, hoNDArray<T>& x)
{
    if ( this->use_non_centered_fft_ )
    {
        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2(im, x));
    }
    else
    {
        if ( !kspace_Managed_.dimensions_equal(&im) )
        {
            kspace_Managed_.create(im.get_dimensions());
        }

        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(im, x, kspace_Managed_));
    }

    return true;
}

template <typename T> 
inline bool gtPlusSPIRIT2DOperator<T>::forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    try
    {
        Gadgetron::multiply(unacquired_points_indicator_, x, y);

        // x to image domain
        this->convertToImage(y, complexIm_);

        size_t ro = x.get_size(0);
        size_t e1 = x.get_size(1);
        size_t CHA = x.get_size(2);

        if ( res_after_apply_kernel_sum_over_.get_number_of_elements() < ro*e1*CHA )
        {
            res_after_apply_kernel_sum_over_.create(ro, e1, CHA);
        }

        hoNDArray<T>* kerArray;
        if ( use_symmetric_spirit_ )
        {
            kerArray = this->adjoint_forward_kernel_.get();
        }
        else
        {
            kerArray = this->forward_kernel_.get();
        }

        Gadgetron::imageDomainUnwrapping2D(complexIm_, *kerArray, res_after_apply_kernel_sum_over_, y);

        /*
        //long long dCha;

        //#pragma omp parallel
        {
            //#ifdef WIN32
            //    int tid = omp_get_thread_num();
            //    DWORD_PTR mask = (1 << tid);
            //    // GADGET_MSG("thread id : " << tid << " - mask : " << mask);
            //    SetThreadAffinityMask( GetCurrentThread(), mask );
            //#endif // WIN32

            //#pragma omp for
            if ( typeid(T)==typeid(GT_Complex8) )
            {
                for ( dCha=0; dCha<CHA; dCha++ )
                {
                    vcMul(ro*e1*CHA, reinterpret_cast<MKL_Complex8*>(pIm), 
                        reinterpret_cast<MKL_Complex8*>(ker+dCha*ro*e1*CHA), 
                        reinterpret_cast<MKL_Complex8*>(ptt));

                    memcpy(pY+dCha*ro*e1, ptt, sizeof(T)*ro*e1);
                    for ( size_t sCha=1; sCha<CHA; sCha++ )
                    {
                        vcAdd(ro*e1, reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1), 
                            reinterpret_cast<MKL_Complex8*>(ptt+sCha*ro*e1), 
                            reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1));
                    }
                }
            }
            else if ( typeid(T)==typeid(GT_Complex16) )
            {
                for ( dCha=0; dCha<CHA; dCha++ )
                {
                    vzMul(ro*e1*CHA, reinterpret_cast<MKL_Complex16*>(pIm), 
                        reinterpret_cast<MKL_Complex16*>(ker+dCha*ro*e1*CHA), 
                        reinterpret_cast<MKL_Complex16*>(ptt));

                    memcpy(pY+dCha*ro*e1, ptt, sizeof(T)*ro*e1);
                    for ( size_t sCha=1; sCha<CHA; sCha++ )
                    {
                        vzAdd(ro*e1, reinterpret_cast<MKL_Complex16*>(pY+dCha*ro*e1), 
                            reinterpret_cast<MKL_Complex16*>(ptt+sCha*ro*e1), 
                            reinterpret_cast<MKL_Complex16*>(pY+dCha*ro*e1));
                    }
                }
            }
        }*/

        this->convertToKSpace(y, res_after_apply_kernel_sum_over_);

        // apply Dc
        if ( use_symmetric_spirit_ )
        {
            Gadgetron::multiply(unacquired_points_indicator_, res_after_apply_kernel_sum_over_, y);
        }
        else
        {
            memcpy(y.begin(), res_after_apply_kernel_sum_over_.begin(), sizeof(T)*ro*e1*CHA);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRIT2DOperator<T>::forwardOperator(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
inline bool gtPlusSPIRIT2DOperator<T>::adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
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
            //Gadgetron::multipleMultiply(complexIm_, *adjoint_kernel_, res_after_apply_kernel_);
            //Gadgetron::sumOverSecondLastDimension(res_after_apply_kernel_, res_after_apply_kernel_sum_over_);

            size_t ro = x.get_size(0);
            size_t e1 = x.get_size(1);
            size_t CHA = x.get_size(2);

            if ( res_after_apply_kernel_sum_over_.get_number_of_elements() < ro*e1*CHA )
            {
                res_after_apply_kernel_sum_over_.create(ro, e1, CHA);
            }

            Gadgetron::imageDomainUnwrapping2D(complexIm_, *adjoint_kernel_, res_after_apply_kernel_sum_over_, y);

            //long long dCha;

            ////#pragma omp parallel default(shared)
            //{
            //    //#ifdef WIN32
            //    //    int tid = omp_get_thread_num();
            //    //    DWORD_PTR mask = (1 << tid);
            //    //    // GADGET_MSG("thread id : " << tid << " - mask : " << mask);
            //    //    SetThreadAffinityMask( GetCurrentThread(), mask );
            //    //#endif // WIN32

            //    //#pragma omp for

            //    if ( typeid(T)==typeid(GT_Complex8) )
            //    {
            //        for ( dCha=0; dCha<CHA; dCha++ )
            //        {
            //            vcMul(ro*e1*CHA, reinterpret_cast<MKL_Complex8*>(pIm), 
            //                reinterpret_cast<MKL_Complex8*>(ker+dCha*ro*e1*CHA), 
            //                reinterpret_cast<MKL_Complex8*>(ptt));

            //            memcpy(pY+dCha*ro*e1, ptt, sizeof(T)*ro*e1);
            //            for ( size_t sCha=1; sCha<CHA; sCha++ )
            //            {
            //                vcAdd(ro*e1, reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1), 
            //                    reinterpret_cast<MKL_Complex8*>(ptt+sCha*ro*e1), 
            //                    reinterpret_cast<MKL_Complex8*>(pY+dCha*ro*e1));
            //            }
            //        }
            //    }

            //}

            // go back to kspace 
            this->convertToKSpace(y, res_after_apply_kernel_sum_over_);

            // apply Dc
            Gadgetron::multiply(unacquired_points_indicator_, res_after_apply_kernel_sum_over_, y);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusSPIRITOperator<T>::adjointOperator(...) ... ");
        return false;
    }
    return true;
}

}}
