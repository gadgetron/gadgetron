/** \file       gtPlusOperator.h
    \brief      Base class for gtPlus operators
    \author     Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"
#include "GadgetronTimer.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusIOAnalyze.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusOperator
{
public:

    typedef typename realType<T>::Type value_type;

    gtPlusOperator();
    virtual ~gtPlusOperator();

    virtual void printInfo(std::ostream& os);

    // forward operator
    virtual bool forwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y) = 0;

    // adjoint operator
    virtual bool adjointOperator(const hoNDArray<T>& x, hoNDArray<T>& y) = 0;

    // adjoint - forward operator
    virtual bool adjointforwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y);

    // compute gradient
    virtual bool grad(const hoNDArray<T>& x, hoNDArray<T>& g) = 0;

    // compute cost value
    virtual bool obj(const hoNDArray<T>& x, T& obj) = 0;

    // perform the proximity operation
    virtual bool proximity(hoNDArray<T>& x, value_type thres);

    // indicate the operator is unitary or not
    // unitary operator, AA' = I
    virtual bool unitary() const = 0;

    // restore acquired kspace points to x
    virtual bool restoreAcquiredKSpace(const hoNDArray<T>& acquired, hoNDArray<T>& y);
    virtual bool restoreAcquiredKSpace(hoNDArray<T>& y);

    // set the acquired kspace, unacquired points are set to be zero
    virtual bool setAcquiredPoints(boost::shared_ptr< hoNDArray<T> >& kspace);

    // set the coil sensivity map
    virtual bool setCoilSenMap(boost::shared_ptr< hoNDArray<T> >& senMap);

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer1_;
    Gadgetron::GadgetronTimer gt_timer2_;
    Gadgetron::GadgetronTimer gt_timer3_;

    bool performTiming_;

    // exporter
    Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

    // debug folder
    std::string debugFolder_;

    // util
    gtPlusISMRMRDReconUtil<T> gtPlus_util_;
    gtPlusISMRMRDReconUtilComplex<T> gtPlus_util_complex_;

public:

    // acquired kspace (unacquired points are zeros)
    boost::shared_ptr< hoNDArray<T> > acquired_points_;
    // acquired point indicator array, acquired points as 1, otherwise, 0
    hoNDArray<T> acquired_points_indicator_;
    // unacquired point indicator array
    hoNDArray<T> unacquired_points_indicator_;

    // coil map
    boost::shared_ptr< hoNDArray<T> > coil_senMap_;

    // helper memory
    hoNDArray<T> kspace_;
    hoNDArray<T> complexIm_;
    hoNDArray<T> res_after_apply_kernel_;
    hoNDArray<T> res_after_apply_kernel_sum_over_;

    hoNDArray<T> kspace_Managed_;
    hoNDArray<T> complexIm_Managed_;
    hoNDArray<T> res_after_apply_kernel_Managed_;
    hoNDArray<T> res_after_apply_kernel_sum_over_Managed_;

    bool performSumOverSrcChannel(const hoNDArray<T>& x, hoNDArray<T>& r);
};

template <typename T> 
gtPlusOperator<T>::gtPlusOperator() : performTiming_(false)
{
    gt_timer1_.set_timing_in_destruction(false);
    gt_timer2_.set_timing_in_destruction(false);
    gt_timer3_.set_timing_in_destruction(false);
}

template <typename T> 
gtPlusOperator<T>::~gtPlusOperator()
{
}

template <typename T> 
void gtPlusOperator<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD operator ------------------" << endl;
    os << "Operator for gtPlus ISMRMRD package" << endl;
    os << "----------------------------------------------------------------------" << endl;
}

template <typename T> 
bool gtPlusOperator<T>::adjointforwardOperator(const hoNDArray<T>& x, hoNDArray<T>& y)
{
    hoNDArray<T> a(x);
    GADGET_CHECK_RETURN_FALSE(this->forwardOperator(x, a));
    GADGET_CHECK_RETURN_FALSE(this->adjointOperator(a, y));
    return true;
}

template <typename T> 
bool gtPlusOperator<T>::restoreAcquiredKSpace(hoNDArray<T>& y)
{
    return this->restoreAcquiredKSpace(*acquired_points_, y);
}

template <typename T> 
bool gtPlusOperator<T>::restoreAcquiredKSpace(const hoNDArray<T>& acquired, hoNDArray<T>& y)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(acquired.get_number_of_elements()==y.get_number_of_elements());

        size_t N = acquired.get_number_of_elements();

        const T* pA = acquired.get_data_ptr();
        T* pY = y.get_data_ptr();

        int n;
        #pragma omp parallel for default(none) private(n) shared(N, pA, pY)
        for ( n=0; n<(int)N; n++ )
        {
            if ( std::abs(pA[n]) > 0 )
            {
                pY[n] = pA[n];
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in gtPlusOperator<T>::restoreAcquiredKSpace(const hoNDArray<T>& acquired, hoNDArray<T>& y) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusOperator<T>::
setAcquiredPoints(boost::shared_ptr< hoNDArray<T> >& kspace)
{
    try
    {
        acquired_points_ = kspace;

        acquired_points_indicator_.create(kspace->get_dimensions());
        Gadgetron::clear(acquired_points_indicator_);

        unacquired_points_indicator_.create(kspace->get_dimensions());
        Gadgetron::clear(unacquired_points_indicator_);

        size_t N = kspace->get_number_of_elements();

        long long ii;

        #pragma omp parallel for default(shared) private(ii) shared(N, kspace)
        for ( ii=0; ii<(long long)N; ii++ )
        {
            if ( std::abs( (*kspace)(ii) ) < DBL_EPSILON )
            {
                unacquired_points_indicator_(ii) = T(1.0);
            }
            else
            {
                acquired_points_indicator_(ii) = T(1.0);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusOperator<T>::setAcquiredPoints(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusOperator<T>::
setCoilSenMap(boost::shared_ptr< hoNDArray<T> >& senMap)
{
    try
    {
        coil_senMap_ = senMap;
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusOperator<T>::setCoilSenMap(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusOperator<T>::
proximity(hoNDArray<T>& /*x*/, value_type /*thres*/)
{
    return true;
}

template<typename T>
bool gtPlusOperator<T>::
performSumOverSrcChannel(const hoNDArray<T>& x, hoNDArray<T>& r)
{
    try
    {
        boost::shared_ptr< std::vector<size_t> > dim = x.get_dimensions();
        size_t NDim = dim->size();

        if (NDim < 2) return true;

        std::vector<size_t> dimR(NDim - 1);
        std::vector<size_t> dimRInternal = *dim;
        dimRInternal[NDim - 2] = 1;

        size_t d;
        for (d = 0; d<NDim - 2; d++)
        {
            dimR[d] = (*dim)[d];
        }
        dimR[NDim - 2] = (*dim)[NDim - 1];

        if (!r.dimensions_equal(&dimR))
        {
            r.create(&dimR);
        }

        if (x.get_size(NDim - 2) <= 1)
        {
            memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
            return true;
        }

        hoNDArray<T> rSum(dimRInternal, r.begin());

        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(x, rSum, NDim - 2));
    }
    catch (...)
    {
        GERROR_STREAM("Errors in performSumOverSrcChannel(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
        return false;
    }
    return true;
}

}}
