/** \file       gtPlusSolver.h
    \brief      Define the base class for GtPlus solvers
    \author     Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"
#include "GadgetronTimer.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusIOAnalyze.h"
#include "gtPlusMemoryManager.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron { namespace gtPlus {

template <typename Array_Type_I, typename Array_Type_O>
class gtPlusSolverCallBack
{
public:

    gtPlusSolverCallBack();
    virtual ~gtPlusSolverCallBack();

    // if true, current solver will exit
    virtual bool exit();

    virtual bool callBack(size_t iter, Array_Type_O& x);

    virtual void printInfo(std::ostream& os) const;

    bool exit_;
};

template <typename Array_Type_I, typename Array_Type_O>
class gtPlusSolver
{
public:

    typedef gtPlusSolver<Array_Type_I, Array_Type_O> Self;

    typedef typename Array_Type_I::value_type ValueType;

    typedef gtPlusSolverCallBack<Array_Type_I, Array_Type_O> CBType;

    gtPlusSolver();
    virtual ~gtPlusSolver();

    CBType* getCallBack();
    void setCallBack(CBType* pCB);

    virtual bool solve(const Array_Type_I& b, Array_Type_O& x) = 0;

    virtual void printInfo(std::ostream& os) const;

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer1_;
    Gadgetron::GadgetronTimer gt_timer2_;
    Gadgetron::GadgetronTimer gt_timer3_;

    bool performTiming_;

    bool printIter_;

    // exporter
    Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

    // debug folder
    std::string debugFolder_;

    // util
    gtPlusISMRMRDReconUtil<ValueType> gtPlus_util_;
    gtPlusISMRMRDReconUtilComplex<ValueType> gtPlus_util_complex_;

    // memory manager
    boost::shared_ptr<gtPlusMemoryManager> gtPlus_mem_manager_;

protected:

    CBType* callback_;
};

template <typename Array_Type_I, typename Array_Type_O>
gtPlusSolverCallBack<Array_Type_I, Array_Type_O>::
gtPlusSolverCallBack() : exit_(true)
{

}

template <typename Array_Type_I, typename Array_Type_O>
gtPlusSolverCallBack<Array_Type_I, Array_Type_O>::
~gtPlusSolverCallBack() 
{

}

template <typename Array_Type_I, typename Array_Type_O>
bool gtPlusSolverCallBack<Array_Type_I, Array_Type_O>::
exit()
{
    return exit_;
}

template <typename Array_Type_I, typename Array_Type_O>
bool gtPlusSolverCallBack<Array_Type_I, Array_Type_O>::
callBack(size_t /*iter*/, Array_Type_O& /*x*/) 
{
    return true;
}

template <typename Array_Type_I, typename Array_Type_O>
void gtPlusSolverCallBack<Array_Type_I, Array_Type_O>::
printInfo(std::ostream& os) const
{
    using namespace std;
    os << "-------------- GTPlus ISMRMRD solver callback ------------------" << endl;
    os << "A callback scheme for ISMRMRD solvers ... " << endl;
    os << "----------------------------------------------------------------" << endl;
}

template <typename Array_Type_I, typename Array_Type_O>
gtPlusSolver<Array_Type_I, Array_Type_O>::
gtPlusSolver() : callback_(NULL), performTiming_(false), printIter_(false)
{
    gt_timer1_.set_timing_in_destruction(false);
    gt_timer2_.set_timing_in_destruction(false);
    gt_timer3_.set_timing_in_destruction(false);
}

template <typename Array_Type_I, typename Array_Type_O>
gtPlusSolver<Array_Type_I, Array_Type_O>::
~gtPlusSolver() 
{

}

template <typename Array_Type_I, typename Array_Type_O>
void gtPlusSolver<Array_Type_I, Array_Type_O>::
printInfo(std::ostream& os) const
{
    using namespace std;
    os << "-------------- GTPlus ISMRMRD solver ------------------" << endl;
    os << "GtPlus ISMRMRD solvers ... " << endl;
    os << "-------------------------------------------------------" << endl;
}

}}
