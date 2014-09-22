/** \file       gtPlusAlgorithmBase.h
    \brief      Base class for GtPlus algorithm
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

template <typename T> 
class gtPlusAlgorithmBase
{
public:

    gtPlusAlgorithmBase();
    virtual ~gtPlusAlgorithmBase();

    virtual void printInfo(std::ostream& os);

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

    // memory manager
    boost::shared_ptr<gtPlusMemoryManager> gtPlus_mem_manager_;
};

template <typename T> 
gtPlusAlgorithmBase<T>::gtPlusAlgorithmBase() : performTiming_(false)
{
    gt_timer1_.set_timing_in_destruction(false);
    gt_timer2_.set_timing_in_destruction(false);
    gt_timer3_.set_timing_in_destruction(false);
}

template <typename T> 
gtPlusAlgorithmBase<T>::~gtPlusAlgorithmBase()
{
}

template <typename T> 
void gtPlusAlgorithmBase<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD Algorithm ------------------" << endl;
    os << "Implementation of algorithms for ISMRMRD package" << endl;
    os << "----------------------------------------------------------" << endl;
}

}}
