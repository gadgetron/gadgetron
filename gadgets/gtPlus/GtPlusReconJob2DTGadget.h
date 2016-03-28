/** \file   GtPlusReconJob2DTGadget.h
    \brief  This is a cloud gadget performing the computation for 2DT job data package.

            This gadget can either serve as the working gadget for the signle layer cloud, or it can work as the 
            second layer gadget for the dual layer cloud.

            Ref to: 

            Hui Xue, Souheil Inati, Thomas Sangild Sorensen, Peter Kellman, Michael S. Hansen. 
            Distributed MRI Reconstruction using Gadgetron based Cloud Computing. 
            Magenetic Resonance in Medicine, doi: 10.1002/mrm.25213.

    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "Gadget.h"
#include "GadgetStreamController.h"
#include "GadgetCloudJobMessageReadWrite.h"
#include "GadgetronTimer.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorkOrder.h"
#include "gtPlusISMRMRDReconWorkFlowCartesian2DT.h"
#include "gtPlusISMRMRDReconWorker2DTGRAPPA.h"
#include "gtPlusISMRMRDReconWorker2DTNoAcceleration.h"
#include "gtPlusISMRMRDReconWorker2DTSPIRIT.h"
#include "GtPlusReconGadgetUtil.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

class EXPORTGTPLUSGADGET GtPlusReconJob2DTGadget : public Gadgetron::Gadget2< int, GtPlusReconJobTypeCPFL >
{
public:
    GADGET_DECLARE(GtPlusReconJob2DTGadget);

    typedef std::complex<float> ValueType;

    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder<ValueType> WorkOrderType;

    typedef Gadget2< int, GtPlusReconJobTypeCPFL > BaseClass;

    GtPlusReconJob2DTGadget();
    ~GtPlusReconJob2DTGadget();

    // debug folder
    std::string debugFolder_;
    std::string debugFolder_fullPath_;

    // whether to perform timing
    bool performTiming_;

    GADGET_PROPERTY(verboseMode, bool, "Whether to print more information", false);
    GADGET_PROPERTY(debugFolder, std::string, "If set, the debug output will be written out", "");
    GADGET_PROPERTY(performTiming, bool, "Whether to perform timing on some computational steps", false);

protected:

    // --------------------------------------------------
    // functional functions
    // --------------------------------------------------

    // default interface function
    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(Gadgetron::GadgetContainerMessage< int >* m1, Gadgetron::GadgetContainerMessage< GtPlusReconJobTypeCPFL > * m2);

    // process config is only to be called once
    bool process_config_called_;

    // read in parameters
    virtual bool readParameters();

    // send the completed job
    bool sendOutJob(int jobID, GtPlusReconJobTypeCPFL* job);

    // worker
    Gadgetron::gtPlus::gtPlusReconWorker2DTGRAPPA<ValueType> worker_grappa_;
    Gadgetron::gtPlus::gtPlusReconWorker2DTNoAcceleration<ValueType> worker_noacceleration_;
    Gadgetron::gtPlus::gtPlusReconWorker2DTSPIRIT<ValueType> worker_spirit_;

    Gadgetron::gtPlus::gtPlusISMRMRDReconUtil<ValueType> gtPlus_util_;
    Gadgetron::gtPlus::gtPlusISMRMRDReconUtilComplex<ValueType> gtPlus_util_complex_;

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer1_;
    Gadgetron::GadgetronTimer gt_timer2_;
    Gadgetron::GadgetronTimer gt_timer3_;

    // exporter
    Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

    // in verbose mode, more info is printed out
    bool verboseMode_;
};

}
