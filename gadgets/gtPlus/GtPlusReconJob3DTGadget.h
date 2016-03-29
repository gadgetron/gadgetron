/** \file   GtPlusReconJob3DTGadget.h
    \brief  This gadget serves as the working gadget for the single layer cloud for 3DT reconstruction.

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
#include "gtPlusISMRMRDReconWorkFlowCartesian3DT.h"
#include "gtPlusISMRMRDReconWorker3DTGRAPPA.h"
#include "gtPlusISMRMRDReconWorker3DTNoAcceleration.h"
#include "gtPlusISMRMRDReconWorker3DTSPIRIT.h"
#include "GtPlusReconGadgetUtil.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

class EXPORTGTPLUSGADGET GtPlusReconJob3DTGadget : public Gadgetron::Gadget2< int, GtPlusReconJobTypeCPFL >
{
public:
    GADGET_DECLARE(GtPlusReconJob3DTGadget);

    typedef std::complex<float> ValueType;

    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder<ValueType> WorkOrderType;

    typedef Gadget2< int, GtPlusReconJobTypeCPFL > BaseClass;

    GtPlusReconJob3DTGadget();
    ~GtPlusReconJob3DTGadget();

    // debug folder
    std::string debugFolder_;
    std::string debugFolder_fullPath_;

    std::string debugFolder2_;
    std::string debugFolder2_fullPath_;

    // whether to perform timing
    bool performTiming_;

    GADGET_PROPERTY(verboseMode, bool, "Whether to print more information", false);
    GADGET_PROPERTY(debugFolder, std::string, "If set, the debug output will be written out", "");
    GADGET_PROPERTY(debugFolder2, std::string, "If set, the debug output will be written out", "");
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
    Gadgetron::gtPlus::gtPlusReconWorker3DTGRAPPA<ValueType> worker_grappa_;
    Gadgetron::gtPlus::gtPlusReconWorker3DTNoAcceleration<ValueType> worker_noacceleration_;
    Gadgetron::gtPlus::gtPlusReconWorker3DTSPIRIT<ValueType> worker_spirit_;

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
