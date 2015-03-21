/** \file   GtPlusRecon3DTGadget.h
    \brief  This gadget encapsulates the reconstruction for 3DT cases.
    \author Hui Xue
*/

#pragma once

#include "GtPlusReconGadget.h"
#include "gtPlusISMRMRDReconWorkFlowCartesian3DT.h"
#include "gtPlusISMRMRDReconWorker3DTGRAPPA.h"
#include "gtPlusISMRMRDReconWorker3DTNoAcceleration.h"
#include "gtPlusISMRMRDReconWorker3DTSPIRIT.h"
#include "gtPlusISMRMRDReconWorker3DTL1SPIRITNCG.h"

namespace Gadgetron
{

struct EXPORTGTPLUSGADGET GtPlusRecon3DTPara
{
    size_t reconSizeRO_;
    size_t reconSizeE1_;
    size_t reconSizeE2_;

    float encodingFOV_RO_;
    float encodingFOV_E1_;
    float encodingFOV_E2_;

    float reconFOV_RO_;
    float reconFOV_E1_;
    float reconFOV_E2_;

    Gadgetron::ISMRMRDDIM dim_5th_;
    Gadgetron::ISMRMRDDIM workOrder_ShareDim_;

    bool no_acceleration_averageall_ref_;
    bool no_acceleration_same_combinationcoeff_allN_;
    int no_acceleration_whichN_combinationcoeff_;

    bool interleaved_same_combinationcoeff_allN_;
    int interleaved_whichN_combinationcoeff_;

    bool embedded_averageall_ref_;
    bool embedded_fullres_coilmap_;
    bool embedded_same_combinationcoeff_allN_;
    int embedded_whichN_combinationcoeff_;
    bool embedded_ref_fillback_;

    bool separate_averageall_ref_;
    bool separate_fullres_coilmap_;
    bool separate_same_combinationcoeff_allN_;
    int separate_whichN_combinationcoeff_;

    bool same_coil_compression_coeff_allN_;

    bool recon_kspace_needed_;

    Gadgetron::gtPlus::gtPlusReconWorkOrderPara workOrderPara_;
};

class EXPORTGTPLUSGADGET GtPlusRecon3DTGadget : public GtPlusReconGadget
{
public:
    GADGET_DECLARE(GtPlusRecon3DTGadget);

    typedef GtPlusReconGadget BaseClass;

    typedef BaseClass::ValueType ValueType;

    typedef BaseClass::WorkOrderType WorkOrderType;
    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder3DT<ValueType> WorkOrder3DTType;

    typedef BaseClass::DimensionRecordType DimensionRecordType;

    GtPlusRecon3DTGadget();
    ~GtPlusRecon3DTGadget();

    GtPlusRecon3DTPara para_;

protected:

    virtual bool readParameters();
    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(Gadgetron::GadgetContainerMessage< GtPlusGadgetImageArray >* m1, Gadgetron::GadgetContainerMessage< WorkOrderType > * m2);

    // set 3DT specific work order parameters
    bool setWorkOrder3DTParameters(WorkOrder3DTType* workOrder);

    // work flow
    Gadgetron::gtPlus::gtPlusISMRMRDReconWorkFlowCartesian3DT<ValueType> workflow_;

    // worker
    Gadgetron::gtPlus::gtPlusReconWorker3DTGRAPPA<ValueType> worker_grappa_;
    Gadgetron::gtPlus::gtPlusReconWorker3DTNoAcceleration<ValueType> worker_noacceleration_;
    Gadgetron::gtPlus::gtPlusReconWorker3DTSPIRIT<ValueType> worker_spirit_;
    Gadgetron::gtPlus::gtPlusReconWorker3DTL1SPIRITNCG<ValueType> worker_spirit_L1_ncg_;

    // workOrder for recon
    WorkOrder3DTType workOrder_recon_;

    // workOrder for recon 'other' data
    WorkOrder3DTType workOrder_recon_other_;
};

}
