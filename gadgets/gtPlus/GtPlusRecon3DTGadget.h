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

    GADGET_PROPERTY_LIMITS(dim_5th, std::string, "The fifth dimension for 3DT recon", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY_LIMITS(workOrder_ShareDim, std::string, "The workOrder coefficients sharing dimension for 3DT recon", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY(no_acceleration_averageall_ref, bool, "Whether to use average-all ref strategy, no accleration mode", true);
    GADGET_PROPERTY(no_acceleration_same_combinationcoeff_allN, bool, "Whether to use same coefficients for all N dimension, no accleration mode", true);
    GADGET_PROPERTY(no_acceleration_whichN_combinationcoeff, int, "Which N is used for coefficient computation, no accleration mode", 0);

    GADGET_PROPERTY(interleaved_same_combinationcoeff_allN, bool, "Whether to use same coefficients for all N dimension, interleaved mode", true);
    GADGET_PROPERTY(interleaved_whichN_combinationcoeff, int, "Which N is used for coefficient computation, interleaved mode", 0);

    GADGET_PROPERTY(embedded_averageall_ref, bool, "Whether to use average-all ref strategy, embedded mode", true);
    GADGET_PROPERTY(embedded_fullres_coilmap, bool, "Whether to compute full resolution coil map, embedded mode", false);
    GADGET_PROPERTY(embedded_same_combinationcoeff_allN, bool, "Whether to use same coefficients for all N dimension, embedded mode", true);
    GADGET_PROPERTY(embedded_whichN_combinationcoeff, int, "Which N is used for coefficient computation, embedded mode", 0);
    GADGET_PROPERTY(embedded_ref_fillback, bool, "Whether to fill back ref kspace lines back to the full kspace reconstructed, embedded mode", true);

    GADGET_PROPERTY(separate_averageall_ref, bool, "Whether to use average-all ref strategy, seperate mode", true);
    GADGET_PROPERTY(separate_fullres_coilmap, bool, "Whether to compute full resolution coil map, seperate mode", false);
    GADGET_PROPERTY(separate_same_combinationcoeff_allN, bool, "Whether to use same coefficients for all N dimension, seperate mode", true);
    GADGET_PROPERTY(separate_whichN_combinationcoeff, int, "Which N is used for coefficient computation, seperate mode", 0);

    GADGET_PROPERTY(same_coil_compression_coeff_allN, bool, "Whether to use same coil compression coefficients for all N dimension", true);

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

    // workOrder for recon
    WorkOrder3DTType workOrder_recon_;

    // workOrder for recon 'other' data
    WorkOrder3DTType workOrder_recon_other_;
};

}
