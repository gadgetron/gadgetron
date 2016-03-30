/** \file   GtPlusRecon2DTGadget.h
    \brief  This gadget encapsulates the reconstruction for 2DT cases.
    \author Hui Xue
*/

#pragma once

#include "GtPlusReconGadget.h"
#include "gtPlusISMRMRDReconWorkFlowCartesian2DT.h"
#include "gtPlusISMRMRDReconWorker2DTGRAPPA.h"
#include "gtPlusISMRMRDReconWorker2DTSPIRIT.h"
#include "gtPlusISMRMRDReconWorker2DTNoAcceleration.h"
#include "GtPlusRecon2DTCloudPackage.h"

namespace Gadgetron
{

class EXPORTGTPLUSGADGET GtPlusRecon2DTGadget : public GtPlusReconGadget
{
public:
    GADGET_DECLARE(GtPlusRecon2DTGadget);

    typedef GtPlusReconGadget BaseClass;

    typedef BaseClass::ValueType ValueType;

    typedef BaseClass::WorkOrderType WorkOrderType;
    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder2DT<ValueType> WorkOrder2DTType;

    typedef BaseClass::DimensionRecordType DimensionRecordType;

    GtPlusRecon2DTGadget();
    ~GtPlusRecon2DTGadget();

    GtPlusRecon2DTPara para_;

    GADGET_PROPERTY_LIMITS(dim_4th, std::string, "The fourth dimension for 2DT recon", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY_LIMITS(dim_5th, std::string, "The fifth dimension for 2DT recon", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY_LIMITS(workOrder_ShareDim, std::string, "The workOrder coefficients sharing dimension for 2DT recon", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY(no_acceleration_averageall_ref, bool, "Whether to use average-all ref strategy, no accleration mode", true);
    GADGET_PROPERTY(no_acceleration_ref_numOfModes, int, "Number of moded kept for average-all ref strategy, no accleration mode", 0);
    GADGET_PROPERTY(no_acceleration_same_combinationcoeff_allS, bool, "Whether to use same coefficients for all S dimension, no accleration mode", true);
    GADGET_PROPERTY(no_acceleration_whichS_combinationcoeff, int, "Which S is used for coefficient computation, no accleration mode", 0);

    GADGET_PROPERTY(interleaved_same_combinationcoeff_allS, bool, "Whether to use same coefficients for all S dimension, interleaved mode", true);
    GADGET_PROPERTY(interleaved_ref_numOfModes, int, "Number of moded kept for interleaved ref, interleaved mode", 0);
    GADGET_PROPERTY(interleaved_whichS_combinationcoeff, int, "Which S is used for coefficient computation, interleaved mode", 0);

    GADGET_PROPERTY(embedded_averageall_ref, bool, "Whether to use average-all ref strategy, embedded mode", true);
    GADGET_PROPERTY(embedded_ref_numOfModes, int, "Number of moded kept for average-all ref strategy, embedded mode", 0);
    GADGET_PROPERTY(embedded_fullres_coilmap, bool, "Whether to compute full resolution coil map, embedded mode", false);
    GADGET_PROPERTY(embedded_fullres_coilmap_useHighestSignal, bool, "Whether to compute full resolution coil map using the highest signal image, embedded mode", true);
    GADGET_PROPERTY(embedded_same_combinationcoeff_allS, bool, "Whether to use same coefficients for all S dimension, embedded mode", true);
    GADGET_PROPERTY(embedded_whichS_combinationcoeff, int, "Which S is used for coefficient computation, embedded mode", 0);
    GADGET_PROPERTY(embedded_ref_fillback, bool, "Whether to fill back ref kspace lines back to the full kspace reconstructed, embedded mode", true);

    GADGET_PROPERTY(separate_averageall_ref, bool, "Whether to use average-all ref strategy, seperate mode", true);
    GADGET_PROPERTY(separate_ref_numOfModes, int, "Number of moded kept for average-all ref strategy, seperate mode", 0);
    GADGET_PROPERTY(separate_fullres_coilmap, bool, "Whether to compute full resolution coil map, seperate mode", false);
    GADGET_PROPERTY(separate_same_combinationcoeff_allS, bool, "Whether to use same coefficients for all S dimension, seperate mode", true);
    GADGET_PROPERTY(separate_whichS_combinationcoeff, int, "Which S is used for coefficient computation, seperate mode", 0);

    GADGET_PROPERTY(same_coil_compression_coeff_allS, bool, "Whether to use same coil compression coefficients for all S dimension", true);

protected:

    virtual bool readParameters();
    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(Gadgetron::GadgetContainerMessage< GtPlusGadgetImageArray >* m1, Gadgetron::GadgetContainerMessage< WorkOrderType > * m2);

    // set 2DT specific work order parameters
    bool setWorkOrder2DTParameters(WorkOrder2DTType* workOrder);

    // work flow
    Gadgetron::gtPlus::gtPlusISMRMRDReconWorkFlowCartesian2DT<ValueType> workflow_;

    // worker
    Gadgetron::gtPlus::gtPlusReconWorker2DTGRAPPA<ValueType> worker_grappa_;
    Gadgetron::gtPlus::gtPlusReconWorker2DTNoAcceleration<ValueType> worker_noacceleration_;
    Gadgetron::gtPlus::gtPlusReconWorker2DTSPIRIT<ValueType> worker_spirit_;

    // workOrder for recon
    WorkOrder2DTType workOrder_recon_;

    // workOrder for recon 'other' data
    WorkOrder2DTType workOrder_recon_other_;
};

}
