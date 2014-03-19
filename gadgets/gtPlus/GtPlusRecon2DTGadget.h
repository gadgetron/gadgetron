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
#include "gtPlusISMRMRDReconWorker2DTL1SPIRITNCG.h"
#include "GtPlusRecon2DTCloudPackage.h"

namespace Gadgetron
{

class EXPORTGTPLUSGADGET GtPlusRecon2DTGadget : public GtPlusReconGadget
{
public:
    typedef GtPlusReconGadget BaseClass;

    typedef BaseClass::ValueType ValueType;

    typedef BaseClass::WorkOrderType WorkOrderType;
    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder2DT<ValueType> WorkOrder2DTType;

    typedef BaseClass::DimensionRecordType DimensionRecordType;

    GtPlusRecon2DTGadget();
    ~GtPlusRecon2DTGadget();

    GtPlusRecon2DTPara para_;

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
    Gadgetron::gtPlus::gtPlusReconWorker2DTL1SPIRITNCG<ValueType> worker_spirit_L1_ncg_;

    // workOrder for recon
    WorkOrder2DTType workOrder_recon_;

    // workOrder for recon 'other' data
    WorkOrder2DTType workOrder_recon_other_;
};

}
