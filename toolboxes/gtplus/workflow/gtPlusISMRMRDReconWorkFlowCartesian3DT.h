/** \file   gtPlusISMRMRDReconWorkFlowCartesian3DT.h
    \brief  Define the base class for the GtPlus 3DT reconstruction workflow for cartesian sampling
    \author Hui Xue
*/

#pragma once

#include "gtPlusISMRMRDReconWorkFlowCartesian.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusISMRMRDReconWorkFlowCartesian3DT : public gtPlusISMRMRDReconWorkFlowCartesian<T>
{
public:

    typedef gtPlusISMRMRDReconWorkFlowCartesian<T> BaseClass;
    typedef typename BaseClass::DimensionRecordType DimensionRecordType;

    gtPlusISMRMRDReconWorkFlowCartesian3DT();
    virtual ~gtPlusISMRMRDReconWorkFlowCartesian3DT();

    void printInfo(std::ostream& os);

    virtual bool recon();

    virtual bool predictDimensions();

    using BaseClass::data_;
    using BaseClass::ref_;
    using BaseClass::noise_;
    using BaseClass::noiseBW_;
    using BaseClass::receriverBWRatio_;
    using BaseClass::overSamplingRatioRO_;
    using BaseClass::reconSizeRO_;
    using BaseClass::reconSizeE1_;
    using BaseClass::reconSizeE2_;
    using BaseClass::encodingFOV_RO_;
    using BaseClass::encodingFOV_E1_;
    using BaseClass::encodingFOV_E2_;
    using BaseClass::reconFOV_RO_;
    using BaseClass::reconFOV_E1_;
    using BaseClass::reconFOV_E2_;
    using BaseClass::res_;

    using BaseClass::worker_;
    using BaseClass::workOrder_;

    using BaseClass::dataDimStartingIndexes_;

    using BaseClass::WorkOrderShareDim_;

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;

    using BaseClass::ref_remove_oversampling_RO_;
    using BaseClass::ref_apply_noisePreWhitening_;

    // the workOrder3DT needs 5 dimensions [RO E1 E2 CHA S]
    ISMRMRDDIM dim5th_;

protected:

    using BaseClass::dataCurr_;
    using BaseClass::refCurr_;

    using BaseClass::RO_;
    using BaseClass::E1_;
    using BaseClass::CHA_;
    using BaseClass::SLC_;
    using BaseClass::E2_;
    using BaseClass::CON_;
    using BaseClass::PHS_;
    using BaseClass::REP_;
    using BaseClass::SET_;
    using BaseClass::SEG_;

    using BaseClass::RO_ref_;
    using BaseClass::E1_ref_;
    using BaseClass::CHA_ref_;
    using BaseClass::SLC_ref_;
    using BaseClass::E2_ref_;
    using BaseClass::CON_ref_;
    using BaseClass::PHS_ref_;
    using BaseClass::REP_ref_;
    using BaseClass::SET_ref_;
    using BaseClass::SEG_ref_;

    using BaseClass::gtPlus_util_;
};

template <typename T> 
gtPlusISMRMRDReconWorkFlowCartesian3DT<T>::
gtPlusISMRMRDReconWorkFlowCartesian3DT() : BaseClass(), dim5th_(DIM_NONE)
{
}

template <typename T> 
gtPlusISMRMRDReconWorkFlowCartesian3DT<T>::~gtPlusISMRMRDReconWorkFlowCartesian3DT() 
{
}

template <typename T> 
void gtPlusISMRMRDReconWorkFlowCartesian3DT<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD Recon workflow Cartesian 3D/3DT -------------" << endl;
    os << "Implementation of general reconstruction workflow for cartesian sampling of 3D and 3D+T use cases" << endl;
    os << "The workOrder needs 5 dimensions [RO E1 E2 CHA S]" << endl;
    os << "----------------------------------------------------------" << endl;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian3DT<T>::predictDimensions()
{
    // if interleaved mode
    if ( workOrder_->CalibMode_ == ISMRMRD_interleaved )
    {
        if ( workOrder_->InterleaveDim_ == DIM_Phase )
        {
            dim5th_ = DIM_Phase;
        }

        if ( workOrder_->InterleaveDim_ == DIM_Repetition )
        {
            dim5th_ = DIM_Repetition;
        }

        if ( workOrder_->InterleaveDim_ == DIM_Contrast )
        {
            dim5th_ = DIM_Contrast;
        }
    }
    else if ( (workOrder_->CalibMode_ == ISMRMRD_embedded) 
        || (workOrder_->CalibMode_ == ISMRMRD_separate)
        || (workOrder_->CalibMode_ == ISMRMRD_noacceleration) ) 
    {
        if ( SLC_.second == 1 )
        {
            std::vector<DimensionRecordType> dimSizes(4);
            dimSizes[0] = CON_;
            dimSizes[1] = PHS_;
            dimSizes[2] = REP_;
            dimSizes[3] = SET_;

            std::sort(dimSizes.begin(), dimSizes.end(), DimensionRecordCompare() );
            dim5th_ = dimSizes[0].first;
        }

        if (SLC_.second > 1 )
        {
            dim5th_ = DIM_Slice; // multiple slab acquisition
        }
    }

    if ( dim5th_==DIM_NONE )
    {
        GADGET_ERROR_MSG("gtPlusISMRMRDReconWorkFlowCartesian3DT<T>::predictDimensions() : cannot find 5th dimensions ... ");
        return false;
    }

    workOrder_->enforceConsistency(dim5th_);

    GADGET_CONDITION_MSG(true, "predictDimensions - dim5th : " << gtPlus_util_.getISMRMRDDimName(dim5th_) );
    GADGET_CONDITION_MSG(true, "predictDimensions - WorkOrderShareDim : " << gtPlus_util_.getISMRMRDDimName(WorkOrderShareDim_) );

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian3DT<T>::recon()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data_!=NULL);
        GADGET_CHECK_RETURN_FALSE(worker_!=NULL);

        if ( dim5th_==DIM_NONE )
        {
            GADGET_CHECK_RETURN_FALSE(this->predictDimensions());
        }

        GADGET_CHECK_RETURN_FALSE(WorkOrderShareDim_!=DIM_ReadOut);
        GADGET_CHECK_RETURN_FALSE(WorkOrderShareDim_!=DIM_Encoding1);
        GADGET_CHECK_RETURN_FALSE(WorkOrderShareDim_!=DIM_Encoding2);
        GADGET_CHECK_RETURN_FALSE(WorkOrderShareDim_!=DIM_Channel);
        GADGET_CHECK_RETURN_FALSE(WorkOrderShareDim_!=dim5th_);

        // find recon dimensions
        std::vector<ISMRMRDDIM> dims;
        dims.push_back(DIM_ReadOut);
        dims.push_back(DIM_Encoding1);
        dims.push_back(DIM_Encoding2);
        dims.push_back(DIM_Channel);
        dims.push_back(dim5th_);

        int dim;
        size_t dd;

        int indWorkOrderSharingDim = -1;
        for ( dim=DIM_Slice; dim<=DIM_Set; dim++ )
        {
            bool exist = false;
            for ( dd=0; dd<dims.size(); dd++ )
            {
                if ( dims[dd] == dim )
                {
                    exist = true;
                    break;
                }
            }

            if ( !exist )
            {
                dims.push_back((ISMRMRDDIM)dim);

                if ( dim == WorkOrderShareDim_ )
                {
                    indWorkOrderSharingDim = dims.size()-1;
                }
            }
        }

        if ( (indWorkOrderSharingDim!=-1) && (indWorkOrderSharingDim > 5) )
        {
            ISMRMRDDIM dim6th = dims[5];
            dims[5] = WorkOrderShareDim_;
            dims[indWorkOrderSharingDim] = dim6th;
        }

        GADGET_CHECK_RETURN_FALSE(this->configureWorkOrder(dims));
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusISMRMRDReconWorkFlowCartesian3DT<T>::recon() ... ");
        return false;
    }

    return true;
}

}}
