/** \file   GtPlusReconJob2DTGadgetCloud.h
    \brief  This gadget serves as the first layer gadget for the dual layer cloud.

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

#include "hoNDArray_utils.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorkOrder.h"
#include "gtPlusISMRMRDReconWorkFlowCartesian2DT.h"
#include "gtPlusISMRMRDReconWorker2DTGRAPPA.h"
#include "gtPlusISMRMRDReconWorker2DTNoAcceleration.h"
#include "gtPlusISMRMRDReconWorker2DTSPIRIT.h"
#include "gtPlusISMRMRDReconWorker2DTL1SPIRITNCG.h"
#include "gtPlusMemoryManager.h"

#include "GtPlusRecon2DTCloudPackage.h"
#include "GtPlusReconGadgetUtil.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

class EXPORTGTPLUSGADGET GtPlusReconJob2DTGadgetCloud : public Gadgetron::Gadget2< int, GtPlusRecon2DTCloudPackageCPFL >
{
public:
    GADGET_DECLARE(GtPlusReconJob2DTGadgetCloud);

    typedef std::complex<float> ValueType;
    typedef Gadgetron::Gadget2< int, GtPlusRecon2DTCloudPackageCPFL > BaseClass;

    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder2DT<ValueType> WorkOrderType;
    typedef WorkOrderType WorkOrder2DTType;

    GtPlusReconJob2DTGadgetCloud();
    ~GtPlusReconJob2DTGadgetCloud();

    // kspace filter
    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterRO_type_;
    double filterRO_sigma_;
    double filterRO_width_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE1_type_;
    double filterE1_sigma_;
    double filterE1_width_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE2_type_;
    double filterE2_sigma_;
    double filterE2_width_;

    // ref data filter
    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterRO_ref_type_;
    double filterRO_ref_sigma_;
    double filterRO_ref_width_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE1_ref_type_;
    double filterE1_ref_sigma_;
    double filterE1_ref_width_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE2_ref_type_;
    double filterE2_ref_sigma_;
    double filterE2_ref_width_;

    // partial fourier filter
    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterRO_pf_type_;
    double filterRO_pf_sigma_;
    double filterRO_pf_width_;
    bool filterRO_pf_densityComp_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE1_pf_type_;
    double filterE1_pf_sigma_;
    double filterE1_pf_width_;
    bool filterE1_pf_densityComp_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE2_pf_type_;
    double filterE2_pf_sigma_;
    double filterE2_pf_width_;
    bool filterE2_pf_densityComp_;

    bool job_split_by_S_;
    size_t job_num_of_N_;
    size_t job_max_Megabytes_;
    size_t job_overlap_;

    /// cloud related definition
    bool CloudComputing_;
    unsigned int CloudSize_;

    typedef boost::tuple<std::string, std::string, std::string, unsigned int> CloudNodeType;
    typedef std::vector<CloudNodeType> CloudType;

    CloudType gt_cloud_;

    // debug folder
    std::string debugFolder_;
    std::string debugFolder_fullPath_;

    // debug folder 2
    std::string debugFolder2_;
    std::string debugFolder2_fullPath_;

    // cloud node file
    std::string cloud_node_file_;

    // whether to perform timing
    bool performTiming_;

protected:

    // --------------------------------------------------
    // functional functions
    // --------------------------------------------------

    // default interface function
    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(Gadgetron::GadgetContainerMessage< int >* m1, Gadgetron::GadgetContainerMessage< GtPlusRecon2DTCloudPackageCPFL > * m2);

    bool parseGTCloudNodeFile(const std::string& filename, CloudType& gtCloud);

    // process config is only to be called once
    bool process_config_called_;

    // read in parameters
    virtual bool readParameters();

    // send the completed job
    bool sendOutJob(int jobID, GtPlusRecon2DTCloudPackageCPFL* job);

    // set 2DT specific work order parameters
    bool setWorkOrder2DTParameters(GtPlusRecon2DTPara& para, WorkOrder2DTType* workOrder);

    // compute the kspace filter
    bool generateKSpaceFilter(WorkOrderType& workOrder);

    // work flow
    Gadgetron::gtPlus::gtPlusISMRMRDReconWorkFlowCartesian2DT<ValueType> workflow_;

    // worker
    Gadgetron::gtPlus::gtPlusReconWorker2DTGRAPPA<ValueType> worker_grappa_;
    Gadgetron::gtPlus::gtPlusReconWorker2DTNoAcceleration<ValueType> worker_noacceleration_;
    Gadgetron::gtPlus::gtPlusReconWorker2DTSPIRIT<ValueType> worker_spirit_;
    Gadgetron::gtPlus::gtPlusReconWorker2DTL1SPIRITNCG<ValueType> worker_spirit_L1_ncg_;

    // workOrder for recon
    WorkOrder2DTType workOrder_recon_;

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

    // memory manager
    boost::shared_ptr<Gadgetron::gtPlus::gtPlusMemoryManager> mem_manager_;
};

}
