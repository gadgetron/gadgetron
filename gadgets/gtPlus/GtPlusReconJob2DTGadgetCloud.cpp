
#include "GtPlusReconJob2DTGadgetCloud.h"
#include "GtPlusGadgetOpenMP.h"
#include "gadgetron_paths.h"
#include "log.h"

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

GtPlusReconJob2DTGadgetCloud::GtPlusReconJob2DTGadgetCloud()
{
    debugFolder_ = "DebugOutput";

    performTiming_ = true;

    verboseMode_ = false;

    gt_timer1_.set_timing_in_destruction(false);
    gt_timer2_.set_timing_in_destruction(false);
    gt_timer3_.set_timing_in_destruction(false);

    filterRO_type_ = ISMRMRD_FILTER_GAUSSIAN;
    filterRO_sigma_ = 1.5;
    filterRO_width_ = 0.15;

    filterE1_type_ = ISMRMRD_FILTER_GAUSSIAN;
    filterE1_sigma_ = 1.5;
    filterE1_width_ = 0.15;

    filterE2_type_ = ISMRMRD_FILTER_GAUSSIAN;
    filterE2_sigma_ = 1.5;
    filterE2_width_ = 0.15;

    filterRO_ref_type_ = ISMRMRD_FILTER_HANNING;
    filterRO_ref_sigma_ = 1.5;
    filterRO_ref_width_ = 0.15;

    filterE1_ref_type_ = ISMRMRD_FILTER_HANNING;
    filterE1_ref_sigma_ = 1.5;
    filterE1_ref_width_ = 0.15;

    filterE2_ref_type_ = ISMRMRD_FILTER_HANNING;
    filterE2_ref_sigma_ = 1.5;
    filterE2_ref_width_ = 0.15;

    filterRO_pf_type_ = ISMRMRD_FILTER_HANNING;
    filterRO_pf_sigma_ = 1.5;
    filterRO_pf_width_ = 0.15;
    filterRO_pf_densityComp_ = false;

    filterE1_pf_type_ = ISMRMRD_FILTER_HANNING;
    filterE1_pf_sigma_ = 1.5;
    filterE1_pf_width_ = 0.15;
    filterE1_pf_densityComp_ = false;

    filterE2_pf_type_ = ISMRMRD_FILTER_HANNING;
    filterE2_pf_sigma_ = 1.5;
    filterE2_pf_width_ = 0.15;
    filterE2_pf_densityComp_ = false;

    process_config_called_ = false;

    Gadgetron::prepOpenMP();
}

GtPlusReconJob2DTGadgetCloud::~GtPlusReconJob2DTGadgetCloud()
{

}

bool GtPlusReconJob2DTGadgetCloud::readParameters()
{
    try
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "------> GtPlusReconJob2DTGadgetCloud parameters <------");

        debugFolder_ = debugFolder.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "debugFolder_ is " << debugFolder_);

        debugFolder2_ = debugFolder2.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "debugFolder2_ is " << debugFolder2_);

        cloud_node_file_ = cloudNodeFile.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "cloud_node_file_ is " << cloud_node_file_);

        performTiming_ = performTiming.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "performTiming_ is " << performTiming_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        // kspace filter parameters
        filterRO_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(filterRO.value());
        filterRO_sigma_ = filterRO_sigma.value();
        filterRO_width_ = filterRO_width.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_type_ is " << filterRO.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_sigma_ is " << filterRO_sigma_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_width_ is " << filterRO_width_);

        filterE1_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(filterE1.value());
        filterE1_sigma_ = filterE1_sigma.value();
        filterE1_width_ = filterE1_width.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_type_ is " << filterE1.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_sigma_ is " << filterE1_sigma_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_width_ is " << filterE1_width_);

        filterE2_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(filterE2.value());
        filterE2_sigma_ = filterE2_sigma.value();
        filterE2_width_ = filterE2_width.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_type_ is " << filterE2.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_sigma_ is " << filterE2_sigma_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_width_ is " << filterE2_width_);

        filterRO_ref_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(filterRefRO.value());
        filterRO_ref_sigma_ = filterRefRO_sigma.value();
        filterRO_ref_width_ = filterRefRO_width.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_ref_type_ is " << filterRefRO.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_ref_sigma_ is " << filterRO_ref_sigma_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_ref_width_ is " << filterRO_ref_width_);

        filterE1_ref_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(filterRefE1.value());
        filterE1_ref_sigma_ = filterRefE1_sigma.value();
        filterE1_ref_width_ = filterRefE1_width.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_ref_type_ is " << filterRefE1.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_ref_sigma_ is " << filterE1_ref_sigma_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_ref_width_ is " << filterE1_ref_width_);

        filterE2_ref_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(filterRefE2.value());
        filterE2_ref_sigma_ = filterRefE2_sigma.value();
        filterE2_ref_width_ = filterRefE2_width.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_ref_type_ is " << filterRefE2.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_ref_sigma_ is " << filterE2_ref_sigma_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_ref_width_ is " << filterE2_ref_width_);

        filterRO_pf_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(filterPartialFourierRO.value());
        filterRO_pf_sigma_ = filterPartialFourierRO_sigma.value();
        filterRO_pf_width_ = filterPartialFourierRO_width.value();
        filterRO_pf_densityComp_ = filterPartialFourierRO_densityComp.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_type_ is " << filterPartialFourierRO.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_sigma_ is " << filterRO_pf_sigma_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_width_ is " << filterRO_pf_width_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterRO_pf_densityComp_ is " << filterRO_pf_densityComp_);

        filterE1_pf_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(filterPartialFourierE1.value());
        filterE1_pf_sigma_ = filterPartialFourierE1_sigma.value();
        filterE1_pf_width_ = filterPartialFourierE1_width.value();
        filterE1_pf_densityComp_ = filterPartialFourierE1_densityComp.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_type_ is " << filterPartialFourierE1.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_sigma_ is " << filterE1_pf_sigma_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_width_ is " << filterE1_pf_width_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE1_pf_densityComp_ is " << filterE1_pf_densityComp_);

        filterE2_pf_type_ = gtPlus_util_.getISMRMRDKSpaceFilterFromName(filterPartialFourierE2.value());
        filterE2_pf_sigma_ = filterPartialFourierE2_sigma.value();
        filterE2_pf_width_ = filterPartialFourierE2_width.value();
        filterE2_pf_densityComp_ = filterPartialFourierE2_densityComp.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_type_ is " << filterPartialFourierE2.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_sigma_ is " << filterE2_pf_sigma_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_width_ is " << filterE2_pf_width_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "filterE2_pf_densityComp_ is " << filterE2_pf_densityComp_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        job_split_by_S_ = this->get_bool_value("job_split_by_S");
        job_num_of_N_ = (size_t)(this->get_int_value("job_num_of_N"));
        job_max_Megabytes_ = (size_t)(this->get_int_value("job_max_Megabytes"));
        job_overlap_ = (size_t)(this->get_int_value("job_overlap"));

        GDEBUG_CONDITION_STREAM(verboseMode_, "job_split_by_S_ is " << job_split_by_S_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "job_num_of_N_ is " << job_num_of_N_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "job_max_Megabytes_ is " << job_max_Megabytes_);
        GDEBUG_CONDITION_STREAM(verboseMode_, "job_overlap_ is " << job_overlap_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        CloudComputing_ = CloudComputing.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "CloudComputing_ is " << CloudComputing_);

        cloud_node_file_ = cloudNodeFile.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "cloud_node_file_ is " << cloud_node_file_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "CloudNodeXMLConfiguration is " << CloudNodeXMLConfiguration.value());

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusReconJob2DTGadgetCloud::readParameters() ... ");
        return false;
    }

    return true;
}

int GtPlusReconJob2DTGadgetCloud::process_config(ACE_Message_Block* mb)
{
    // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
    //   0  1  2   3    4   5    6     7  8   9

    verboseMode_ = verboseMode.value();

    // read in parameters from the xml
    GADGET_CHECK_RETURN(this->readParameters(), GADGET_FAIL);

    // generate the destination folder
    if ( !debugFolder_.empty() )
    {
        getDebugFolderPath(debugFolder_, debugFolder_fullPath_, verboseMode_);
    }
    else
    {
        GDEBUG_STREAM("GtPlusRecon, debugFolder is not set ...");
    }

    if ( !debugFolder2_.empty() )
    {
        getDebugFolderPath(debugFolder2_, debugFolder2_fullPath_, verboseMode_);
    }
    else
    {
        GDEBUG_STREAM("GtPlusRecon, debugFolder2 is not set ...");
    }

    return GADGET_OK;
}

bool GtPlusReconJob2DTGadgetCloud::setWorkOrder2DTParameters(GtPlusRecon2DTPara& para, WorkOrder2DTType* workOrder)
{
    workOrder->recon_kspace_needed_ = para.recon_kspace_needed_;

    if ( para.workOrderPara_.coil_compression_thres_>0 || para.workOrderPara_.coil_compression_num_modesKept_>0 )
    {
        workOrder->coil_compression_ = true;
    }
    else
    {
        workOrder->coil_compression_ = false;
    }

    workOrder->same_coil_compression_coeff_allS_ = para.same_coil_compression_coeff_allS_;

    workOrder->embedded_averageall_ref_ = para.embedded_averageall_ref_;
    workOrder->embedded_ref_numOfModes_ = para.embedded_ref_numOfModes_;
    workOrder->embedded_fullres_coilmap_ = para.embedded_fullres_coilmap_;
    workOrder->embedded_fullres_coilmap_useHighestSignal_ = para.embedded_fullres_coilmap_useHighestSignal_;
    workOrder->embedded_same_combinationcoeff_allS_ = para.embedded_same_combinationcoeff_allS_;
    workOrder->embedded_whichS_combinationcoeff_ = para.embedded_whichS_combinationcoeff_;
    workOrder->embedded_ref_fillback_ = para.embedded_ref_fillback_;

    workOrder->separate_averageall_ref_ = para.separate_averageall_ref_;
    workOrder->separate_ref_numOfModes_ = para.separate_ref_numOfModes_;
    workOrder->separate_fullres_coilmap_ = para.separate_fullres_coilmap_;
    workOrder->separate_same_combinationcoeff_allS_ = para.separate_same_combinationcoeff_allS_;
    workOrder->separate_whichS_combinationcoeff_ = para.separate_whichS_combinationcoeff_;

    workOrder->interleaved_same_combinationcoeff_allS_ = para.interleaved_same_combinationcoeff_allS_;
    workOrder->interleaved_whichS_combinationcoeff_ = para.interleaved_whichS_combinationcoeff_;
    workOrder->interleaved_ref_numOfModes_ = para.interleaved_ref_numOfModes_;

    workOrder->no_acceleration_averageall_ref_ = para.no_acceleration_averageall_ref_;
    workOrder->no_acceleration_ref_numOfModes_ = para.no_acceleration_ref_numOfModes_;
    workOrder->no_acceleration_same_combinationcoeff_allS_ = para.no_acceleration_same_combinationcoeff_allS_;
    workOrder->no_acceleration_whichS_combinationcoeff_ = para.no_acceleration_whichS_combinationcoeff_;

    return true;
}

bool GtPlusReconJob2DTGadgetCloud::parseGTCloudNodeFile(const std::string& filename, CloudType& gtCloud)
{
    std::string nodeFileName = get_gadgetron_home();
    nodeFileName.append("/share/gadgetron/config/gtCloud/");
    nodeFileName.append(filename);
    GDEBUG_CONDITION_STREAM(verboseMode_, "Cloud node file name is " << nodeFileName);

    std::ifstream fs(nodeFileName.c_str(), std::ios::in);
    if (!fs.is_open()) 
    {
        GWARN_STREAM("Cannot open GT CloudNodeFile; use the local setting instead ... ");
        return false;
    }

    // control node hostname
    std::string controlNode;
    fs >> controlNode;

    std::string portControlNode;
    fs >> portControlNode;

    // number of GadgetLevel nodes
    unsigned int num;
    fs >> num;

    gtCloud.resize(num);

    unsigned int n;
    for ( n=0; n<num; n++ )
    {
        std::string gadgetNode;
        fs >> gadgetNode;

        std::string portGadgetNode;
        fs >> portGadgetNode;

        std::string xmlGadgetNode;
        fs >> xmlGadgetNode;

        unsigned int computingPowerIndex;
        fs >> computingPowerIndex;

        gtCloud[n].get<0>() = gadgetNode;
        gtCloud[n].get<1>() = portGadgetNode;
        gtCloud[n].get<2>() = xmlGadgetNode;
        gtCloud[n].get<3>() = computingPowerIndex;

        GDEBUG_CONDITION_STREAM(verboseMode_, "Gadget Node " << n << " : " << gt_cloud_[n]);
    }

    fs.close();

    return true;
}

int GtPlusReconJob2DTGadgetCloud::process(Gadgetron::GadgetContainerMessage< int >* m1, Gadgetron::GadgetContainerMessage< GtPlusRecon2DTCloudPackageCPFL > * m2)
{
    // because the parameter configuration will not be sent, we need to call process_config explicitly
    if ( !process_config_called_ )
    {
        GADGET_CHECK_RETURN( (this->process_config(m1)==0), GADGET_FAIL);
        process_config_called_ = true;

        if ( CloudComputing_ )
        {
            bool parseSuccess = this->parseGTCloudNodeFile(cloud_node_file_, gt_cloud_);
            if ( parseSuccess )
            {
                CloudComputing_ = true;
                CloudSize_ = (int)gt_cloud_.size();

                if ( CloudSize_ == 0 )
                {
                    CloudComputing_ = false;
                    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusReconJob2DTGadgetCloud : cannot find algorithm nodes ... ");
                }
            }
        }
    }
    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusReconJob2DTGadgetCloud::process(...) starts ... ");

    int* jobID = m1->getObjectPtr();
    GDEBUG_CONDITION_STREAM(verboseMode_, "--> arriving job : " << *jobID << " ... ");

    GtPlusRecon2DTCloudPackageCPFL* job = m2->getObjectPtr();

    boost::shared_ptr< std::vector<size_t> > dims = job->kspace.get_dimensions();

    GDEBUG_CONDITION_STREAM(verboseMode_, "job array size : [Ro E1 Cha Slice E2 Con Phase Rep Set Seg] = [" 
        << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " << (*dims)[3] << " " << (*dims)[4] 
        << " " << (*dims)[5] << " " << (*dims)[6] << " " << (*dims)[7] << " " << (*dims)[8] << " " << (*dims)[9] << "]");

    GtPlusRecon2DTPara& para = job->para;

    // ---------------------------------------------------------
    // set the work flow
    // ---------------------------------------------------------
    workflow_.reconSizeRO_ = para.reconSizeRO_;
    workflow_.reconSizeE1_ = para.reconSizeE1_;
    workflow_.reconSizeE2_ = para.reconSizeE2_;
    workflow_.encodingFOV_RO_ = para.encodingFOV_RO_;
    workflow_.encodingFOV_E1_ = para.encodingFOV_E1_;
    workflow_.encodingFOV_E2_ = para.encodingFOV_E2_;
    workflow_.reconFOV_RO_ = para.reconFOV_RO_;
    workflow_.reconFOV_E1_ = para.reconFOV_E1_;
    workflow_.reconFOV_E2_ = para.reconFOV_E2_;

    // workflow_.dataDimStartingIndexes_ = workOrder->dataDimStartingIndexes_;
    workflow_.dim4th_ = para.dim_4th_;
    workflow_.dim5th_ = para.dim_5th_;
    workflow_.WorkOrderShareDim_ = para.workOrder_ShareDim_;
    workflow_.performTiming_ = performTiming_;

    // ---------------------------------------------------------
    // set work order
    // ---------------------------------------------------------
    WorkOrder2DTType workOrder;

    workOrder.copyFromPara(para.workOrderPara_);

    workOrder.job_split_by_S_ = job_split_by_S_;
    workOrder.job_num_of_N_ = job_num_of_N_;
    workOrder.job_max_Megabytes_ = job_max_Megabytes_;
    workOrder.job_overlap_ = job_overlap_;

    workOrder.CloudComputing_ = CloudComputing_;
    workOrder.CloudSize_ = CloudSize_;
    workOrder.gt_cloud_ = gt_cloud_;

    if ( workOrder.acceFactorE1_>1 && workOrder.CalibMode_==Gadgetron::ISMRMRD_interleaved )
    {
        Gadgetron::fillSampledLinesUpTo11DArray(job->kspace, workOrder.data_, job->timeStamp);
    }
    else
    {
        workOrder.data_ = job->kspace;
    }

    workOrder.time_stamp_ = job->timeStamp;
    workOrder.physio_time_stamp_ = job->physioTimeStamp;
    workOrder.ref_ = job->ref;

    // ---------------------------------------------------------
    // set the worker
    // ---------------------------------------------------------
    worker_grappa_.performTiming_ = performTiming_;
    if ( !debugFolder_fullPath_.empty() ) worker_grappa_.debugFolder_ = debugFolder_fullPath_;

    worker_noacceleration_.performTiming_ = performTiming_;
    if ( !debugFolder_fullPath_.empty() ) worker_noacceleration_.debugFolder_ = debugFolder_fullPath_;

    worker_spirit_.performTiming_ = performTiming_;
    if ( !debugFolder_fullPath_.empty() ) worker_spirit_.debugFolder_ = debugFolder_fullPath_;

    worker_spirit_L1_ncg_.performTiming_ = performTiming_;
    if ( !debugFolder_fullPath_.empty() ) worker_spirit_L1_ncg_.debugFolder_ = debugFolder_fullPath_;

    if ( !debugFolder_fullPath_.empty() ) workflow_.debugFolder_ = debugFolder_fullPath_;

    // perform the recon
    if ( performTiming_ ) { gt_timer1_.start("Recon 2DT workorder on cloud node ... "); }

    GADGET_CHECK_RETURN(this->generateKSpaceFilter(workOrder), GADGET_FAIL);

    workOrder.duplicate(workOrder_recon_);
    setWorkOrder2DTParameters(para, &workOrder_recon_);

    workflow_.workOrder_ = &workOrder_recon_;

    if ( verboseMode_ )
    {
        workflow_.workOrder_->print(std::cout);
    }

    workflow_.setDataArray(workOrder.data_, workOrder.time_stamp_, workOrder.physio_time_stamp_);

    if ( workOrder.ref_.get_number_of_elements() > 0 )
    {
        workflow_.setRefArray(workOrder.ref_);
    }
    else if ( para.workOrderPara_.CalibMode_==Gadgetron::ISMRMRD_interleaved )
    {
        workOrder.ref_ = workOrder.data_;
        workflow_.setRefArray(workOrder.ref_);
    }

    // set the work flow for worker and workOrder
    if ( workOrder.acceFactorE1_ > 1 )
    {
        if ( para.workOrderPara_.recon_algorithm_ == Gadgetron::ISMRMRD_SPIRIT )
        {
            workflow_.worker_ = &worker_spirit_;
        }
        else if ( para.workOrderPara_.recon_algorithm_ == Gadgetron::ISMRMRD_L1SPIRIT )
        {
            workflow_.worker_ = &worker_spirit_L1_ncg_;
        }
        else
        {
            workflow_.worker_ = &worker_grappa_;
        }
    }
    else
    {
        workflow_.worker_ = &worker_noacceleration_;
    }

    bool succeed = true;
    succeed = workflow_.preProcessing();
    if ( succeed )
    {
        succeed = workflow_.recon();
        if ( succeed )
        {
            succeed = workflow_.postProcessing();
        }
    }

    if ( !succeed )
    {
        GERROR_STREAM("GtPlusReconJob2DTGadgetCloud::process(...) failed... ");
    }

    if ( performTiming_ ) { gt_timer1_.stop(); }

    if ( !debugFolder2_fullPath_.empty() )
    {
        std::ostringstream ostr;
        ostr << "Node_Recon2DT_" << *jobID;

        hoNDArray< std::complex<float> > res = workflow_.res_;
        res.squeeze();
        if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder2_fullPath_+ostr.str()); }

        if ( workflow_.res_second_.get_number_of_elements() > 0 )
        {
            hoNDArray< std::complex<float> > res = workflow_.res_second_;
            res.squeeze();

            std::ostringstream ostr;
            ostr << "Node_Recon2DT_second_" << *jobID;

            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder2_fullPath_+ostr.str()); }
        }
    }

    // clean the kspace and ker and coil map
    job->kspace.clear();
    job->timeStamp.clear();
    job->physioTimeStamp.clear();
    job->ref.clear();

    if ( succeed )
    {
        job->complexIm = workflow_.res_;
        job->complexImSecond = workflow_.res_second_;
        job->resTimeStampSecond = workflow_.res_time_stamp_second_;
        job->resPhysioTimeStampSecond = workflow_.res_physio_time_stamp_second_;
    }
    else
    {
        job->complexIm.clear();
        job->res.clear();

        job->complexImSecond.clear();
        job->resTimeStampSecond.clear();
        job->resPhysioTimeStampSecond.clear();
    }

    // send out the results
    GADGET_CHECK_RETURN(this->sendOutJob(*jobID, job), GADGET_FAIL);

    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusReconJob2DTGadgetCloud::process(...) ends ... ");

    // reset the status
    workflow_.data_ = NULL;
    workflow_.time_stamp_ = NULL;
    workflow_.physio_time_stamp_ = NULL;
    workflow_.ref_ = NULL;
    workflow_.noise_ = NULL;
    workflow_.workOrder_ = NULL;
    // Gadgetron::clear(&workflow_.res_);

    m1->release();

    if ( this->verboseMode_ )
    {
        std::string procTime;
        gtPlus_util_.getCurrentMoment(procTime);

        GDEBUG_STREAM("* ============================================================================== *");
        GDEBUG_STREAM("---> MR recon 2DT gadget cloud, Currnt processing time : " << procTime << " <---");
        GDEBUG_STREAM("* ============================================================================== *");
    }

    return GADGET_OK;
}

bool GtPlusReconJob2DTGadgetCloud::
sendOutJob(int jobID, GtPlusRecon2DTCloudPackageCPFL* job)
{
    try
    {
      GDEBUG("GtPlusReconJob2DTGadgetCloud sendOutJob ...\n");

        if (!this->controller_)
        {
	  GERROR("Cannot return result to controller, no controller set\n");
	  return false;
        }

        GadgetContainerMessage<GadgetMessageIdentifier>* mb =
            new GadgetContainerMessage<GadgetMessageIdentifier>();

        mb->getObjectPtr()->id = GADGET_MESSAGE_GADGETCLOUD_JOB;

        GadgetContainerMessage<int>* m1 = new GadgetContainerMessage<int>();
        *(m1->getObjectPtr()) = jobID;

        GadgetContainerMessage<GtPlusRecon2DTCloudPackageCPFL>* m2 = new GadgetContainerMessage<GtPlusRecon2DTCloudPackageCPFL>();

        *(m2->getObjectPtr()) = *job;

        m1->cont(m2);
        mb->cont(m1);

        int ret =  this->controller_->output_ready(mb);
        if (ret < 0)
        {
            GDEBUG("Failed to return GtPlusReconJob2DTGadgetCloud job massage to controller\n");
            return false;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusReconJob2DTGadgetCloud::sendOutJob(...) ... ");
        return false;
    }

    return true;
}

bool GtPlusReconJob2DTGadgetCloud::
generateKSpaceFilter(WorkOrderType& workOrder)
{
    try
    {
        size_t RO = workOrder.data_.get_size(0);
        size_t E1 = workOrder.data_.get_size(1);
        size_t E2 = workOrder.data_.get_size(4);

        size_t RO_ref = workOrder.ref_.get_size(0);
        size_t E1_ref = workOrder.ref_.get_size(1);
        size_t E2_ref = workOrder.ref_.get_size(4);

        if ( workOrder.CalibMode_ == Gadgetron::ISMRMRD_interleaved )
        {
            RO_ref = RO;
            E1_ref = E1;
            E2_ref = E2;
        }

        // image data filter
        if ( RO>1 && filterRO_type_ != ISMRMRD_FILTER_NONE )
        {
            workOrder.filterRO_.create(RO);
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(RO, workOrder.start_RO_, workOrder.end_RO_, workOrder.filterRO_, filterRO_type_, filterRO_sigma_, (size_t)std::ceil(filterRO_width_*RO)));
            if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterRO_, debugFolder_fullPath_+"filterRO"); }
        }

        if ( E1>1 && filterE1_type_ != ISMRMRD_FILTER_NONE )
        {
            workOrder.filterE1_.create(E1);
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(E1, workOrder.start_E1_, workOrder.end_E1_, workOrder.filterE1_, filterE1_type_, filterE1_sigma_, (size_t)std::ceil(filterE1_width_*E1)));
            if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_, debugFolder_fullPath_+"filterE1"); }
        }

        if ( E2>1 && filterE2_type_ != ISMRMRD_FILTER_NONE )
        {
            workOrder.filterE2_.create(E2);
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(E2, workOrder.start_E2_, workOrder.end_E2_, workOrder.filterE2_, filterE2_type_, filterE2_sigma_, (size_t)std::ceil(filterE2_width_*E2)));
            if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE2_, debugFolder_fullPath_+"filterE2"); }
        }

        // ref data filter
        if ( workOrder.ref_.get_number_of_elements() > 0 )
        {
            size_t startRO(0), endRO(0), startE1(0), endE1(0), startE2(0), endE2(0);
            if ( E2_ref == 1 )
            {
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.detectSampledRegion2D(workOrder.ref_, startRO, endRO, startE1, endE1));
            }
            else
            {
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_complex_.detectSampledRegion3D(workOrder.ref_, startRO, endRO, startE1, endE1, startE2, endE2));
            }

            if ( (workOrder.CalibMode_ == ISMRMRD_interleaved) || (workOrder.CalibMode_ == ISMRMRD_embedded) )
            {
                // use the image data sample range
                startRO = workOrder.start_RO_; if ( startRO < 0 ) startRO=0;
                endRO = workOrder.end_RO_; if ( endRO < 0 ) endRO = RO_ref-1;
            }

            if ( RO_ref > 1 && filterRO_ref_type_ != ISMRMRD_FILTER_NONE )
            {
                workOrder.filterRO_ref_.create(RO_ref);
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilterForRef(RO_ref, startRO, endRO, workOrder.filterRO_ref_, filterRO_ref_type_, filterRO_ref_sigma_, (size_t)std::ceil(filterRO_ref_width_*RO_ref)));
                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterRO_ref_, debugFolder_fullPath_+"filterRO_ref"); }
            }

            if ( (workOrder.CalibMode_ == ISMRMRD_separate) || (workOrder.CalibMode_ == ISMRMRD_external) )
            {
                if ( E1_ref > 1 && filterE1_ref_type_ != ISMRMRD_FILTER_NONE )
                {
                    size_t len = endE1-startE1+1;
                    workOrder.filterE1_ref_.create(len);
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(len, 0, len-1, workOrder.filterE1_ref_, filterE1_ref_type_, filterE1_ref_sigma_, (size_t)std::ceil(filterE1_ref_width_*len)));
                    if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_ref_, debugFolder_fullPath_+"filterE1_ref"); }
                }

                if ( E2_ref > 1 && filterE2_ref_type_ != ISMRMRD_FILTER_NONE )
                {
                    size_t len = endE2-startE2+1;
                    workOrder.filterE2_ref_.create(len);
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilter(len, 0, len-1, workOrder.filterE2_ref_, filterE2_ref_type_, filterE2_ref_sigma_, (size_t)std::ceil(filterE2_ref_width_*len)));
                    if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE2_ref_, debugFolder_fullPath_+"filterE2_ref"); }
                }
            }
            else
            {
                // this makes sure for interleaved and embedded, the kspace filter is applied at correct lines
                if ( E1_ref > 1 && filterE1_ref_type_ != ISMRMRD_FILTER_NONE )
                {
                    size_t len = E1_ref;
                    workOrder.filterE1_ref_.create(len);
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilterForRef(len, startE1, endE1, workOrder.filterE1_ref_, filterE1_ref_type_, filterE1_ref_sigma_, (size_t)std::ceil(filterE1_ref_width_*len)));
                    if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_ref_, debugFolder_fullPath_+"filterE1_ref"); }
                }

                if ( E2_ref > 1 && filterE2_ref_type_ != ISMRMRD_FILTER_NONE )
                {
                    size_t len = E2_ref;
                    workOrder.filterE2_ref_.create(len);
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateSymmetricFilterForRef(len, startE2, endE2, workOrder.filterE2_ref_, filterE2_ref_type_, filterE2_ref_sigma_, (size_t)std::ceil(filterE2_ref_width_*len)));
                    if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE2_ref_, debugFolder_fullPath_+"filterE2_ref"); }
                }
            }
        }

        // partial fourier handling filter
        if ( RO>1 && workOrder.start_RO_>=0 && workOrder.end_RO_>0 )
        {
            workOrder.filterRO_partialfourier_.create(RO);
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateAsymmetricFilter(RO, workOrder.start_RO_, workOrder.end_RO_, workOrder.filterRO_partialfourier_, filterRO_pf_type_, (size_t)std::ceil(filterRO_pf_width_*RO), filterRO_pf_densityComp_));
            if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterRO_partialfourier_, debugFolder_fullPath_+"filterRO_partialfourier"); }
        }

        if ( E1>1 && workOrder.start_E1_>=0 && workOrder.end_E1_>0 )
        {
            workOrder.filterE1_partialfourier_.create(E1);
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateAsymmetricFilter(E1, workOrder.start_E1_, workOrder.end_E1_, workOrder.filterE1_partialfourier_, filterE1_pf_type_, (size_t)std::ceil(filterE1_pf_width_*E1), filterE1_pf_densityComp_));
            if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE1_partialfourier_, debugFolder_fullPath_+"filterE1_partialfourier"); }
        }

        if ( E2>1 && workOrder.start_E2_>=0 && workOrder.end_E2_>0 )
        {
            workOrder.filterE2_partialfourier_.create(E2);
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.generateAsymmetricFilter(E2, workOrder.start_E2_, workOrder.end_E2_, workOrder.filterE2_partialfourier_, filterE2_pf_type_, (size_t)std::ceil(filterE2_pf_width_*E2), filterE2_pf_densityComp_));
            if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(workOrder.filterE2_partialfourier_, debugFolder_fullPath_+"filterE2_partialfourier"); }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusReconJob2DTGadgetCloud::generateKSpaceFilter(...) ... ");
        return false;
    }

    return true;
}

GADGET_FACTORY_DECLARE(GtPlusReconJob2DTGadgetCloud)

}
