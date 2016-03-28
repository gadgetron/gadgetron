
#include "GtPlusRecon2DTGadgetCloud.h"

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

GtPlusRecon2DTGadgetCloud::GtPlusRecon2DTGadgetCloud() : BaseClass(), curr_node_(0), num_of_jobs_(0)
{
    packages_sent_.resize(1024);
    packages_received_.resize(1024);
    packages_passed_to_next_gadget_.resize(1024);
    gt_timer_2DT_cloud_.set_timing_in_destruction(false);
}

GtPlusRecon2DTGadgetCloud::~GtPlusRecon2DTGadgetCloud()
{

}

int GtPlusRecon2DTGadgetCloud::process_config(ACE_Message_Block* mb)
{
    // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
    //   0  1  2   3    4   5    6     7  8   9
    GADGET_CHECK_RETURN(BaseClass::process_config(mb)==GADGET_OK, GADGET_FAIL);

    if ( CloudComputing_ )
    {
        bool parseSuccess = this->parseGTCloudNodeFile(cloud_node_file_, gt_cloud_);
        if ( parseSuccess )
        {
            CloudSize_ = (unsigned int)gt_cloud_.size();
            if ( CloudSize_ == 0 ) CloudComputing_ = false;
        }

        if ( CloudComputing_ )
        {
            // set up the cloud
            if (controller_.open () == -1)
            {
                GERROR_STREAM("Cloud controller cannot open the cloud ...");
                controller_.handle_close (ACE_INVALID_HANDLE, 0);
                CloudComputing_ = false;
            }
            else
            {
                readers_.resize(CloudSize_, NULL);
                writers_.resize(CloudSize_, NULL);

                unsigned int j;
                for ( j=0; j<CloudSize_; j++ )
                {
                    readers_[j] = new GtPlus2DTGadgetCloudJobMessageReaderCPFL();
                    writers_[j] = new GtPlus2DTGadgetCloudJobMessageWriterCPFL();
                }

                if ( controller_.createConnector(gt_cloud_, GADGET_MESSAGE_GADGETCLOUD_JOB, readers_, GADGET_MESSAGE_GADGETCLOUD_JOB, writers_) != 0 )
                {
                    GERROR_STREAM("Cloud controller_ creates connectors failed ...");
                    controller_.handle_close (ACE_INVALID_HANDLE, 0);
                    CloudComputing_ = false;
                }
                else if ( controller_.connectToCloud(gt_cloud_) != 0 )
                {
                    GERROR_STREAM("Cloud controller_ cannot connect to the cloud ...");
                    controller_.handle_close (ACE_INVALID_HANDLE, 0);
                    CloudComputing_ = false;
                }
            }
        }
    }

    return GADGET_OK;
}

int GtPlusRecon2DTGadgetCloud::process(Gadgetron::GadgetContainerMessage< GtPlusGadgetImageArray >* m1, Gadgetron::GadgetContainerMessage< WorkOrderType > * m2)
{
    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusRecon2DTGadgetCloud::process(...) starts ... ");

    processed_called_times_++;

    // start a gadget level timer
    if ( processed_called_times_ == 1 )
    {
        gt_timer_2DT_cloud_.start("GtPlusRecon2DTGadgetCloud::process(...) gadegt level timer ... ");
    }

    // send out the package to current node
    if ( CloudComputing_ )
    {
        GtPlusGadgetImageArray* images = m1->getObjectPtr();

        WorkOrderType* workOrder = m2->getObjectPtr();

        boost::shared_ptr< std::vector<size_t> > dims = workOrder->data_.get_dimensions();

        GDEBUG_CONDITION_STREAM(verboseMode_, "[Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave] = [" 
            << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " << (*dims)[3] << " " << (*dims)[4] 
            << " " << (*dims)[5] << " " << (*dims)[6] << " " << (*dims)[7] << " " << (*dims)[8] << " " << (*dims)[9] << " " << (*dims)[10] << "]");

        dimensions_ = *dims;

        // fill in more parameters
        para_.reconSizeRO_ = std::max(matrix_size_recon_[0], (*dims)[0]);
        para_.reconSizeE1_ = reconE1_;
        para_.reconSizeE2_ = reconE2_;
        para_.encodingFOV_RO_ = field_of_view_encoding_[0];
        para_.encodingFOV_E1_ = field_of_view_encoding_[1];
        para_.encodingFOV_E2_ = field_of_view_encoding_[2];
        para_.reconFOV_RO_ = field_of_view_recon_[0];
        para_.reconFOV_E1_ = field_of_view_recon_[1];
        para_.reconFOV_E2_ = field_of_view_recon_[2];

        para_.workOrderPara_.CalibMode_ = workOrder->CalibMode_;
        para_.workOrderPara_.InterleaveDim_ = workOrder->InterleaveDim_;

        para_.workOrderPara_.acceFactorE1_ = workOrder->acceFactorE1_;
        para_.workOrderPara_.acceFactorE2_ = workOrder->acceFactorE2_;

        para_.workOrderPara_.kSpaceCenterRO_ = workOrder->kSpaceCenterRO_;
        para_.workOrderPara_.kSpaceCenterEncode1_ = workOrder->kSpaceCenterEncode1_;
        para_.workOrderPara_.kSpaceCenterEncode2_ = workOrder->kSpaceCenterEncode2_;

        para_.workOrderPara_.kSpaceMaxRO_ = workOrder->kSpaceMaxRO_;
        para_.workOrderPara_.kSpaceMaxEncode1_ = workOrder->kSpaceMaxEncode1_;
        para_.workOrderPara_.kSpaceMaxEncode2_ = workOrder->kSpaceMaxEncode2_;

        para_.workOrderPara_.start_RO_ = workOrder->start_RO_;
        para_.workOrderPara_.end_RO_ = workOrder->end_RO_;

        para_.workOrderPara_.start_E1_ = workOrder->start_E1_;
        para_.workOrderPara_.end_E1_ = workOrder->end_E1_;

        para_.workOrderPara_.start_E2_ = workOrder->start_E2_;
        para_.workOrderPara_.end_E2_ = workOrder->end_E2_;

        para_.workOrderPara_.workFlow_BufferKernel_ = workOrder->workFlow_BufferKernel_;
        para_.workOrderPara_.workFlow_use_BufferedKernel_ = workOrder->workFlow_use_BufferedKernel_;
        para_.workOrderPara_.num_channels_res_ = workOrder->num_channels_res_;

        // set up a cloud package
        CloudPackageType package;
        package.para = para_;

        packages_sent_[num_of_jobs_] = package;
        packages_sent_[num_of_jobs_].kspace = workOrder->data_;

        packages_received_[num_of_jobs_] = package;

        packages_passed_to_next_gadget_[num_of_jobs_].first = num_of_jobs_;
        packages_passed_to_next_gadget_[num_of_jobs_].second = false;

        // store image headers
        GtPlusGadgetImageArray imArray;
        image_headers_.push_back(imArray);
        image_headers_[image_headers_.size()-1].copy(*images);

        // send the package to current node
        std::vector<CloudPackageType* > jobListCloud(1);
        std::vector<CloudPackageType* > completedJobListCloud(1);
        std::vector<int> node_ids(1, curr_node_);

        jobListCloud[0] = &packages_sent_[num_of_jobs_];
        completedJobListCloud[0] = &packages_received_[num_of_jobs_];

        // set the data and ref arrays

        // get the data to be compressed format
        if ( workOrder->acceFactorE1_>1 && workOrder->CalibMode_==Gadgetron::ISMRMRD_interleaved )
        {
            Gadgetron::extractSampledLinesUpTo11DArray(workOrder->data_, jobListCloud[0]->kspace, workOrder->time_stamp_, workOrder->acceFactorE1_, workOrder->acceFactorE2_);
        }
        else
        {
            jobListCloud[0]->kspace = workOrder->data_;
        }

        jobListCloud[0]->timeStamp = workOrder->time_stamp_;
        jobListCloud[0]->physioTimeStamp = workOrder->physio_time_stamp_;
        if ( workOrder->ref_.get_number_of_elements() > 0 )
        {
            jobListCloud[0]->ref = workOrder->ref_;
        }
        else if ( CalibMode_==Gadgetron::ISMRMRD_interleaved )
        {
            // jobListCloud[0]->ref = workOrder->data_;
            jobListCloud[0]->ref.clear();
        }

        num_of_jobs_++;

        if ( controller_.runJobsOnCloud(jobListCloud, completedJobListCloud, node_ids) != 0 )
        {
            GERROR_STREAM("Cloud controller runs jobs on the cloud failed ...");
            controller_.handle_close (ACE_INVALID_HANDLE, 0);

            // run locally
            int retval = BaseClass::process(m1, m2);
            packages_passed_to_next_gadget_[num_of_jobs_].second = true;

            return retval;
        }

        curr_node_++;
        if ( curr_node_ >= CloudSize_ ) curr_node_ = 0;

        m1->release();
    }
    else
    {
        return BaseClass::process(m1, m2);
    }

    return GADGET_OK;
}

bool GtPlusRecon2DTGadgetCloud::processJob(CloudPackageType& jobSent, CloudPackageType& jobReceived)
{
    try
    {
        GtPlusRecon2DTCloudPackageCPFL* job = &jobSent;

        boost::shared_ptr< std::vector<size_t> > dims = job->kspace.get_dimensions();

        GDEBUG_CONDITION_STREAM(verboseMode_, "job array size : [Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave] = [" 
            << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " << (*dims)[3] << " " << (*dims)[4] 
            << " " << (*dims)[5] << " " << (*dims)[6] << " " << (*dims)[7] << " " << (*dims)[8] << " " << (*dims)[9] << " " << (*dims)[10] << "]");

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

        workOrder.CloudComputing_ = CloudComputing_;
        workOrder.CloudSize_ = CloudSize_;
        workOrder.gt_cloud_ = gt_cloud_;

        if ( workOrder.acceFactorE1_ <= 1 )
        {
            workOrder.data_ = job->kspace;
        }
        else
        {
            Gadgetron::fillSampledLinesUpTo11DArray(job->kspace, workOrder.data_, job->timeStamp);
        }

        workOrder.time_stamp_ = job->timeStamp;
        workOrder.physio_time_stamp_ = job->physioTimeStamp;
        workOrder.ref_ = job->ref;

        // ---------------------------------------------------------
        // set the worker
        // ---------------------------------------------------------
        worker_grappa_.verbose_ = verboseMode_;
        worker_grappa_.performTiming_ = performTiming_;
        if ( !debugFolder_fullPath_.empty() ) worker_grappa_.debugFolder_ = debugFolder_fullPath_;

        worker_noacceleration_.verbose_ = verboseMode_;
        worker_noacceleration_.performTiming_ = performTiming_;
        if ( !debugFolder_fullPath_.empty() ) worker_noacceleration_.debugFolder_ = debugFolder_fullPath_;

        worker_spirit_.verbose_ = verboseMode_;
        worker_spirit_.performTiming_ = performTiming_;
        if ( !debugFolder_fullPath_.empty() ) worker_spirit_.debugFolder_ = debugFolder_fullPath_;

        if ( !debugFolder_fullPath_.empty() ) workflow_.debugFolder_ = debugFolder_fullPath_;

        if ( verboseMode_ )
        {
            workOrder.print(std::cout);
        }

        // perform the recon
        if ( performTiming_ ) { gt_timer1_.start("Recon 2DT workorder on master node ... "); }

        GADGET_CHECK_RETURN_FALSE(this->generateKSpaceFilter(workOrder));

        workOrder.duplicate(workOrder_recon_);
        this->setWorkOrder2DTParameters(&workOrder_recon_);

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

        if ( performTiming_ ) { gt_timer1_.stop(); }

        if ( !debugFolder_fullPath_.empty() )
        {
            std::ostringstream ostr;
            ostr << "Recon2DT";

            hoNDArray< std::complex<float> > res = workflow_.res_;
            res.squeeze();
            if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_fullPath_+ostr.str()); }

            if ( workflow_.res_second_.get_number_of_elements() > 0 )
            {
                hoNDArray< std::complex<float> > res = workflow_.res_second_;
                res.squeeze();

                std::ostringstream ostr;
                ostr << "Recon2DT_Second_" << processed_called_times_;

                if ( !debugFolder_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_fullPath_+ostr.str()); }

                if ( workflow_.res_time_stamp_second_.get_number_of_elements() > 0 )
                {
                    std::ostringstream ostr;
                    ostr << "Recon2DT_Second_TimeStamp_" << processed_called_times_;

                    hoNDArray<float> res = workflow_.res_time_stamp_second_;
                    res.squeeze();
                    if ( debugFolder_fullPath_.empty() ) { gt_exporter_.exportArray(res, debugFolder_fullPath_+ostr.str()); }
                }

                if ( workflow_.res_physio_time_stamp_second_.get_number_of_elements() > 0 )
                {
                    std::ostringstream ostr;
                    ostr << "Recon2DT_Second_PhysioTimeStamp_" << processed_called_times_;

                    hoNDArray<float> res = workflow_.res_physio_time_stamp_second_;
                    res.squeeze();
                    if ( debugFolder_fullPath_.empty() ) { gt_exporter_.exportArray(res, debugFolder_fullPath_+ostr.str()); }
                }
            }
        }

        if ( succeed )
        {
            jobReceived.complexIm = workflow_.res_;
            jobReceived.complexImSecond = workflow_.res_second_;
            jobReceived.resTimeStampSecond = workflow_.res_time_stamp_second_;
            jobReceived.resPhysioTimeStampSecond = workflow_.res_physio_time_stamp_second_;
        }
        else
        {
            jobReceived.complexIm.clear();
            jobReceived.complexImSecond.clear();
            jobReceived.resTimeStampSecond.clear();
            jobReceived.resPhysioTimeStampSecond.clear();
            jobReceived.res.clear();
        }

        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusRecon2DTGadgetCloud::process(...) ends ... ");

        // reset the status
        workflow_.data_ = NULL;
        workflow_.time_stamp_ = NULL;
        workflow_.physio_time_stamp_ = NULL;
        workflow_.ref_ = NULL;
        workflow_.noise_ = NULL;
        workflow_.workOrder_ = NULL;
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in GtPlusRecon2DTGadgetCloud::processJob(CloudPackageType& jobSent, CloudPackageType& jobReceived) ... ");
        return false;
    }

    return true;
}

int GtPlusRecon2DTGadgetCloud::close(unsigned long flags)
{
    GDEBUG_CONDITION_STREAM(true, "GtPlusRecon2DTGadgetCloud - close(flags) : " << flags);

    if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

    if ( flags!=0 )
    {
        GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusRecon2DTGadgetCloud number of total jobs : " << num_of_jobs_ << " ... ");

        if ( CloudComputing_ )
        {
            controller_.closeCloudNode();

            // register a job handler
            GtPlusRecon2DTGadgetCloudSender gadgetJobHandler;
            gadgetJobHandler.gadget_ = this;
            controller_.job_handler_ = &gadgetJobHandler;

            controller_.waitForJobToComplete();

            // if some jobs are not completed successfully, reprocess them; otherwise, send out images
            std::vector<DimensionRecordType> dataDimStartingIndexes;
            unsigned int N = (unsigned int)image_headers_.size();
            unsigned int ii;
            for ( ii=0; ii<N; ii++ )
            {
                bool jobIsOk = true;

                bool recomputeJob = (packages_received_[ii].complexIm.get_number_of_elements() == 0);

                // special check if the second set of recon results is needed
                if ( recon_res_second_required_ )
                {
                    GDEBUG_STREAM("Check received recon results (second set) ... ");

                    if (packages_received_[ii].complexImSecond.get_number_of_elements() == 0)
                    {
                        recomputeJob = true;
                    }
                    else
                    {
                        // check the images are not empty
                        real_value_type v(0);
                        Gadgetron::norm2(packages_received_[ii].complexImSecond, v);

                        if ( std::abs(v) < FLT_EPSILON )
                        {
                            recomputeJob = true;
                            GWARN_STREAM("Received recon results (second set) contain no content ... ");
                        }
                    }
                }

                if ( recomputeJob )
                {
                    // if the cloud goes wrong, do not try again
                    CloudComputing_ = false;
                    jobIsOk = this->processJob(packages_sent_[ii], packages_received_[ii]);
                }

                if ( jobIsOk )
                {
                    if ( !packages_passed_to_next_gadget_[ii].second )
                    {
                        GADGET_CHECK_RETURN(this->scalingImages(packages_received_[ii].complexIm), GADGET_FAIL);

                        if ( this->send_out_recon_ )
                        {
                            GADGET_CHECK_RETURN(this->sendOutRecon(&image_headers_[ii], packages_received_[ii].complexIm, image_series_, dataDimStartingIndexes, "Image", GADGETRON_IMAGE_REGULAR), GADGET_FAIL);
                        }

                        if ( this->send_out_recon_second_ )
                        {
                            if ( packages_received_[ii].complexImSecond.get_number_of_elements() > 0 )
                            {
                                Gadgetron::scal((float)scalingFactor_, packages_received_[ii].complexImSecond);

                                if ( this->para_.workOrderPara_.retro_gated_images_>0 )
                                {
                                    GADGET_CHECK_RETURN(this->sendOutRecon(&image_headers_[ii], 
                                                                            packages_received_[ii].complexImSecond, 
                                                                            packages_received_[ii].resTimeStampSecond, 
                                                                            packages_received_[ii].resPhysioTimeStampSecond, 
                                                                            image_series_+1, dataDimStartingIndexes, 
                                                                            "ImageRetro", GADGETRON_IMAGE_RETRO), GADGET_FAIL);
                                }
                                else
                                {
                                    GADGET_CHECK_RETURN(this->sendOutRecon(&image_headers_[ii], 
                                                                            packages_received_[ii].complexImSecond, 
                                                                            packages_received_[ii].resTimeStampSecond, 
                                                                            packages_received_[ii].resPhysioTimeStampSecond, 
                                                                            image_series_+1, dataDimStartingIndexes, 
                                                                            "Image", GADGETRON_IMAGE_REGULAR), GADGET_FAIL);
                                }
                            }
                        }
                    }
                }

                if ( !debugFolder2_fullPath_.empty() )
                {
                    std::ostringstream ostr;
                    ostr << "GadgetCloud_Recon2DT_" << ii;

                    hoNDArray< std::complex<float> > res = packages_received_[ii].complexIm;
                    res.squeeze();
                    if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder2_fullPath_+ostr.str()); }

                    if (packages_received_[ii].complexImSecond.get_number_of_elements() > 0 )
                    {
                        hoNDArray< std::complex<float> > res = packages_received_[ii].complexImSecond;
                        res.squeeze();

                        std::ostringstream ostr;
                        ostr << "GadgetCloud_Recon2DT_Second_" << ii;

                        if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder2_fullPath_+ostr.str()); }

                        if ( packages_received_[ii].resTimeStampSecond.get_number_of_elements() > 0 )
                        {
                            std::ostringstream ostr;
                            ostr << "GadgetCloud_Recon2DT_Second_TimeStamp_" << ii;

                            hoNDArray<float> res = packages_received_[ii].resTimeStampSecond;
                            res.squeeze();
                            if ( debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArray(res, debugFolder2_fullPath_+ostr.str()); }
                        }

                        if ( packages_received_[ii].resPhysioTimeStampSecond.get_number_of_elements() > 0 )
                        {
                            std::ostringstream ostr;
                            ostr << "GadgetCloud_Recon2DT_Second_PhysioTimeStamp_" << ii;

                            hoNDArray<float> res = packages_received_[ii].resPhysioTimeStampSecond;
                            res.squeeze();
                            if ( debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArray(res, debugFolder2_fullPath_+ostr.str()); }
                        }
                    }
                }
            }
        }

        gt_timer_2DT_cloud_.stop();
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GtPlusRecon2DTGadgetCloud)

// -------------------------------------------------------------------------------------------
// GtPlusRecon2DTGadgetCloudSender
// -------------------------------------------------------------------------------------------

GtPlusRecon2DTGadgetCloudSender::GtPlusRecon2DTGadgetCloudSender()
{
}

GtPlusRecon2DTGadgetCloudSender::~GtPlusRecon2DTGadgetCloudSender()
{
}

bool GtPlusRecon2DTGadgetCloudSender::processJob(int jobID, GtPlusRecon2DTCloudPackage< std::complex<float> >& ajob)
{
    try
    {
        bool jobIsOk = true;
        if ( (gadget_->packages_received_[jobID].complexIm.get_number_of_elements() == 0) 
            && (gadget_->packages_received_[jobID].res.get_number_of_elements() == 0) )
        {
            jobIsOk = false;
            return true;
        }

        if ( jobIsOk )
        {
            std::vector<DimensionRecordType> dataDimStartingIndexes;

            if ( !gadget_->packages_passed_to_next_gadget_[jobID].second )
            {
                gadget_->packages_passed_to_next_gadget_[jobID].second = true;

                GADGET_CHECK_RETURN(gadget_->scalingImages(gadget_->packages_received_[jobID].complexIm), false);

                if ( gadget_->send_out_recon_ )
                {
                    GADGET_CHECK_RETURN(gadget_->sendOutRecon(&gadget_->image_headers_[jobID], 
                        gadget_->packages_received_[jobID].complexIm, gadget_->image_series_, dataDimStartingIndexes, "Image", GADGETRON_IMAGE_REGULAR), false);
                }

                if ( gadget_->send_out_recon_second_ )
                {
                    if ( gadget_->packages_received_[jobID].complexImSecond.get_number_of_elements() > 0 )
                    {
                        GDEBUG_STREAM("Check received recon results (second set) in cloud sender ... ");

                        // check the images are not empty
                        float v(0);
                        Gadgetron::norm2(gadget_->packages_received_[jobID].complexImSecond, v);

                        bool reconResSecondValid = true;
                        if ( std::abs(v) < FLT_EPSILON )
                        {
                            reconResSecondValid = false;
                            GWARN_STREAM("Received recon results (second set) contain no content ... ");
                        }

                        if ( reconResSecondValid )
                        {
                            Gadgetron::scal((float)gadget_->scalingFactor_, gadget_->packages_received_[jobID].complexImSecond);
                            if ( gadget_->para_.workOrderPara_.retro_gated_images_ > 0 )
                            {
                                GADGET_CHECK_RETURN(gadget_->sendOutRecon(&gadget_->image_headers_[jobID], 
                                                                        gadget_->packages_received_[jobID].complexImSecond, 
                                                                        gadget_->packages_received_[jobID].resTimeStampSecond,
                                                                        gadget_->packages_received_[jobID].resPhysioTimeStampSecond,
                                                                        gadget_->image_series_+1, dataDimStartingIndexes, 
                                                                        "ImageRetro", GADGETRON_IMAGE_RETRO), false);
                            }
                            else
                            {
                                GADGET_CHECK_RETURN(gadget_->sendOutRecon(&gadget_->image_headers_[jobID], 
                                                                        gadget_->packages_received_[jobID].complexImSecond, 
                                                                        gadget_->packages_received_[jobID].resTimeStampSecond,
                                                                        gadget_->packages_received_[jobID].resPhysioTimeStampSecond,
                                                                        gadget_->image_series_+1, dataDimStartingIndexes, 
                                                                        "Image", GADGETRON_IMAGE_REGULAR), false);
                            }
                        }
                    }
                }

                if ( !gadget_->debugFolder2_fullPath_.empty() )
                {
                    std::ostringstream ostr;
                    ostr << "Recon2DT_" << jobID;

                    hoNDArray< std::complex<float> > res = gadget_->packages_received_[jobID].complexIm;
                    res.squeeze();
                    if ( !gadget_->debugFolder2_fullPath_.empty() ) { gadget_->gt_exporter_.exportArrayComplex(res, gadget_->debugFolder2_fullPath_+ostr.str()); }

                    if ( gadget_->packages_received_[jobID].complexImSecond.get_number_of_elements() > 0 )
                    {
                        std::ostringstream ostr;
                        ostr << "Recon2DT_Second_" << jobID;

                        hoNDArray< std::complex<float> > res = gadget_->packages_received_[jobID].complexImSecond;
                        res.squeeze();
                        if ( !gadget_->debugFolder2_fullPath_.empty() ) { gadget_->gt_exporter_.exportArrayComplex(res, gadget_->debugFolder2_fullPath_+ostr.str()); }

                        if ( gadget_->packages_received_[jobID].resTimeStampSecond.get_number_of_elements() > 0 )
                        {
                            std::ostringstream ostr;
                            ostr << "Recon2DT_Second_TimeStamp_" << jobID;

                            hoNDArray<float> res = gadget_->packages_received_[jobID].resTimeStampSecond;
                            res.squeeze();
                            if ( gadget_->debugFolder2_fullPath_.empty() ) { gadget_->gt_exporter_.exportArray(res, gadget_->debugFolder2_fullPath_+ostr.str()); }
                        }

                        if ( gadget_->packages_received_[jobID].resPhysioTimeStampSecond.get_number_of_elements() > 0 )
                        {
                            std::ostringstream ostr;
                            ostr << "Recon2DT_Second_PhysioTimeStamp_" << jobID;

                            hoNDArray<float> res = gadget_->packages_received_[jobID].resPhysioTimeStampSecond;
                            res.squeeze();
                            if ( gadget_->debugFolder2_fullPath_.empty() ) { gadget_->gt_exporter_.exportArray(res, gadget_->debugFolder2_fullPath_+ostr.str()); }
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GDEBUG("GtPlusRecon2DTGadgetCloudSender handling close...\n");
        return false;
    }

    return true;
}

}
