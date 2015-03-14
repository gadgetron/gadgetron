
#include "GtPlusRecon3DTGadget.h"

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

GtPlusRecon3DTGadget::GtPlusRecon3DTGadget() : BaseClass()
{

}

GtPlusRecon3DTGadget::~GtPlusRecon3DTGadget()
{

}

bool GtPlusRecon3DTGadget::readParameters()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(BaseClass::readParameters());

        GDEBUG_CONDITION_STREAM(verboseMode_, "------> GtPlusRecon3DTGadget parameters <------");

        boost::shared_ptr<std::string> str = this->get_string_value("dim_5th");
        para_.dim_5th_ = gtPlus_util_.getISMRMRDDimFromName(*str);
        GDEBUG_CONDITION_STREAM(verboseMode_, "dim_5th_ is " << *str);

        str = this->get_string_value("workOrder_ShareDim");
        para_.workOrder_ShareDim_ = gtPlus_util_.getISMRMRDDimFromName(*str);
        GDEBUG_CONDITION_STREAM(verboseMode_, "workOrder_ShareDim_ is " << *str);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.no_acceleration_averageall_ref_ = this->get_bool_value("no_acceleration_averageall_ref");
        GDEBUG_CONDITION_STREAM(verboseMode_, "no_acceleration_averageall_ref_ is " << para_.no_acceleration_averageall_ref_);

        para_.no_acceleration_same_combinationcoeff_allN_ = this->get_bool_value("no_acceleration_same_combinationcoeff_allN");
        GDEBUG_CONDITION_STREAM(verboseMode_, "no_acceleration_same_combinationcoeff_allN_ is " << para_.no_acceleration_same_combinationcoeff_allN_);

        para_.no_acceleration_whichN_combinationcoeff_ = this->get_int_value("no_acceleration_whichN_combinationcoeff");
        GDEBUG_CONDITION_STREAM(verboseMode_, "no_acceleration_whichN_combinationcoeff_ is " << para_.no_acceleration_whichN_combinationcoeff_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.interleaved_same_combinationcoeff_allN_ = this->get_bool_value("interleaved_same_combinationcoeff_allN");
        GDEBUG_CONDITION_STREAM(verboseMode_, "interleaved_same_combinationcoeff_allN_ is " << para_.interleaved_same_combinationcoeff_allN_);

        para_.interleaved_whichN_combinationcoeff_ = this->get_int_value("interleaved_whichN_combinationcoeff");
        GDEBUG_CONDITION_STREAM(verboseMode_, "interleaved_whichN_combinationcoeff_ is " << para_.interleaved_whichN_combinationcoeff_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.embedded_averageall_ref_ = this->get_bool_value("embedded_averageall_ref");
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_averageall_ref_ is " << para_.embedded_averageall_ref_);

        para_.embedded_fullres_coilmap_ = this->get_bool_value("embedded_fullres_coilmap");
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_fullres_coilmap_ is " << para_.embedded_fullres_coilmap_);

        para_.embedded_same_combinationcoeff_allN_ = this->get_bool_value("embedded_same_combinationcoeff_allN");
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_same_combinationcoeff_allN_ is " << para_.embedded_same_combinationcoeff_allN_);

        para_.embedded_whichN_combinationcoeff_ = this->get_int_value("embedded_whichN_combinationcoeff");
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_whichN_combinationcoeff_ is " << para_.embedded_whichN_combinationcoeff_);

        para_.embedded_ref_fillback_ = this->get_bool_value("embedded_ref_fillback");
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_ref_fillback_ is " << para_.embedded_ref_fillback_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.separate_averageall_ref_ = this->get_bool_value("separate_averageall_ref");
        GDEBUG_CONDITION_STREAM(verboseMode_, "separate_averageall_ref_ is " << para_.separate_averageall_ref_);

        para_.separate_fullres_coilmap_ = this->get_bool_value("separate_fullres_coilmap");
        GDEBUG_CONDITION_STREAM(verboseMode_, "separate_fullres_coilmap_ is " << para_.separate_fullres_coilmap_);

        para_.separate_same_combinationcoeff_allN_ = this->get_bool_value("separate_same_combinationcoeff_allN");
        GDEBUG_CONDITION_STREAM(verboseMode_, "separate_same_combinationcoeff_allN_ is " << para_.separate_same_combinationcoeff_allN_);

        para_.separate_whichN_combinationcoeff_ = this->get_int_value("separate_whichN_combinationcoeff");
        GDEBUG_CONDITION_STREAM(verboseMode_, "separate_whichN_combinationcoeff_ is " << para_.separate_whichN_combinationcoeff_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.same_coil_compression_coeff_allN_ = this->get_bool_value("same_coil_compression_coeff_allN");
        GDEBUG_CONDITION_STREAM(verboseMode_, "same_coil_compression_coeff_allN_ is " << para_.same_coil_compression_coeff_allN_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        // get the parameters from base class
        // BaseClass::readParameters();

        para_.recon_kspace_needed_ = recon_kspace_needed_;
        para_.workOrderPara_ = workOrderPara_;

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusRecon3DTGadget::readParameters() ... ");
        return false;
    }

    return true;
}

bool GtPlusRecon3DTGadget::setWorkOrder3DTParameters(WorkOrder3DTType* workOrder)
{
    workOrder->recon_kspace_needed_ = recon_kspace_needed_;

    if ( para_.workOrderPara_.coil_compression_thres_>0 || para_.workOrderPara_.coil_compression_num_modesKept_>0 )
    {
        workOrder->coil_compression_ = true;
    }
    else
    {
        workOrder->coil_compression_ = false;
    }

    workOrder->same_coil_compression_coeff_allN_ = para_.same_coil_compression_coeff_allN_;

    workOrder->embedded_averageall_ref_ = para_.embedded_averageall_ref_;
    workOrder->embedded_fullres_coilmap_ = para_.embedded_fullres_coilmap_;
    workOrder->embedded_same_combinationcoeff_allN_ = para_.embedded_same_combinationcoeff_allN_;
    workOrder->embedded_whichN_combinationcoeff_ = para_.embedded_whichN_combinationcoeff_;
    workOrder->embedded_ref_fillback_ = para_.embedded_ref_fillback_;

    workOrder->separate_averageall_ref_ = para_.separate_averageall_ref_;
    workOrder->separate_fullres_coilmap_ = para_.separate_fullres_coilmap_;
    workOrder->separate_same_combinationcoeff_allN_ = para_.separate_same_combinationcoeff_allN_;
    workOrder->separate_whichN_combinationcoeff_ = para_.separate_whichN_combinationcoeff_;

    //workOrder->interleaved_same_combinationcoeff_allN_ = interleaved_same_combinationcoeff_allN_;
    //workOrder->interleaved_whichN_combinationcoeff_ = interleaved_whichN_combinationcoeff_;

    workOrder->no_acceleration_averageall_ref_ = para_.no_acceleration_averageall_ref_;
    workOrder->no_acceleration_same_combinationcoeff_allN_ = para_.no_acceleration_same_combinationcoeff_allN_;
    workOrder->no_acceleration_whichN_combinationcoeff_ = para_.no_acceleration_whichN_combinationcoeff_;

    return true;
}

int GtPlusRecon3DTGadget::process_config(ACE_Message_Block* mb)
{
    // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
    //   0  1  2   3    4   5    6     7  8   9
    GADGET_CHECK_RETURN(BaseClass::process_config(mb)==GADGET_OK, GADGET_FAIL);

    // pre-allocate memory
    size_t numOfBytes;
    if ( para_.workOrderPara_.coil_compression_num_modesKept_ > 0 )
    {
        if ( num_acq_channels_ > 2*para_.workOrderPara_.coil_compression_num_modesKept_ )
        {
            numOfBytes = (size_t)( (double)matrix_size_encoding_[0]*kSpaceMaxAcqE1No_*kSpaceMaxAcqE2No_*num_acq_channels_*para_.workOrderPara_.coil_compression_num_modesKept_*sizeof(ValueType));
        }
        else
        {
            numOfBytes = (size_t)( (double)matrix_size_encoding_[0]*kSpaceMaxAcqE1No_*kSpaceMaxAcqE2No_*num_acq_channels_*para_.workOrderPara_.coil_compression_num_modesKept_*sizeof(ValueType) );
        }

        if ( para_.workOrderPara_.recon_algorithm_ == Gadgetron::ISMRMRD_GRAPPA && para_.workOrderPara_.job_num_of_N_>0 )
        {
            numOfBytes = (size_t)( (double)para_.workOrderPara_.job_num_of_N_*kSpaceMaxAcqE1No_*kSpaceMaxAcqE2No_*num_acq_channels_*para_.workOrderPara_.coil_compression_num_modesKept_*sizeof(ValueType)*1.5 );
        }
    }
    else
    {
        if ( para_.workOrderPara_.recon_algorithm_ == Gadgetron::ISMRMRD_SPIRIT || para_.workOrderPara_.recon_algorithm_ == Gadgetron::ISMRMRD_L1SPIRIT )
        {
            numOfBytes = (size_t)((double)matrix_size_encoding_[0]*kSpaceMaxAcqE1No_*kSpaceMaxAcqE2No_*num_acq_channels_*num_acq_channels_*sizeof(ValueType)*0.8);
        }
        else
        {
            numOfBytes = (size_t)((double)matrix_size_encoding_[0]*kSpaceMaxAcqE1No_*kSpaceMaxAcqE2No_*num_acq_channels_*num_acq_channels_*sizeof(ValueType)*0.6);
        }
    }

    if ( (num_acq_channels_<=12) || (para_.workOrderPara_.coil_compression_num_modesKept_>0 && 2*para_.workOrderPara_.coil_compression_num_modesKept_>num_acq_channels_) )
    {
        numOfBytes *= 2;
    }

    if ( numOfBytes > 1024*1024*1024*128.0 )
    {
        numOfBytes = (size_t)(1024*1024*1024*4.0);
    }

    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusRecon3DTGadget::Pre allocate : " << numOfBytes/1024.0/1024.0 << " Megabytes ... ");

    if ( performTiming_ ) { gt_timer1_.start("Pre-allocate memory ... "); }
    mem_manager_->increase(numOfBytes);
    if ( performTiming_ ) { gt_timer1_.stop(); }

    worker_grappa_.gtPlus_mem_manager_ = mem_manager_;
    worker_noacceleration_.gtPlus_mem_manager_ = mem_manager_;
    worker_spirit_.gtPlus_mem_manager_ = mem_manager_;
    worker_spirit_L1_ncg_.gtPlus_mem_manager_ = mem_manager_;

    if ( CloudComputing_ )
    {
        bool parseSuccess = this->parseGTCloudNodeFile(cloud_node_file_, gt_cloud_);
        if ( parseSuccess )
        {
            CloudSize_ = (unsigned int)gt_cloud_.size();
            if ( CloudSize_ == 0 ) CloudComputing_ = false;
        }
        else
        {
            CloudComputing_ = false;
        }
    }

    return GADGET_OK;
}

int GtPlusRecon3DTGadget::process(Gadgetron::GadgetContainerMessage< GtPlusGadgetImageArray >* m1, Gadgetron::GadgetContainerMessage< WorkOrderType > * m2)
{
    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusRecon3DTGadget::process(...) starts ... ");

    processed_called_times_++;

    GtPlusGadgetImageArray* images = m1->getObjectPtr();

    WorkOrderType* workOrder = m2->getObjectPtr();

    boost::shared_ptr< std::vector<size_t> > dims = workOrder->data_.get_dimensions();

    GDEBUG_CONDITION_STREAM(verboseMode_, "[Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave] = [" 
        << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " << (*dims)[3] << " " << (*dims)[4] 
        << " " << (*dims)[5] << " " << (*dims)[6] << " " << (*dims)[7] << " " << (*dims)[8] << " " << (*dims)[9] << " " << (*dims)[10] << "]");

    dimensions_ = *dims;

    // fill in more parameters
    para_.reconSizeRO_ = matrix_size_recon_[0];
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

    // ---------------------------------------------------------
    // set the work flow
    // ---------------------------------------------------------
    workflow_.reconSizeRO_ = para_.reconSizeRO_;
    workflow_.reconSizeE1_ = para_.reconSizeE1_;
    workflow_.reconSizeE2_ = para_.reconSizeE2_;
    workflow_.encodingFOV_RO_ = para_.encodingFOV_RO_;
    workflow_.encodingFOV_E1_ = para_.encodingFOV_E1_;
    workflow_.encodingFOV_E2_ = para_.encodingFOV_E2_;
    workflow_.reconFOV_RO_ = para_.reconFOV_RO_;
    workflow_.reconFOV_E1_ = para_.reconFOV_E1_;
    workflow_.reconFOV_E2_ = para_.reconFOV_E2_;

    workflow_.dataDimStartingIndexes_ = workOrder->dataDimStartingIndexes_;
    workflow_.dim5th_ = para_.dim_5th_;
    workflow_.WorkOrderShareDim_ = para_.workOrder_ShareDim_;
    workflow_.performTiming_ = performTiming_;

    // ---------------------------------------------------------
    // set work order
    // ---------------------------------------------------------
    workOrder->copyFromPara(para_.workOrderPara_);

    workOrder->CloudComputing_ = CloudComputing_;
    workOrder->CloudSize_ = CloudSize_;
    workOrder->gt_cloud_ = gt_cloud_;

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

    worker_spirit_L1_ncg_.verbose_ = verboseMode_;
    worker_spirit_L1_ncg_.performTiming_ = performTiming_;
    if ( !debugFolder_fullPath_.empty() ) worker_spirit_L1_ncg_.debugFolder_ = debugFolder_fullPath_;

    if ( !debugFolder_fullPath_.empty() ) workflow_.debugFolder_ = debugFolder_fullPath_;

    // if 'other' data is coming in
    if ( workOrder->other_.get_number_of_elements() > 0 )
    {
        workOrder->duplicate(workOrder_recon_other_);
        setWorkOrder3DTParameters(&workOrder_recon_other_);
        workflow_.workOrder_ = &workOrder_recon_other_;

        // perform a simple FFT recon
        workOrder_recon_other_.CalibMode_ = ISMRMRD_noacceleration;
        workOrder_recon_other_.acceFactorE1_ = 1;
        workOrder_recon_other_.acceFactorE2_ = 1;

        workOrder_recon_other_.start_RO_ = -1;
        workOrder_recon_other_.end_RO_ = -1;
        workOrder_recon_other_.start_E1_ = -1;
        workOrder_recon_other_.end_E1_ = -1;
        workOrder_recon_other_.start_E2_ = -1;
        workOrder_recon_other_.end_E2_ = -1;

        workflow_.worker_ = &worker_noacceleration_;
        workflow_.setDataArray(workOrder->other_);
        GADGET_CHECK_RETURN(workflow_.recon(), GADGET_FAIL);

        GADGET_CHECK_RETURN(this->scalingImages(workflow_.res_), GADGET_FAIL);
        GADGET_CHECK_RETURN(this->sendOutRecon(images, workflow_.res_, image_series_+1, workOrder->dataDimStartingIndexes_, "Other", GADGETRON_IMAGE_OTHER), GADGET_FAIL);

        workflow_.res_.clear();
        workflow_.data_ = NULL;
        workflow_.ref_ = NULL;
        workflow_.workOrder_ = NULL;

        workOrder_recon_other_.reset();
    }

    // perform the recon
    if ( performTiming_ ) { gt_timer1_.start("Recon 3DT workorder ... "); }

    GADGET_CHECK_RETURN(this->generateKSpaceFilter(*workOrder), GADGET_FAIL);

    workOrder->duplicate(workOrder_recon_);
    setWorkOrder3DTParameters(&workOrder_recon_);

    workflow_.workOrder_ = &workOrder_recon_;
    if ( verboseMode_ )
    {
        workflow_.workOrder_->print(std::cout);
    }

    workflow_.setDataArray(workOrder->data_, workOrder->time_stamp_, workOrder->physio_time_stamp_);

    if ( workOrder->ref_.get_number_of_elements() > 0 )
    {
        workflow_.setRefArray(workOrder->ref_);
    }
    else if ( CalibMode_==Gadgetron::ISMRMRD_interleaved )
    {
        workOrder->ref_ = workOrder->data_;
        workflow_.setRefArray(workOrder->ref_);
    }

    // set the work flow for worker and workOrder
    if ( workOrder->acceFactorE1_>1 || workOrder->acceFactorE2_>1 )
    {
        if ( para_.workOrderPara_.recon_algorithm_ == Gadgetron::ISMRMRD_SPIRIT )
        {
            workflow_.worker_ = &worker_spirit_;
        }
        else if ( para_.workOrderPara_.recon_algorithm_ == Gadgetron::ISMRMRD_L1SPIRIT )
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

    if ( workflow_.worker_ != &worker_grappa_ )
    {
        GWARN_STREAM("The gfactor computation is currently only avaialbe for grappa reconstruction ... ");
        workflow_.workOrder_->gfactor_needed_ = false;
    }

    GADGET_CHECK_RETURN(workflow_.preProcessing(), GADGET_FAIL);
    GADGET_CHECK_RETURN(workflow_.recon(), GADGET_FAIL);
    GADGET_CHECK_RETURN(workflow_.postProcessing(), GADGET_FAIL);

    if ( performTiming_ ) { gt_timer1_.stop(); }

    if ( !debugFolder2_fullPath_.empty() )
    {
        std::ostringstream ostr;
        ostr << "Recon3DT";

        hoNDArray< std::complex<float> > res = workflow_.res_;
        res.squeeze();
        if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder2_fullPath_+ostr.str()); }

        if ( workflow_.workOrder_->gfactor_needed_ )
        {
            std::ostringstream ostr;
            ostr << "Recon3DT_GFactor";

            hoNDArray< std::complex<float> > gfactor = workflow_.gfactor_;
            gfactor.squeeze();
            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(gfactor, debugFolder2_fullPath_+ostr.str()); }
        }
    }

    // compute SNR image and stdmap
    hoNDArray<ValueType> snrImage, stdMap;
    bool snrImageComputed = false;
    bool stdMapComputed = false;

    if ( workflow_.workOrder_->gfactor_needed_ || workOrder->acceFactorE1_*workOrder->acceFactorE2_==1 )
    {
        if ( scalingFactor_snr_image_>0 || scalingFactor_std_map_>0)
        {
            bool withAcceleration = (workOrder->acceFactorE1_*workOrder->acceFactorE2_>1);

            if ( !this->computeSNRImage(workflow_.res_, workflow_.gfactor_, 
                    start_frame_for_std_map_, withAcceleration, snrImage, stdMap) )
            {
                snrImage.clear();
                stdMap.clear();
            }
            else
            {
                snrImageComputed = true;
                stdMapComputed = true;
            }

            if ( workOrder->acceFactorE1_*workOrder->acceFactorE2_==1 ) snrImageComputed = false;
        }
    }

    // send out the results
    GADGET_CHECK_RETURN(this->scalingImages(workflow_.res_), GADGET_FAIL);
    GADGET_CHECK_RETURN(this->sendOutRecon(images, workflow_.res_, image_series_, workOrder->dataDimStartingIndexes_, "Image", GADGETRON_IMAGE_REGULAR), GADGET_FAIL);

    if ( workflow_.workOrder_->gfactor_needed_ )
    {
        Gadgetron::scal((float)scalingFactor_gfactor_, workflow_.gfactor_);
        GADGET_CHECK_RETURN(this->sendOutRecon(images, workflow_.gfactor_, image_series_+1, workOrder->dataDimStartingIndexes_, "gfactor", GADGETRON_IMAGE_GFACTOR), GADGET_FAIL);
    }

    if ( scalingFactor_snr_image_>0 && snrImage.get_number_of_elements()>0 && snrImageComputed )
    {
        Gadgetron::scal((float)scalingFactor_snr_image_, snrImage);
        GADGET_CHECK_RETURN(this->sendOutRecon(images, snrImage, image_series_+2, workOrder->dataDimStartingIndexes_, "snr_map", GADGETRON_IMAGE_SNR_MAP), GADGET_FAIL);
    }

    if ( scalingFactor_std_map_>0 && stdMap.get_number_of_elements()>0 && stdMapComputed )
    {
        Gadgetron::scal((float)scalingFactor_std_map_, stdMap);
        GADGET_CHECK_RETURN(this->sendOutRecon(images, stdMap, image_series_+3, workOrder->dataDimStartingIndexes_, "std_map", GADGETRON_IMAGE_STD_MAP), GADGET_FAIL);
    }

    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusRecon3DTGadget::process(...) ends ... ");

    // reset the status
    workflow_.data_ = NULL;
    workflow_.ref_ = NULL;
    workflow_.noise_ = NULL;
    workflow_.workOrder_ = NULL;
    // Gadgetron::clear(&workflow_.res_);

    m1->release();
    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GtPlusRecon3DTGadget)

}
