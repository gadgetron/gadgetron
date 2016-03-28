
#include "GtPlusRecon2DTGadget.h"

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

GtPlusRecon2DTGadget::GtPlusRecon2DTGadget() : BaseClass()
{

}

GtPlusRecon2DTGadget::~GtPlusRecon2DTGadget()
{

}

bool GtPlusRecon2DTGadget::readParameters()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(BaseClass::readParameters());

        GDEBUG_CONDITION_STREAM(verboseMode_, "------> GtPlusRecon2DTGadget parameters <------");

        para_.dim_4th_ = gtPlus_util_.getISMRMRDDimFromName(dim_4th.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "dim_4th_ is " << dim_4th.value());

        para_.dim_5th_ = gtPlus_util_.getISMRMRDDimFromName(dim_5th.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "dim_5th_ is " << dim_5th.value());

        para_.workOrder_ShareDim_ = gtPlus_util_.getISMRMRDDimFromName(workOrder_ShareDim.value());
        GDEBUG_CONDITION_STREAM(verboseMode_, "workOrder_ShareDim_ is " << workOrder_ShareDim.value());

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.no_acceleration_averageall_ref_ = no_acceleration_averageall_ref.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "no_acceleration_averageall_ref_ is " << para_.no_acceleration_averageall_ref_);

        para_.no_acceleration_ref_numOfModes_ = no_acceleration_ref_numOfModes.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "no_acceleration_ref_numOfModes_ is " << para_.no_acceleration_ref_numOfModes_);

        para_.no_acceleration_same_combinationcoeff_allS_ = no_acceleration_same_combinationcoeff_allS.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "no_acceleration_same_combinationcoeff_allS_ is " << para_.no_acceleration_same_combinationcoeff_allS_);

        para_.no_acceleration_whichS_combinationcoeff_ = no_acceleration_whichS_combinationcoeff.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "no_acceleration_whichS_combinationcoeff_ is " << para_.no_acceleration_whichS_combinationcoeff_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.interleaved_same_combinationcoeff_allS_ = interleaved_same_combinationcoeff_allS.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "interleaved_same_combinationcoeff_allS_ is " << para_.interleaved_same_combinationcoeff_allS_);

        para_.interleaved_ref_numOfModes_ = interleaved_ref_numOfModes.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "interleaved_ref_numOfModes_ is " << para_.interleaved_ref_numOfModes_);

        para_.interleaved_whichS_combinationcoeff_ = interleaved_whichS_combinationcoeff.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "interleaved_whichS_combinationcoeff_ is " << para_.interleaved_whichS_combinationcoeff_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.embedded_averageall_ref_ = embedded_averageall_ref.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_averageall_ref_ is " << para_.embedded_averageall_ref_);

        para_.embedded_ref_numOfModes_ = embedded_ref_numOfModes.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_ref_numOfModes_ is " << para_.embedded_ref_numOfModes_);

        para_.embedded_fullres_coilmap_ = embedded_fullres_coilmap.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_fullres_coilmap_ is " << para_.embedded_fullres_coilmap_);

        para_.embedded_fullres_coilmap_useHighestSignal_ = embedded_fullres_coilmap_useHighestSignal.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_fullres_coilmap_useHighestSignal_ is " << para_.embedded_fullres_coilmap_useHighestSignal_);

        para_.embedded_same_combinationcoeff_allS_ = embedded_same_combinationcoeff_allS.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_same_combinationcoeff_allS_ is " << para_.embedded_same_combinationcoeff_allS_);

        para_.embedded_whichS_combinationcoeff_ = embedded_whichS_combinationcoeff.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_whichS_combinationcoeff_ is " << para_.embedded_whichS_combinationcoeff_);

        para_.embedded_ref_fillback_ = embedded_ref_fillback.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "embedded_ref_fillback_ is " << para_.embedded_ref_fillback_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.separate_averageall_ref_ = separate_averageall_ref.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "separate_averageall_ref_ is " << para_.separate_averageall_ref_);

        para_.separate_ref_numOfModes_ = separate_ref_numOfModes.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "separate_ref_numOfModes_ is " << para_.separate_ref_numOfModes_);

        para_.separate_fullres_coilmap_ = separate_fullres_coilmap.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "separate_fullres_coilmap_ is " << para_.separate_fullres_coilmap_);

        para_.separate_same_combinationcoeff_allS_ = separate_same_combinationcoeff_allS.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "separate_same_combinationcoeff_allS_ is " << para_.separate_same_combinationcoeff_allS_);

        para_.separate_whichS_combinationcoeff_ = separate_whichS_combinationcoeff.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "separate_whichS_combinationcoeff_ is " << para_.separate_whichS_combinationcoeff_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        para_.same_coil_compression_coeff_allS_ = same_coil_compression_coeff_allS.value();
        GDEBUG_CONDITION_STREAM(verboseMode_, "same_coil_compression_coeff_allS_ is " << para_.same_coil_compression_coeff_allS_);

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");

        // get the parameters from base class
        // BaseClass::readParameters();

        para_.recon_kspace_needed_ = recon_kspace_needed_;
        para_.workOrderPara_ = workOrderPara_;

        GDEBUG_CONDITION_STREAM(verboseMode_, "-----------------------------------------------");
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlusRecon2DTGadget::readParameters() ... ");
        return false;
    }

    return true;
}

bool GtPlusRecon2DTGadget::setWorkOrder2DTParameters(WorkOrder2DTType* workOrder)
{
    workOrder->recon_kspace_needed_ = para_.recon_kspace_needed_;

    if ( para_.workOrderPara_.coil_compression_thres_>0 || para_.workOrderPara_.coil_compression_num_modesKept_>0 )
    {
        workOrder->coil_compression_ = true;
    }
    else
    {
        workOrder->coil_compression_ = false;
    }

    workOrder->same_coil_compression_coeff_allS_ = para_.same_coil_compression_coeff_allS_;

    workOrder->embedded_averageall_ref_ = para_.embedded_averageall_ref_;
    workOrder->embedded_ref_numOfModes_ = para_.embedded_ref_numOfModes_;
    workOrder->embedded_fullres_coilmap_ = para_.embedded_fullres_coilmap_;
    workOrder->embedded_fullres_coilmap_useHighestSignal_ = para_.embedded_fullres_coilmap_useHighestSignal_;
    workOrder->embedded_same_combinationcoeff_allS_ = para_.embedded_same_combinationcoeff_allS_;
    workOrder->embedded_whichS_combinationcoeff_ = para_.embedded_whichS_combinationcoeff_;
    workOrder->embedded_ref_fillback_ = para_.embedded_ref_fillback_;

    workOrder->separate_averageall_ref_ = para_.separate_averageall_ref_;
    workOrder->separate_ref_numOfModes_ = para_.separate_ref_numOfModes_;
    workOrder->separate_fullres_coilmap_ = para_.separate_fullres_coilmap_;
    workOrder->separate_same_combinationcoeff_allS_ = para_.separate_same_combinationcoeff_allS_;
    workOrder->separate_whichS_combinationcoeff_ = para_.separate_whichS_combinationcoeff_;

    workOrder->interleaved_same_combinationcoeff_allS_ = para_.interleaved_same_combinationcoeff_allS_;
    workOrder->interleaved_whichS_combinationcoeff_ = para_.interleaved_whichS_combinationcoeff_;
    workOrder->interleaved_ref_numOfModes_ = para_.interleaved_ref_numOfModes_;

    workOrder->no_acceleration_averageall_ref_ = para_.no_acceleration_averageall_ref_;
    workOrder->no_acceleration_ref_numOfModes_ = para_.no_acceleration_ref_numOfModes_;
    workOrder->no_acceleration_same_combinationcoeff_allS_ = para_.no_acceleration_same_combinationcoeff_allS_;
    workOrder->no_acceleration_whichS_combinationcoeff_ = para_.no_acceleration_whichS_combinationcoeff_;

    return true;
}

int GtPlusRecon2DTGadget::process_config(ACE_Message_Block* mb)
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
        else
        {
            CloudComputing_ = false;
        }
    }

    return GADGET_OK;
}

int GtPlusRecon2DTGadget::process(Gadgetron::GadgetContainerMessage< GtPlusGadgetImageArray >* m1, Gadgetron::GadgetContainerMessage< WorkOrderType > * m2)
{
    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusRecon2DTGadget::process(...) starts ... ");

    processed_called_times_++;


    GtPlusGadgetImageArray* images = m1->getObjectPtr();

    WorkOrderType* workOrder = m2->getObjectPtr();

    boost::shared_ptr< std::vector<size_t> > dims = workOrder->data_.get_dimensions();

    size_t SEG = (*dims)[9];

    GDEBUG_CONDITION_STREAM(verboseMode_, "[Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave] = [" 
                                                << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " 
                                                << (*dims)[3] << " " << (*dims)[4] << " " << (*dims)[5] << " " 
                                                << (*dims)[6] << " " << (*dims)[7] << " " << (*dims)[8] << " " 
                                                << (*dims)[9] << " " << (*dims)[10] << "]");

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

    para_.workOrderPara_.retro_gated_images_ = workOrder->retro_gated_images_;
    para_.workOrderPara_.retro_gated_segment_size_ = workOrder->retro_gated_segment_size_;

    para_.workOrderPara_.workFlow_BufferKernel_ = workOrder->workFlow_BufferKernel_;
    para_.workOrderPara_.workFlow_use_BufferedKernel_ = workOrder->workFlow_use_BufferedKernel_;
    para_.workOrderPara_.num_channels_res_ = workOrder->num_channels_res_;

    bool perform_retro_gating = (para_.workOrderPara_.retro_gated_images_>0);

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
    workflow_.dim4th_ = para_.dim_4th_;
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

    if ( !debugFolder_fullPath_.empty() ) workflow_.debugFolder_ = debugFolder_fullPath_;

    // if 'other' data is coming in
    if ( workOrder->other_.get_number_of_elements() > 0 )
    {
        workOrder->duplicate(workOrder_recon_other_);
        setWorkOrder2DTParameters(&workOrder_recon_other_);
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

        //hoNDArray<ValueType> resResized;
        //GADGET_CHECK_RETURN(gtPlus_util_complex_.zpadResize2D(workflow_.res_, workflow_.reconSizeRO_, workflow_.reconSizeE1_, resResized), GADGET_FAIL);
        //GADGET_CHECK_RETURN(this->sendOutRecon(images, resResized, image_series_+1, workOrder->dataDimStartingIndexes_, "Other"), GADGET_FAIL);

        GADGET_CHECK_RETURN(this->scalingImages(workflow_.res_), GADGET_FAIL);
        GADGET_CHECK_RETURN(this->sendOutRecon(images, workflow_.res_, image_series_+1, workOrder->dataDimStartingIndexes_, "Other", GADGETRON_IMAGE_OTHER), GADGET_FAIL);

        workflow_.res_.clear();
        workflow_.data_ = NULL;
        workflow_.ref_ = NULL;
        workflow_.workOrder_ = NULL;

        workOrder_recon_other_.reset();
    }

    // ------------------------------------------------------------------
    // perform the recon
    // ------------------------------------------------------------------
    if ( performTiming_ ) { gt_timer1_.start("Recon 2DT workorder ..."); }

    GADGET_CHECK_RETURN(this->generateKSpaceFilter(*workOrder), GADGET_FAIL);

    /// set the work order
    workOrder->duplicate(workOrder_recon_);
    setWorkOrder2DTParameters(&workOrder_recon_);

    workflow_.workOrder_ = &workOrder_recon_;

    if ( verboseMode_ )
    {
        workflow_.workOrder_->print(std::cout);
    }

    /// set the data
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
    if ( workOrder->acceFactorE1_ > 1 )
    {
        if ( para_.workOrderPara_.recon_algorithm_ == Gadgetron::ISMRMRD_SPIRIT )
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

    if ( workflow_.worker_ != &worker_grappa_ )
    {
        GWARN_STREAM("The gfactor computation is currently only avaialbe for grappa reconstruction ... ");
        workflow_.workOrder_->gfactor_needed_ = false;

        GWARN_STREAM("The wrap-around map computation is currently only avaialbe for grappa reconstruction ... ");
        workflow_.workOrder_->wrap_around_map_needed_ = false;
    }

    /// perform the recon
    GADGET_CHECK_RETURN(workflow_.preProcessing(), GADGET_FAIL);
    GADGET_CHECK_RETURN(workflow_.recon(), GADGET_FAIL);
    GADGET_CHECK_RETURN(workflow_.postProcessing(), GADGET_FAIL);

    if ( performTiming_ ) { gt_timer1_.stop(); }

    if ( !debugFolder2_fullPath_.empty() )
    {
        std::ostringstream ostr;
        ostr << "Recon2DT_" << processed_called_times_;

        hoNDArray< std::complex<float> > res = workflow_.res_;
        res.squeeze();
        if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder2_fullPath_+ostr.str()); }

        if ( workflow_.workOrder_->gfactor_needed_ )
        {
            std::ostringstream ostr;
            ostr << "Recon2DT_GFactor_" << processed_called_times_;

            hoNDArray< std::complex<float> > gfactor = workflow_.gfactor_;
            gfactor.squeeze();
            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(gfactor, debugFolder2_fullPath_+ostr.str()); }
        }

        if ( workflow_.workOrder_->wrap_around_map_needed_ )
        {
            std::ostringstream ostr;
            ostr << "Recon2DT_WrapAroundMap_" << processed_called_times_;

            hoNDArray< std::complex<float> > wrap_around_map = workflow_.wrap_around_map_;
            wrap_around_map.squeeze();
            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(wrap_around_map, debugFolder2_fullPath_+ostr.str()); }
        }

        if ( workflow_.res_second_.get_number_of_elements() > 0 )
        {
            hoNDArray< std::complex<float> > res = workflow_.res_second_;
            res.squeeze();

            std::ostringstream ostr;
            ostr << "Recon2DT_second_" << processed_called_times_;

            if ( !debugFolder2_fullPath_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder2_fullPath_+ostr.str()); }
        }
    }

    // compute SNR image and stdmap
    hoNDArray<ValueType> snrImage, stdMap;
    bool snrImageComputed = false;
    bool stdMapComputed = false;

    if ( workflow_.workOrder_->gfactor_needed_ || workOrder->acceFactorE1_==1 )
    {
        if ( scalingFactor_snr_image_>0 || scalingFactor_std_map_>0)
        {
            bool withAcceleration = (workOrder->acceFactorE1_>1);

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

            if ( workOrder->acceFactorE1_==1 ) snrImageComputed = false;
        }
    }

    // send out the results
    GADGET_CHECK_RETURN(this->scalingImages(workflow_.res_), GADGET_FAIL);

    if ( send_out_recon_ )
    {
        if ( perform_retro_gating )
        {
            GADGET_CHECK_RETURN(this->sendOutRecon(images, workflow_.res_, workflow_.res_time_stamp_, workflow_.res_physio_time_stamp_, image_series_, workOrder->dataDimStartingIndexes_, "ImageRetro", GADGETRON_IMAGE_RETRO), GADGET_FAIL);
        }
        else
        {
            GADGET_CHECK_RETURN(this->sendOutRecon(images, workflow_.res_, workflow_.res_time_stamp_, workflow_.res_physio_time_stamp_, image_series_, workOrder->dataDimStartingIndexes_, "Image", GADGETRON_IMAGE_REGULAR), GADGET_FAIL);
        }

        if ( workflow_.workOrder_->gfactor_needed_ )
        {
            Gadgetron::scal((float)scalingFactor_gfactor_, workflow_.gfactor_);
            GADGET_CHECK_RETURN(this->sendOutRecon(images, workflow_.gfactor_, workflow_.res_time_stamp_, workflow_.res_physio_time_stamp_, image_series_+1, workOrder->dataDimStartingIndexes_, "gfactor", GADGETRON_IMAGE_GFACTOR), GADGET_FAIL);
        }

        if ( workflow_.workOrder_->wrap_around_map_needed_ )
        {
            Gadgetron::scal((float)scalingFactor_wrap_around_map_, workflow_.wrap_around_map_);
            GADGET_CHECK_RETURN(this->sendOutRecon(images, workflow_.wrap_around_map_, workflow_.res_time_stamp_, workflow_.res_physio_time_stamp_, image_series_+2, workOrder->dataDimStartingIndexes_, "wrap_around_map", GADGETRON_IMAGE_WRAPAROUNDMAP), GADGET_FAIL);
        }

        if ( scalingFactor_snr_image_>0 && snrImage.get_number_of_elements()>0 && snrImageComputed )
        {
            Gadgetron::scal((float)scalingFactor_snr_image_, snrImage);
            GADGET_CHECK_RETURN(this->sendOutRecon(images, snrImage, workflow_.res_time_stamp_, workflow_.res_physio_time_stamp_, image_series_+3, workOrder->dataDimStartingIndexes_, "snr_map", GADGETRON_IMAGE_SNR_MAP), GADGET_FAIL);
        }

        if ( scalingFactor_std_map_>0 && stdMap.get_number_of_elements()>0 && stdMapComputed )
        {
            Gadgetron::scal((float)scalingFactor_std_map_, stdMap);
            GADGET_CHECK_RETURN(this->sendOutRecon(images, stdMap, workflow_.res_time_stamp_, workflow_.res_physio_time_stamp_, image_series_+4, workOrder->dataDimStartingIndexes_, "std_map", GADGETRON_IMAGE_STD_MAP), GADGET_FAIL);
        }
    }

    if ( send_out_recon_second_ )
    {
        if ( workflow_.res_second_.get_number_of_elements() > 0 )
        {
            Gadgetron::scal((float)scalingFactor_, workflow_.res_second_);
            GADGET_CHECK_RETURN(this->sendOutRecon(images, workflow_.res_second_, workflow_.res_time_stamp_second_, workflow_.res_physio_time_stamp_second_, image_series_+5, workOrder->dataDimStartingIndexes_, "ImageSecond", GADGETRON_IMAGE_REGULAR), GADGET_FAIL);
        }
    }

    GDEBUG_CONDITION_STREAM(verboseMode_, "GtPlusRecon2DTGadget::process(...) ends ... ");

    // reset the status
    workflow_.data_ = NULL;
    workflow_.time_stamp_ = NULL;
    workflow_.physio_time_stamp_ = NULL;
    workflow_.ref_ = NULL;
    workflow_.noise_ = NULL;
    workflow_.workOrder_ = NULL;
    // Gadgetron::clear(&workflow_.res_);

    m1->release();
    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GtPlusRecon2DTGadget)

}
