
#include "GtPlusAccumulatorImageTriggerGadget.h"
#include "GtPlusReconGadgetUtil.h"

using namespace Gadgetron::gtPlus;

namespace Gadgetron
{

GtPlusAccumulatorImageTriggerGadget::GtPlusAccumulatorImageTriggerGadget() : image_counter_(0), triggered_in_close_(false), verboseMode_(false)
{
    cha_trigger_ = false;
    slc_trigger_ = false;
    e2_trigger_ = false;
    con_trigger_ = false;
    phs_trigger_ = false;
    rep_trigger_ = false;
    set_trigger_ = false;
    ave_trigger_ = false;

    num_of_dimensions_ = 8; // [CHA SLC E2 CON PHS REP SET AVE]

    // this may be changed later if multi-channel image workflow are used
    meas_max_channel_ = 1;

    pass_image_immediate_ = false;
}

GtPlusAccumulatorImageTriggerGadget::~GtPlusAccumulatorImageTriggerGadget()
{

}

// extract necessary configuration information from the xml
int GtPlusAccumulatorImageTriggerGadget::process_config(ACE_Message_Block* mb)
{
    // gadget parameters
    verboseMode_ = this->get_bool_value("verboseMode");

    cha_trigger_ = this->get_bool_value("TriggerChannel");
    slc_trigger_ = this->get_bool_value("TriggerSlice");
    e2_trigger_  = this->get_bool_value("TriggerE2");
    con_trigger_ = this->get_bool_value("TriggerContrast");
    phs_trigger_ = this->get_bool_value("TriggerPhase");
    rep_trigger_ = this->get_bool_value("TriggerRepetition");
    set_trigger_ = this->get_bool_value("TriggerSet");
    ave_trigger_ = this->get_bool_value("TriggerAverage");

    pass_image_immediate_ = this->get_bool_value("PassImageImmediately");

    // ---------------------------------------------------------------------------------------------------------
    // pass the xml file
    ISMRMRD::IsmrmrdHeader h;
    try {
      deserialize(mb->rd_ptr(),h);
    } catch (...) {
      GDEBUG("Error parsing ISMRMRD Header");
      throw;
      return GADGET_FAIL;
    }

    // seq object
    if (h.encoding.size() != 1)
    {
        GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
        GDEBUG("This simple GtPlusAccumulatorImageTriggerGadget only supports one encoding space\n");
        return GADGET_FAIL;
    }

    // ---------------------------------------------------------------------------------------------------------

    // find out the encoding space 
    findMatrixSizeEncoding(h, matrix_size_encoding_);
    findFOVEncoding(h, field_of_view_encoding_);

    findMatrixSizeRecon(h, matrix_size_recon_);
    findFOVRecon(h, field_of_view_recon_);

    GDEBUG_CONDITION_STREAM(verboseMode_, "Encoding matrix size: " << matrix_size_encoding_[0] << " " << matrix_size_encoding_[1] << " " << matrix_size_encoding_[2]);
    GDEBUG_CONDITION_STREAM(verboseMode_, "Encoding field_of_view : " << field_of_view_encoding_[0] << " " << field_of_view_encoding_[1] << " " << field_of_view_encoding_[2]);
    GDEBUG_CONDITION_STREAM(verboseMode_, "Recon matrix size : " << matrix_size_recon_[0] << " " << matrix_size_recon_[1] << " " << matrix_size_recon_[2]);
    GDEBUG_CONDITION_STREAM(verboseMode_, "Recon field_of_view :  " << field_of_view_recon_[0] << " " << field_of_view_recon_[1] << " " << field_of_view_recon_[2]);

    // ---------------------------------------------------------------------------------------------------------
    // encoding limits
    GADGET_CHECK_RETURN(findEncodingLimits(h, meas_max_idx_, verboseMode_), GADGET_FAIL);

    //ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

    //if (e_limits.kspace_encoding_step_1) 
    //{
    //    meas_max_idx_.kspace_encode_step_1 = (uint16_t)(matrix_size_encoding_[1]-1); // e_limits.kspace_encoding_step_1().get().maximum();
    //}
    //else
    //{
    //    meas_max_idx_.kspace_encode_step_1 = 0;
    //    GDEBUG_STREAM("Setting number of kspace_encode_step_1 to 0" << std::endl);
    //    return GADGET_FAIL;
    //}

    //if (e_limits.set)
    //{
    //    if ( e_limits.set->maximum > 0 )
    //        meas_max_idx_.set = e_limits.set->maximum - 1;
    //    else
    //        meas_max_idx_.set = 0;

    //    if ( meas_max_idx_.set < 0 ) meas_max_idx_.set = 0;
    //}
    //else
    //{
    //    meas_max_idx_.set = 0;
    //}

    //if (e_limits.phase)
    //{
    //    if ( e_limits.phase->maximum > 0 )
    //        meas_max_idx_.phase = e_limits.phase->maximum-1;
    //    else
    //        meas_max_idx_.phase = 0;

    //    if ( meas_max_idx_.phase < 0 ) meas_max_idx_.phase = 0;
    //}
    //else
    //{
    //    meas_max_idx_.phase = 0;
    //}

    //if (e_limits.kspace_encoding_step_2)
    //{
    //    meas_max_idx_.kspace_encode_step_2 = (uint16_t)(matrix_size_encoding_[2] - 1); // e_limits.kspace_encoding_step_2().get().maximum();
    //}
    //else
    //{
    //    meas_max_idx_.kspace_encode_step_2 = 0;
    //}
    //meas_max_idx_.kspace_encode_step_2 = (uint16_t)(matrix_size_recon_[2]);

    //if (e_limits.contrast)
    //{
    //    if ( e_limits.contrast->maximum > 0 )
    //        meas_max_idx_.contrast = e_limits.contrast->maximum-1;
    //    else
    //        meas_max_idx_.contrast = 0;

    //    if ( meas_max_idx_.contrast < 0 ) meas_max_idx_.contrast = 0;
    //}
    //else
    //{
    //    meas_max_idx_.contrast = 0;
    //}

    //if (e_limits.slice)
    //{
    //    meas_max_idx_.slice = e_limits.slice->maximum;
    //}
    //else
    //{
    //    meas_max_idx_.slice = 0;
    //}

    //if (e_limits.repetition)
    //{
    //    meas_max_idx_.repetition = e_limits.repetition->maximum;
    //}
    //else
    //{
    //    meas_max_idx_.repetition = 0;
    //}

    //if (e_limits.average)
    //{
    //    meas_max_idx_.average = e_limits.average->maximum-1;
    //}
    //else
    //{
    //    meas_max_idx_.average = 0;
    //}

    //if (e_limits.segment)
    //{
    //    // meas_max_idx_.segment = e_limits.segment().get().maximum()-1;
    //    meas_max_idx_.segment = 0;
    //}
    //else
    //{
    //    meas_max_idx_.segment = 0;
    //}

    // allocate the image buffers
    // [Cha Slice E2 Con Phase Rep Set Ave]
    //   0    1    2   3   4    5   6   7

    meas_max_idx_.kspace_encode_step_2 = (uint16_t)(matrix_size_recon_[2]);

    dimensions_.resize(GT_DIM_NUM_IMAGE, 0);
    dimensions_[0] = meas_max_channel_;
    dimensions_[1] = meas_max_idx_.slice+1;
    dimensions_[2] = meas_max_idx_.kspace_encode_step_2;
    dimensions_[3] = meas_max_idx_.contrast+1;
    dimensions_[4] = meas_max_idx_.phase+1;
    dimensions_[5] = meas_max_idx_.repetition+1;
    dimensions_[6] = meas_max_idx_.set+1;
    dimensions_[7] = meas_max_idx_.average+1;

    imageBuffer_.create(dimensions_);
    imageSent_.create(dimensions_);

    otherBuffer_.create(dimensions_);
    otherSent_.create(dimensions_);

    size_t nElem = imageBuffer_.get_number_of_elements();
    size_t ii;
    for ( ii=0; ii<nElem; ii++ )
    {
        imageBuffer_(ii) = NULL;
        otherBuffer_(ii) = NULL;
        imageSent_(ii) = false;
        otherSent_(ii) = false;
    }

    // set the dimensions under/not under trigger
    this->setDimensionsUnderTrigger();

    GDEBUG_CONDITION_STREAM(verboseMode_, "dimension limits                [Cha Slice E2 Con Phase Rep Set Ave] = [" 
                               << " " << dimensions_[0] 
                               << " " << dimensions_[1] 
                               << " " << dimensions_[2] 
                               << " " << dimensions_[3]
                               << " " << dimensions_[4]
                               << " " << dimensions_[5]
                               << " " << dimensions_[6] 
                               << " " << dimensions_[7] << "]");

    GDEBUG_CONDITION_STREAM(verboseMode_, "dimension under trigger         [Cha Slice E2 Con Phase Rep Set Ave] = [" 
                               << " " << dim_under_trigger_[0] 
                               << " " << dim_under_trigger_[1] 
                               << " " << dim_under_trigger_[2] 
                               << " " << dim_under_trigger_[3]
                               << " " << dim_under_trigger_[4]
                               << " " << dim_under_trigger_[5]
                               << " " << dim_under_trigger_[6] 
                               << " " << dim_under_trigger_[7] << "]");

    GDEBUG_CONDITION_STREAM(verboseMode_, "dimension limits under trigger  [Cha Slice E2 Con Phase Rep Set Ave] = [" 
                               << " " << dim_limit_under_trigger_[0] 
                               << " " << dim_limit_under_trigger_[1] 
                               << " " << dim_limit_under_trigger_[2] 
                               << " " << dim_limit_under_trigger_[3]
                               << " " << dim_limit_under_trigger_[4]
                               << " " << dim_limit_under_trigger_[5]
                               << " " << dim_limit_under_trigger_[6] 
                               << " " << dim_limit_under_trigger_[7] << "]");

    GDEBUG_CONDITION_STREAM(verboseMode_, "dimension NOT under trigger     [Cha Slice E2 Con Phase Rep Set Ave] = [" 
                               << " " << dim_not_under_trigger_[0] 
                               << " " << dim_not_under_trigger_[1] 
                               << " " << dim_not_under_trigger_[2] 
                               << " " << dim_not_under_trigger_[3]
                               << " " << dim_not_under_trigger_[4]
                               << " " << dim_not_under_trigger_[5]
                               << " " << dim_not_under_trigger_[6] 
                               << " " << dim_not_under_trigger_[7] << "]");

    GDEBUG_CONDITION_STREAM(verboseMode_, "dimension limits NOT under trigger [Cha Slice E2 Con Phase Rep Set Ave] = [" 
                               << " " << dim_limit_not_under_trigger_[0] 
                               << " " << dim_limit_not_under_trigger_[1] 
                               << " " << dim_limit_not_under_trigger_[2] 
                               << " " << dim_limit_not_under_trigger_[3]
                               << " " << dim_limit_not_under_trigger_[4]
                               << " " << dim_limit_not_under_trigger_[5]
                               << " " << dim_limit_not_under_trigger_[6] 
                               << " " << dim_limit_not_under_trigger_[7] << "]");

    return GADGET_OK;
}

void GtPlusAccumulatorImageTriggerGadget::setDimensionsUnderTrigger()
{
    dim_under_trigger_.resize(num_of_dimensions_, false);
    dim_not_under_trigger_.resize(num_of_dimensions_, false);

    dim_limit_under_trigger_.resize(num_of_dimensions_, 1);
    dim_limit_not_under_trigger_.resize(num_of_dimensions_, 1);

    if (cha_trigger_)
    {
        dim_under_trigger_[0] = true;
        dim_limit_under_trigger_[0] = dimensions_[0];
    }
    else
    {
        dim_not_under_trigger_[0] = true;
        dim_limit_not_under_trigger_[0] = dimensions_[0];
    }

    if (slc_trigger_)
    {
        dim_under_trigger_[1] = true;
        dim_limit_under_trigger_[1] = dimensions_[1];
    }
    else
    {
        dim_not_under_trigger_[1] = true;
        dim_limit_not_under_trigger_[1] = dimensions_[1];
    }

    if (e2_trigger_)
    {
        dim_under_trigger_[2] = true;
        dim_limit_under_trigger_[2] = dimensions_[2];
    }
    else
    {
        dim_not_under_trigger_[2] = true;
        dim_limit_not_under_trigger_[2] = dimensions_[2];
    }

    if (con_trigger_)
    {
        dim_under_trigger_[3] = true;
        dim_limit_under_trigger_[3] = dimensions_[3];
    }
    else
    {
        dim_not_under_trigger_[3] = true;
        dim_limit_not_under_trigger_[3] = dimensions_[3];
    }

    if (phs_trigger_)
    {
        dim_under_trigger_[4] = true;
        dim_limit_under_trigger_[4] = dimensions_[4];
    }
    else
    {
        dim_not_under_trigger_[4] = true;
        dim_limit_not_under_trigger_[4] = dimensions_[4];
    }

    if (rep_trigger_)
    {
        dim_under_trigger_[5] = true;
        dim_limit_under_trigger_[5] = dimensions_[5];
    }
    else
    {
        dim_not_under_trigger_[5] = true;
        dim_limit_not_under_trigger_[5] = dimensions_[5];
    }

    if (set_trigger_)
    {
        dim_under_trigger_[6] = true;
        dim_limit_under_trigger_[6] = dimensions_[6];
    }
    else
    {
        dim_not_under_trigger_[6] = true;
        dim_limit_not_under_trigger_[6] = dimensions_[6];
    }

    if (ave_trigger_)
    {
        dim_under_trigger_[7] = true;
        dim_limit_under_trigger_[7] = dimensions_[7];
    }
    else
    {
        dim_not_under_trigger_[7] = true;
        dim_limit_not_under_trigger_[7] = dimensions_[7];
    }

    imageSentBuffer_.create(dim_limit_under_trigger_);
    imageSentBuffer_.delete_data_on_destruct(false);
}

int GtPlusAccumulatorImageTriggerGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3)
{
    // find the data role
    std::string dataRole;
    dataRole = std::string(m3->getObjectPtr()->as_str(GTPLUS_DATA_ROLE, 0));

    GDEBUG_CONDITION_STREAM(verboseMode_, "--> receive image : " << m1->getObjectPtr()->image_index << " -- " << dataRole);

    if ( dataRole == GTPLUS_IMAGE_REGULAR )
    {
        GADGET_CHECK_RETURN(this->storeImage(*m1->getObjectPtr(), *m2->getObjectPtr(), *m3->getObjectPtr(), imageBuffer_), GADGET_FAIL);
        GADGET_CHECK_RETURN(this->trigger(imageBuffer_, imageSent_, false), GADGET_FAIL);
    }

    if ( dataRole == GTPLUS_IMAGE_OTHER )
    {
        GADGET_CHECK_RETURN(this->storeImage(*m1->getObjectPtr(), *m2->getObjectPtr(), *m3->getObjectPtr(), otherBuffer_), GADGET_FAIL);
        GADGET_CHECK_RETURN(this->trigger(otherBuffer_, otherSent_, false), GADGET_FAIL);
    }

    if ( dataRole == GTPLUS_IMAGE_GFACTOR )
    {
        // pass the image to the next gadget
        Gadgetron::GadgetContainerMessage<ImageBufferType>* cm1 = new Gadgetron::GadgetContainerMessage<ImageBufferType>();

        ImageBufferType& imgBuf = *(cm1->getObjectPtr());

        std::vector<size_t> dim2D(num_of_dimensions_, 1);
        imgBuf.create(dim2D);
        imgBuf(0) = new ImageType();
        GADGET_CHECK_RETURN(imgBuf(0)!=NULL, GADGET_FAIL);

        // set image content
        imgBuf(0)->from_NDArray( *m2->getObjectPtr() );
        // set image attrib
        imgBuf(0)->attrib_ = *m3->getObjectPtr();

        // pass the ISMRMRD header info
        GADGET_CHECK_RETURN(gtPlus_util_.setMetaAttributesFromImageHeaderISMRMRD(*m1->getObjectPtr(), imgBuf(0)->attrib_), GADGET_FAIL);

        if (this->next()->putq(cm1) < 0) 
        {
            m1->release();
            cm1->release();
            return GADGET_FAIL;
        }
    }

    m1->release();
    return GADGET_OK;
}

bool GtPlusAccumulatorImageTriggerGadget::trigger(ImageBufferType& buf, ImageSentFlagBufferType& sentFlagBuf, bool inClose)
{
    try
    {
        // scan the buffered images, if the trigger dimensions are complete, sent out this package

        // not under trigger
        size_t cha, slc, e2, con, phs, rep, set, ave;

        // under trigger
        size_t cha_t, slc_t, e2_t, con_t, phs_t, rep_t, set_t, ave_t;

        std::vector<size_t> image_ind(num_of_dimensions_, 0);
        std::vector<size_t> image_sent_ind(num_of_dimensions_, 0);

        size_t numOfElem = imageSentBuffer_.get_number_of_elements();
        size_t ii;
        for ( ii=0; ii<numOfElem; ii++ ) { imageSentBuffer_(ii) = NULL; }

        for ( ave=0; ave<dim_limit_not_under_trigger_[7]; ave++ )
        {
            if ( dim_not_under_trigger_[7] ) image_ind[7] = ave;
            // -------------------
            for ( set=0; set<dim_limit_not_under_trigger_[6]; set++ )
            {
                if ( dim_not_under_trigger_[6] ) image_ind[6] = set;
                // -------------------
                for ( rep=0; rep<dim_limit_not_under_trigger_[5]; rep++ )
                {
                    if ( dim_not_under_trigger_[5] ) image_ind[5] = rep;
                    // -------------------
                    for ( phs=0; phs<dim_limit_not_under_trigger_[4]; phs++ )
                    {
                        if ( dim_not_under_trigger_[4] ) image_ind[4] = phs;
                        // -------------------
                        for ( con=0; con<dim_limit_not_under_trigger_[3]; con++ )
                        {
                            if ( dim_not_under_trigger_[3] ) image_ind[3] = con;
                            // -------------------
                            for ( e2=0; e2<dim_limit_not_under_trigger_[2]; e2++ )
                            {
                                if ( dim_not_under_trigger_[2] ) image_ind[2] = e2;
                                // -------------------
                                for ( slc=0; slc<dim_limit_not_under_trigger_[1]; slc++ )
                                {
                                    if ( dim_not_under_trigger_[1] ) image_ind[1] = slc;
                                    // -------------------
                                    for ( cha=0; cha<dim_limit_not_under_trigger_[0]; cha++ )
                                    {
                                        if ( dim_not_under_trigger_[0] ) image_ind[0] = cha;
                                        // -------------------

                                        // loop over under triggered dimensions and check whether every images are there
                                        bool needTrigger = true;
                                        if ( inClose )
                                        {
                                            needTrigger = false;
                                        }

                                        {
                                            for ( ii=0; ii<numOfElem; ii++ ) { imageSentBuffer_(ii) = NULL; }

                                            // =================================================

                                            for ( ave_t=0; ave_t<dim_limit_under_trigger_[7]; ave_t++ )
                                            {
                                                if ( dim_under_trigger_[7] ) image_ind[7] = ave_t;
                                                image_sent_ind[7] = ave_t;
                                                // -------------------
                                                for ( set_t=0; set_t<dim_limit_under_trigger_[6]; set_t++ )
                                                {
                                                    if ( dim_under_trigger_[6] ) image_ind[6] = set_t;
                                                    image_sent_ind[6] = set_t;
                                                    // -------------------
                                                    for ( rep_t=0; rep_t<dim_limit_under_trigger_[5]; rep_t++ )
                                                    {
                                                        if ( dim_under_trigger_[5] ) image_ind[5] = rep_t;
                                                        image_sent_ind[5] = rep_t;
                                                        // -------------------
                                                        for ( phs_t=0; phs_t<dim_limit_under_trigger_[4]; phs_t++ )
                                                        {
                                                            if ( dim_under_trigger_[4] ) image_ind[4] = phs_t;
                                                            image_sent_ind[4] = phs_t;
                                                            // -------------------
                                                            for ( con_t=0; con_t<dim_limit_under_trigger_[3]; con_t++ )
                                                            {
                                                                if ( dim_under_trigger_[3] ) image_ind[3] = con_t;
                                                                image_sent_ind[3] = con_t;
                                                                // -------------------
                                                                for ( e2_t=0; e2_t<dim_limit_under_trigger_[2]; e2_t++ )
                                                                {
                                                                    if ( dim_under_trigger_[2] ) image_ind[2] = e2_t;
                                                                    image_sent_ind[2] = e2_t;
                                                                    // -------------------
                                                                    for ( slc_t=0; slc_t<dim_limit_under_trigger_[1]; slc_t++ )
                                                                    {
                                                                        if ( dim_under_trigger_[1] ) image_ind[1] = slc_t;
                                                                        image_sent_ind[1] = slc_t;
                                                                        // -------------------
                                                                        for ( cha_t=0; cha_t<dim_limit_under_trigger_[0]; cha_t++ )
                                                                        {
                                                                            if ( dim_under_trigger_[0] ) image_ind[0] = cha_t;
                                                                            image_sent_ind[0] = cha_t;
                                                                            // -------------------

                                                                            ImageType* pImage = buf(image_ind);
                                                                            bool sentFlag = sentFlagBuf(image_ind);

                                                                            if ( inClose )
                                                                            {
                                                                                // if in close call, send out all unsent images
                                                                                if ( pImage != NULL && !sentFlag )
                                                                                {
                                                                                    imageSentBuffer_(image_sent_ind) = pImage;
                                                                                    buf(image_ind) = NULL;
                                                                                    needTrigger = true;
                                                                                }
                                                                            }
                                                                            else
                                                                            {
                                                                                if ( pImage != NULL && !sentFlag )
                                                                                {
                                                                                    imageSentBuffer_(image_sent_ind) = pImage;
                                                                                    // buf(image_ind) = NULL;
                                                                                }
                                                                                else
                                                                                {
                                                                                    needTrigger = false; // if all images for current under-trigger dimensions are filled, trigger
                                                                                    break;
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }

                                            if ( needTrigger )
                                            {
                                                // if a image has been sent, not sent again
                                                for ( ave_t=0; ave_t<dim_limit_under_trigger_[7]; ave_t++ )
                                                {
                                                    if ( dim_under_trigger_[7] ) image_ind[7] = ave_t;
                                                    for ( set_t=0; set_t<dim_limit_under_trigger_[6]; set_t++ )
                                                    {
                                                        if ( dim_under_trigger_[6] ) image_ind[6] = set_t;
                                                        for ( rep_t=0; rep_t<dim_limit_under_trigger_[5]; rep_t++ )
                                                        {
                                                            if ( dim_under_trigger_[5] ) image_ind[5] = rep_t;
                                                            for ( phs_t=0; phs_t<dim_limit_under_trigger_[4]; phs_t++ )
                                                            {
                                                                if ( dim_under_trigger_[4] ) image_ind[4] = phs_t;
                                                                for ( con_t=0; con_t<dim_limit_under_trigger_[3]; con_t++ )
                                                                {
                                                                    if ( dim_under_trigger_[3] ) image_ind[3] = con_t;
                                                                    for ( e2_t=0; e2_t<dim_limit_under_trigger_[2]; e2_t++ )
                                                                    {
                                                                        if ( dim_under_trigger_[2] ) image_ind[2] = e2_t;
                                                                        for ( slc_t=0; slc_t<dim_limit_under_trigger_[1]; slc_t++ )
                                                                        {
                                                                            if ( dim_under_trigger_[1] ) image_ind[1] = slc_t;
                                                                            for ( cha_t=0; cha_t<dim_limit_under_trigger_[0]; cha_t++ )
                                                                            {
                                                                                if ( dim_under_trigger_[0] ) image_ind[0] = cha_t;

                                                                                bool sentFlag = sentFlagBuf(image_ind);
                                                                                if ( sentFlag )
                                                                                {
                                                                                    imageSentBuffer_(cha_t, slc_t, e2_t, con_t, phs_t, rep_t, set_t) = NULL;
                                                                                }
                                                                                else
                                                                                {
                                                                                    sentFlagBuf(image_ind) = true;
                                                                                }

                                                                                buf(image_ind) = NULL;
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }

                                                GDEBUG_STREAM("--> Accumulator image trigger for [CHA SLC E2 CON PHS REP SET AVE] : [" 
                                                                                                                            << image_ind[0] << " " 
                                                                                                                            << image_ind[1] << " " 
                                                                                                                            << image_ind[2] << " " 
                                                                                                                            << image_ind[3] << " " 
                                                                                                                            << image_ind[4] << " " 
                                                                                                                            << image_ind[5] << " " 
                                                                                                                            << image_ind[6] << " " 
                                                                                                                            << image_ind[7] << "]" );

                                                Gadgetron::GadgetContainerMessage<ImageBufferType>* cm1 = new Gadgetron::GadgetContainerMessage<ImageBufferType>();
                                                ImageBufferType& imgBuf = *(cm1->getObjectPtr());
                                                imgBuf = imageSentBuffer_;
                                                imgBuf.delete_data_on_destruct(true);

                                                if (this->next()->putq(cm1) < 0) 
                                                {
                                                    cm1->release();
                                                    return false;
                                                }
                                            }
                                            else
                                            {
                                                for ( ii=0; ii<numOfElem; ii++ )
                                                {
                                                    imageSentBuffer_(ii) = NULL;
                                                }
                                            }

                                            // =================================================
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happens in GtPlusAccumulatorImageTriggerGadget::trigger(...) ... ");
        return false;
    }

    return true;
}

bool GtPlusAccumulatorImageTriggerGadget::storeImage(const ISMRMRD::ImageHeader& imgHeader, const hoNDArray<ValueType>& img, const ISMRMRD::MetaContainer& attrib, ImageBufferType& buf)
{
    try
    {
        long long cha = attrib.as_long(GTPLUS_CHA, 0);

        size_t slc = imgHeader.slice;

        long long e2 = attrib.as_long(GTPLUS_E2, 0);

        size_t con = imgHeader.contrast;
        size_t phs = imgHeader.phase;
        size_t rep = imgHeader.repetition;
        size_t set = imgHeader.set;
        size_t ave = imgHeader.average;

        // create image
        ImageType* storedImage = new ImageType();
        GADGET_CHECK_RETURN_FALSE(storedImage!=NULL);

        storedImage->from_NDArray(img);
        storedImage->attrib_ = attrib;
        GADGET_CHECK_RETURN_FALSE(gtPlus_util_.setMetaAttributesFromImageHeaderISMRMRD(imgHeader, storedImage->attrib_));

        storedImage->attrib_.set(GTPLUS_PASS_IMMEDIATE, (long)0);
        buf(cha, slc, e2, con, phs, rep, set, ave) = storedImage;

        if ( pass_image_immediate_ )
        {
            Gadgetron::GadgetContainerMessage<ImageBufferType>* cm1 = new Gadgetron::GadgetContainerMessage<ImageBufferType>();

            ImageBufferType& imgBuf = *(cm1->getObjectPtr());

            std::vector<size_t> dim2D(num_of_dimensions_, 1);
            imgBuf.create(dim2D);

            imgBuf(0) = new ImageType();
            *imgBuf(0) = *storedImage;

            // set the pass_image flag, so next gadget knows
            imgBuf(0)->attrib_.set(GTPLUS_PASS_IMMEDIATE, (long)1);

            if (this->next()->putq(cm1) < 0) 
            {
                cm1->release();
                return false;
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happens in GtPlusAccumulatorImageTriggerGadget::storeImage(const ISMRMRD::ImageHeader& imgHeader, const hoNDArray<ValueType>& img, const ISMRMRD::MetaContainer& attrib, ImageBufferType& buf) ... ");
        return false;
    }

    return true;
}

int GtPlusAccumulatorImageTriggerGadget::close(unsigned long flags)
{
    GDEBUG_CONDITION_STREAM(true, "GtPlusAccumulatorImageTriggerGadget - close(flags) : " << flags);

    if ( BaseClass::close(flags) != GADGET_OK ) return GADGET_FAIL;

    if ( flags!=0 && !triggered_in_close_ )
    {
        triggered_in_close_ = true;

        GDEBUG_CONDITION_STREAM(true, "GtPlusAccumulatorImageTriggerGadget - trigger in close(flags) ... ");

        GADGET_CHECK_RETURN(this->trigger(imageBuffer_, imageSent_, true), GADGET_FAIL);
        GADGET_CHECK_RETURN(this->trigger(otherBuffer_, otherSent_, true), GADGET_FAIL);
    }

    return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GtPlusAccumulatorImageTriggerGadget)

}
