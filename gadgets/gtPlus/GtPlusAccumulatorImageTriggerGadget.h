/** \file   GtPlusAccumulatorImageTriggerGadget.h
    \brief  The GtPlus image accmulation and triggering gadget, used after GtPlus reconstruction for image data
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "hoNDObjectArray.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/meta.h"
#include "GadgetIsmrmrdReadWrite.h"

#include "hoNDArray_utils.h"
#include "hoNDImage.h"

#include "GtPlusGadgetImageArray.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"

namespace Gadgetron
{

// the dimensionsal order of buffered images
// [Cha Slice E2 Con Phase Rep Set Ave]
//   0    1    2   3   4    5   6   7
#define GT_DIM_NUM_IMAGE 8

class EXPORTGTPLUSGADGET GtPlusAccumulatorImageTriggerGadget : public Gadget3< ISMRMRD::ImageHeader, hoNDArray< std::complex<float> >, ISMRMRD::MetaContainer >
{
public:
    GADGET_DECLARE(GtPlusAccumulatorImageTriggerGadget);

    typedef std::complex<float> ValueType;

    typedef Gadget3< ISMRMRD::ImageHeader, hoNDArray< ValueType >, ISMRMRD::MetaContainer > BaseClass;

    typedef hoNDImage<ValueType, 2> ImageType;

    typedef hoNDObjectArray<ImageType> ImageBufferType;
    typedef hoNDArray<bool> ImageSentFlagBufferType;

    GtPlusAccumulatorImageTriggerGadget();
    ~GtPlusAccumulatorImageTriggerGadget();

    virtual int close(unsigned long flags);

    /// parameters to control the triggering

    /// for every dimension, user can define whether it is under the trigger
    /// if the dimensional index of buffered images reache maximum for all dimensions under the trigger, 
    /// the image buffer will be send to the next gadget
    /// e.g., if the PHS dimension limit is 40 and the dimension PHS is under the trigger, all 40 images 
    /// will be sent to the next gadget as a data buffer
    /// every buffered images will only  be sent once
    /// GADGETRON_IMAGE_GFACTOR gfactor images will be sent to the next gadget immediately

    /// dimension limits
    /// the dimension limits by default is read from the protocol,but 
    /// user can set them via the input parameters
    ISMRMRD::EncodingCounters meas_max_idx_;

    /// whether a dimension is under the trigger
    /// if no dimension is under the trigger, images will be passed to next gadget right away
    bool cha_trigger_;
    bool slc_trigger_;
    bool e2_trigger_;
    bool con_trigger_;
    bool phs_trigger_;
    bool rep_trigger_;
    bool set_trigger_;
    bool ave_trigger_;

    /// whether to immediately pass the image to the next gadget
    bool pass_image_immediate_;

protected:

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray<ValueType> >* m2, GadgetContainerMessage<ISMRMRD::MetaContainer>* m3);

    // perform the triggering
    virtual bool trigger(ImageBufferType& buf, ImageSentFlagBufferType& sentFlagBuf, bool inClose);

    // store the incoming image
    // if pass_image_immediate_==true, the image will be immediately passed to the next gadget with 
    virtual bool storeImage(const ISMRMRD::ImageHeader& imgHeader, const hoNDArray<ValueType>& img, const ISMRMRD::MetaContainer& attrib, ImageBufferType& buf);

    // set dimensions under trigger
    void setDimensionsUnderTrigger();

    // buffer for regular images whose data role is GADGETRON_IMAGE_REGULAR
    ImageBufferType imageBuffer_;
    ImageSentFlagBufferType imageSent_;

    // buffer for other images whole data role is not GADGETRON_IMAGE_REGULAR and not GADGETRON_IMAGE_GFACTOR
    ImageBufferType otherBuffer_;
    ImageSentFlagBufferType otherSent_;

    // buffer sent to next gadget
    ImageBufferType imageSentBuffer_;

    // number of total dimensions
    size_t num_of_dimensions_;

    // dimensions under trigger
    std::vector<bool> dim_under_trigger_;
    std::vector<size_t> dim_limit_under_trigger_;

    std::vector<bool> dim_not_under_trigger_;
    std::vector<size_t> dim_limit_not_under_trigger_;

    // whether the next gadget has been triggered in close(...)
    bool triggered_in_close_;

    // dimension for image kspace
    std::vector<size_t> dimensions_;

    // encoding matrix size (the real sampled size)
    size_t matrix_size_encoding_[3];

    // encoding space size (the logic kspace size)
    size_t space_size_[3];

    // encoding filed of view [mm]
    float field_of_view_encoding_[3];

    // recon matrix size (the final image size)
    size_t matrix_size_recon_[3];

    // recon filed of view [mm]
    float field_of_view_recon_[3];

    int image_counter_;

    int meas_max_ro_;
    int meas_max_channel_;

    // util for gtplus
    Gadgetron::gtPlus::gtPlusISMRMRDReconUtil< std::complex<float> > gtPlus_util_;

    // in verbose mode, more info is printed out
    bool verboseMode_;
};

}
