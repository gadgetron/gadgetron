/** \file   GenericReconAccumulateImageTriggerGadget.h
    \brief  This is the image accmulation and triggering gadget for the generic recon chain
            This gadget receives the IsmrmrdImageArray [RO E1 E2 CHA N S SLC] and can be triggered by CHA/SLC/CON/PHS/REP/SET/AVE

            For every dimension, user can define whether it is under the trigger. If the dimensional index of buffered images reache maximum for ALL 
            dimensions which are under the trigger, the image buffer will be send to the next gadget; e.g., if the PHS dimension limit is 40 and 
            the dimension PHS is under the trigger, all 40 images will be sent to the next gadget as a data buffer. Every buffered images will only be sent once
            Certain images, such as GADGETRON_IMAGE_GFACTOR gfactor images will be sent to the next gadget immediately

    \author Hui Xue
*/

#pragma once

#include <complex>
#include "gadgetron_mricore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"

#include "GadgetronTimer.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"

#include "mri_core_def.h"
#include "mri_core_data.h"
#include "mri_core_utility.h"

#include "hoNDArray_utils.h"
#include "hoMRImage.h"
#include "hoNDObjectArray.h"

#include "GenericReconBase.h"

#include "ImageIOAnalyze.h"

namespace Gadgetron { 

// the dimensionsal order of buffered images
// [Cha Slice Con Phase Rep Set Ave]
//   0    1    2   3    4   5    6
#define GT_DIM_NUM_IMAGE 7

template <typename T, int D> 
class EXPORTGADGETSMRICORE GenericReconAccumulateImageTriggerGadget : public Gadgetron::GenericReconBase< IsmrmrdImageArray >
{
public:
    GADGET_DECLARE(GenericReconAccumulateImageTriggerGadget);

    typedef std::complex<T> ValueType;

    typedef Gadgetron::GenericReconBase< IsmrmrdImageArray > BaseClass;

    typedef hoMRImage<ValueType, D> ImageType;
    typedef hoNDObjectArray<ImageType> ImageBufferType;
    typedef hoNDArray<bool> ImageSentFlagBufferType;

    GenericReconAccumulateImageTriggerGadget();
    ~GenericReconAccumulateImageTriggerGadget();

    virtual int close(unsigned long flags);

    /// ------------------------------------------------------------------------------------
    /// parameters to control the triggering
    /// ------------------------------------------------------------------------------------

    GADGET_PROPERTY(TriggerChannel, bool, "Whether to trigger along channel", false);
    GADGET_PROPERTY(TriggerSlice, bool, "Whether to trigger along slice", false);
    GADGET_PROPERTY(TriggerContrast, bool, "Whether to trigger along contrast", false);
    GADGET_PROPERTY(TriggerPhase, bool, "Whether to trigger along phase", false);
    GADGET_PROPERTY(TriggerRepetition, bool, "Whether to trigger along repetition", false);
    GADGET_PROPERTY(TriggerSet, bool, "Whether to trigger along set", false);
    GADGET_PROPERTY(TriggerAverage, bool, "Whether to trigger along average", false);
    GADGET_PROPERTY(PassImageImmediately, bool, "Whether to pass image immediately", false);

    // allow user to input image boundary limit
    // if <=0, then the input values are ignored and protocol settings are used
    GADGET_PROPERTY(meas_max_kspace_encode_step_1, int, "Maximal encoding limit for E1", 0);
    GADGET_PROPERTY(meas_max_set, int, "Maximal encoding limit for SET", 0);
    GADGET_PROPERTY(meas_max_phase, int, "Maximal encoding limit for PHS", 0);
    GADGET_PROPERTY(meas_max_kspace_encode_step_2, int, "Maximal encoding limit for E2", 0);
    GADGET_PROPERTY(meas_max_contrast, int, "Maximal encoding limit for CON", 0);
    GADGET_PROPERTY(meas_max_slice, int, "Maximal encoding limit for SLC", 0);
    GADGET_PROPERTY(meas_max_repetition, int, "Maximal encoding limit for REP", 0);
    GADGET_PROPERTY(meas_max_average, int, "Maximal encoding limit for AVE", 0);

    // whether to consider concatenation on repetition
    GADGET_PROPERTY(concatenation_on_repetition, bool, "If multiple concatenation is used, whether to enlarge the repetition limit", false);

protected:

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(GadgetContainerMessage<IsmrmrdImageArray>* m1);

    // perform the triggering
    virtual int trigger(ImageBufferType& buf, ImageSentFlagBufferType& sentFlagBuf, bool inClose);

    // store the incoming image
    // if pass_image_immediate_==true, the image will be immediately passed to the next gadget with 
    virtual int store_image(const IsmrmrdImageArray& img, const std::vector<size_t>& buf_dimension, ImageBufferType& buf);

    // set dimensions under trigger
    void set_dimensions_under_trigger();

    /// the dimension limits by default is read from the protocol,but user overwrites the protocol setting via the input parameters
    ISMRMRD::EncodingCounters meas_max_idx_;

    /// whether a dimension is under the trigger
    /// if no dimension is under the trigger, images will be passed to next gadget right away
    bool cha_trigger_;
    bool slc_trigger_;
    bool con_trigger_;
    bool phs_trigger_;
    bool rep_trigger_;
    bool set_trigger_;
    bool ave_trigger_;

    /// whether to immediately pass the image to the next gadget
    bool pass_image_immediate_;

    // buffer for regular images whose data role is GADGETRON_IMAGE_REGULAR
    ImageBufferType imageBuffer_;
    ImageSentFlagBufferType imageSent_;

    // buffer for other images whole data role is not GADGETRON_IMAGE_REGULAR and not GADGETRON_IMAGE_GFACTOR
    ImageBufferType otherBuffer_;
    ImageSentFlagBufferType otherSent_;

    // buffer sent to next gadget
    ImageBufferType imageSentBuffer_;

    // image waveform
    std::vector<ISMRMRD::Waveform> wave_form_buffer_;

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

    int image_counter_;
};

class EXPORTGADGETSMRICORE GenericReconAccumulateImage2DTriggerGadget : public GenericReconAccumulateImageTriggerGadget<float, 2>
{
public:
    GADGET_DECLARE(GenericReconAccumulateImage2DTriggerGadget);

    typedef GenericReconAccumulateImageTriggerGadget<float, 2> BaseClass;

    GenericReconAccumulateImage2DTriggerGadget();
    ~GenericReconAccumulateImage2DTriggerGadget();
};

class EXPORTGADGETSMRICORE GenericReconAccumulateImage3DTriggerGadget : public GenericReconAccumulateImageTriggerGadget<float, 3>
{
public:
    GADGET_DECLARE(GenericReconAccumulateImage3DTriggerGadget);

    typedef GenericReconAccumulateImageTriggerGadget<float, 3> BaseClass;

    GenericReconAccumulateImage3DTriggerGadget();
    ~GenericReconAccumulateImage3DTriggerGadget();
};

}
