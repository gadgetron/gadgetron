/** \file   GtPlusAccumulatorWorkOrderTriggerGadget.h
    \brief  The GtPlus reconstruction data accmulation and triggering gadget
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"

#include "hoNDArray_utils.h"

#include "GtPlusGadgetImageArray.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorkOrder.h"

namespace Gadgetron
{

// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg AVE]
//   0  1  2   3    4   5    6     7  8   9  10

struct ReadOutBuffer
{
    ISMRMRD::AcquisitionHeader acqHead_;
    hoNDArray< std::complex<float> > data_;
    bool isReflect_;
};

class EXPORTGTPLUSGADGET GtPlusAccumulatorWorkOrderTriggerGadget : public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
{
public:
    GADGET_DECLARE(GtPlusAccumulatorWorkOrderTriggerGadget);

    typedef float real_value_type;
    typedef std::complex<real_value_type> ValueType;

    typedef Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< ValueType > > BaseClass;

    typedef std::vector< ReadOutBuffer > ReadOutBufferType;
    typedef hoNDArray< std::complex<float> > BufferType;
    typedef hoNDArray< int > TimeStampBufferType;
    typedef hoNDArray< unsigned short > ReflectBufferType;

    // typedef Gadgetron::gtPlus::gtPlusReconWorkOrder2DT<ValueType> WorkOrderType;
    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder<ValueType> WorkOrderType;

    GtPlusAccumulatorWorkOrderTriggerGadget();
    ~GtPlusAccumulatorWorkOrderTriggerGadget();

    virtual int close(unsigned long flags);

    /// parameters to control the triggering

    GADGET_PROPERTY(image_series, int, "Image series number", 100);
    GADGET_PROPERTY(verboseMode, bool, "Whether to print more information", false);
    GADGET_PROPERTY(timeStampResolution, float, "Time tick resolution in second", 0.0025f);

    // ----------------------------------------------

    GADGET_PROPERTY_LIMITS(noacceleration_triggerDim1, std::string, "First trigger dimension for no acceleration mode", "DIM_NONE", 
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2", 
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY_LIMITS(noacceleration_triggerDim2, std::string, "Second trigger dimension for no acceleration mode", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY(noacceleration_numOfKSpace_triggerDim1, int, "Number of full kspace for dynamic triggering, no acceleration mode", 0);

    // ----------------------------------------------

    GADGET_PROPERTY_LIMITS(interleaved_triggerDim1, std::string, "First trigger dimension for interleaved mode", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY_LIMITS(interleaved_triggerDim2, std::string, "Second trigger dimension for interleaved mode", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY(interleaved_numOfKSpace_triggerDim1, int, "Number of full kspace for dynamic triggering, interleaved mode", 0);

    // ----------------------------------------------

    GADGET_PROPERTY_LIMITS(embedded_triggerDim1, std::string, "First trigger dimension for embedded mode", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY_LIMITS(embedded_triggerDim2, std::string, "Second trigger dimension for embedded mode", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY(embedded_numOfKSpace_triggerDim1, int, "Number of kspace for dynamic triggering, embedded mode", 0);

    // ----------------------------------------------

    GADGET_PROPERTY_LIMITS(separate_triggerDim1, std::string, "First trigger dimension for separate mode", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY_LIMITS(separate_triggerDim2, std::string, "Second trigger dimension for separate mode", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");

    GADGET_PROPERTY(separate_numOfKSpace_triggerDim1, int, "Number of kspace for dynamic triggering, separate mode", 0);

    // ----------------------------------------------

    GADGET_PROPERTY_LIMITS(other_kspace_matching_Dim, std::string, "Extra kspace dimension to differentiate other readout lines", "DIM_NONE",
        GadgetPropertyLimitsEnumeration, "DIM_NONE", "DIM_ReadOut", "DIM_Encoding1", "DIM_Channel", "DIM_Slice", "DIM_Encoding2",
        "DIM_Contrast", "DIM_Phase", "DIM_Repetition", "DIM_Set", "DIM_Segment", "DIM_Average", "DIM_other1", "DIM_other2", "DIM_other3");


    /// for interleaved mode
    // if DIM_NONE, the trigger is performed in the close function

    // the triggering for interleaved mode is defined as:
    // a) if triggerDim1_==DIM_NONE, and triggerDim2_ != DIM_NONE
    // the trigger is performed whenever the triggerDim2_ is changed

    // b) if triggerDim2_==NONE and triggerDim1_ != NONE, 
    // if numOfKSpace_triggerDim1_==0, 
    // the trigger is performed whenever the triggerDim1_ is changed

    // if numOfKSpace_triggerDim1_>0, 
    // the trigger is first performed when numOfKSpace_triggerDim1_ of triggerDim1_ kspace is buffered, 
    // then trigger is performed whenever a new triggerDim1_ kspace arrives
    // the new triggerDim1_ kspace will be reconed using the kernel estimated from numOfKSpace_triggerDim1_ of triggerDim1_

    // when the resetTriggerStatus(m1)==true , the status is resetted and the numOfKSpace_triggerDim1_ of triggerDim1_ kspace 
    // will be buffered and then trigger the recon

    // c) if both triggerDim1_ and triggerDim2_ are NONE, the trigger is performed 
    // in the close(flags) functions

    // e) if both triggerDim1_ and triggerDim2_ are NOT NONE, the trigger is first performed
    // when numOfKSpace_triggerDim1_ of triggerDim1_ kspace is buffered and then trigger is performed whenever a new 
    // triggerDim1_ kspace arrives
    // the new triggerDim1_ kspace will be reconed using the kernel estimated from numOfKSpace_triggerDim1_ of triggerDim1_
    // when the triggerDim2_ changes or resetTriggerStatus(m1)==true , the status is resetted and the numOfKSpace_triggerDim1_ of triggerDim1_ kspace 
    // will be buffered and then trigger the recon

    // f) if numOfKSpace_triggerDim1_==0 and both triggerDim1_ and triggerDim2_ are NOT NONE,
    // the trigger is performed whenever the triggerDim2_ is changed

    // noacceleration
    Gadgetron::ISMRMRDDIM noacceleration_triggerDim1_;
    Gadgetron::ISMRMRDDIM noacceleration_triggerDim2_;
    int noacceleration_numOfKSpace_triggerDim1_;

    // interleaved
    Gadgetron::ISMRMRDDIM interleaved_triggerDim1_;
    Gadgetron::ISMRMRDDIM interleaved_triggerDim2_;
    int interleaved_numOfKSpace_triggerDim1_;

    // embedded
    Gadgetron::ISMRMRDDIM embedded_triggerDim1_;
    Gadgetron::ISMRMRDDIM embedded_triggerDim2_;
    int embedded_numOfKSpace_triggerDim1_;

    // separate
    Gadgetron::ISMRMRDDIM separate_triggerDim1_;
    Gadgetron::ISMRMRDDIM separate_triggerDim2_;
    int separate_numOfKSpace_triggerDim1_;

    // for other kspace data, if other_kspace_matching_Dim != DIM_NONE, the other data dimension will be made to match the image data at 
    // dimension other_kspace_matching_Dim
    Gadgetron::ISMRMRDDIM other_kspace_matching_Dim_;

    // default behavior is to compare the readout geometry
    // if the imaging slice changes, the trigger status is reset
    virtual bool resetTriggerStatus(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1);

protected:

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1, GadgetContainerMessage< hoNDArray< ValueType > > * m2);

    // check the status of incoming readout
    // bIsKSpace: whether this data is for image
    // bIsRef: whether this data is for calibration signal
    // bIsNoise: whether this data is a noise scan
    // bIsPhaseCorr: whether this data is for phase correction
    // bIsReflect: whether this data is acquired reflectly (for EPI and similar scans)
    // bIsOther: other scans
    virtual bool checkStatus(uint64_t flag, int samples, 
                        bool& bIsKSpace, bool& bIsRef, bool& bIsNoise, bool& bIsPhaseCorr, bool& bIsReflect, bool& bIsOther,
                        bool& bIsNavigator, bool& bIsRTFeedback, bool& bIsHPFeedback, bool& bIsDummyScan);

    // store the image data
    virtual bool storeImageData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, bool isReflect);

    // store the ref data
    virtual bool storeRefData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, bool isReflect);

    // fill the dynamically buffered data
    virtual bool fillBuffer(ReadOutBufferType& readOutBuffer, BufferType& buf, ReflectBufferType& reflectBuf);

    // fill the per 2D image info
    virtual bool fillImageInfo(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GtPlusGadgetImageArray* messageImage, const ISMRMRD::EncodingCounters& idx);

    // compute the encoded size
    size_t computeEncodedSizeE1(size_t centerE1, size_t maxE1);
    size_t computeEncodedSizeE2(size_t centerE2, size_t maxE2);

    // perform the triggering
    virtual bool triggerWorkOrder(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, bool inClose, bool isKSpace);

    // workFlow_BufferKernel_ and workFlow_use_BufferedKernel_ is the command to work flow
    // if workFlow_BufferKernel_ == true, work flow will buffer the kernels computed for this work order
    // if workFlow_use_BufferedKernel_ == true, work flow will recon this work order using buffered kernels
    // if both triggerDim1_ and triggerDim2_ are NOT NONE and numOfKSpace_triggerDim1_ > 0, 
    // the first work order with workFlow_BufferKernel_==true will be sent out when values of triggerDim1_  equals numOfKSpace_triggerDim1_-1
    // the next work orders will be sent out when triggerDim1_ changes with workFlow_BufferKernel_==false and workFlow_use_BufferedKernel_==true
    // when the triggerDim2_ changes, the status will be reset
    virtual bool triggerWorkOrder(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, 
            Gadgetron::ISMRMRDDIM& triggerDim1_, Gadgetron::ISMRMRDDIM& triggerDim2_, int numOfKSpace_triggerDim1_);

    // trigger by extract array with triggerDim being value
    virtual bool triggerByDimEqual(Gadgetron::ISMRMRDDIM& triggerDim, size_t value, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_);
    virtual bool triggerByDimLessEqual(Gadgetron::ISMRMRDDIM& triggerDim, size_t value, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_);
    virtual bool triggerByDimEqual(Gadgetron::ISMRMRDDIM& triggerDim1, size_t value1, Gadgetron::ISMRMRDDIM& triggerDim2, size_t value2, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_);

    // trigger by extract array with triggerDim being <= value
    virtual bool triggerByDim1LessEqualDim2Equal(Gadgetron::ISMRMRDDIM& triggerDim1, size_t value1, Gadgetron::ISMRMRDDIM& triggerDim2, size_t value2, bool workFlow_BufferKernel_, bool workFlow_use_BufferedKernel_);

    // whether need to trigger all buffered data in close()
    bool needTriggerWorkOrderAllInClose();
    // trigger with all buffered data
    virtual bool triggerWorkOrderAllInClose();

    // trigger the last count in the close function
    virtual bool triggerWorkOrderLastCountInClose(Gadgetron::ISMRMRDDIM& triggerDim1_, Gadgetron::ISMRMRDDIM& triggerDim2_, int numOfKSpace_triggerDim1_);

    size_t getDimValue(const ISMRMRD::AcquisitionHeader& acqHeader, Gadgetron::ISMRMRDDIM& dim);
    void setDimValue(ISMRMRD::AcquisitionHeader& acqHeader, Gadgetron::ISMRMRDDIM& dim, size_t value);

    // buffer for per 2D image information
    GtPlusGadgetImageArray* messageImage_;

    // buffer for image kspace data
    // if the partial fourier is used, the kspace center is put at the center of buffer
    // this means zeros will be added accordingly
    BufferType kspaceBuffer_;
    BufferType refBuffer_;

    // dynamic buffer for other kspace data
    ReadOutBufferType noiseBuffer_;
    ReadOutBufferType phaseCorrBuffer_;
    ReadOutBufferType otherBuffer_;

    // dimension for image kspace
    std::vector<size_t> dimensions_;

    // encoding matrix size (the real sampled size)
    size_t matrix_size_encoding_[3];

    // minimal sampled line along E1 and E2
    long long min_sampled_E1_;
    long long min_sampled_E2_;

    // index of center line along E1 and E2
    long long center_line_E1_;
    long long center_line_E2_;

    // maximal sampled line along E1 and E2
    long long max_sampled_E1_;
    long long max_sampled_E2_;

    // encoding space size (the logic kspace size)
    size_t space_size_[3];

    // offset along E1 and E2 directions for incoming readouts
    long long space_matrix_offset_E1_;
    long long space_matrix_offset_E2_;

    // encoding filed of view [mm]
    float field_of_view_encoding_[3];

    // recon matrix size (the final image size)
    size_t matrix_size_recon_[3];

    // recon filed of view [mm]
    float field_of_view_recon_[3];

    // for the embedded mode
    size_t embedded_ref_lines_E1_;
    size_t embedded_ref_lines_E2_;

    size_t image_counter_;
    size_t image_series_;

    // mark the first kspace line
    bool first_kspace_scan_;

    // whether the next gadget has been triggered in close(...)
    bool triggered_in_close_;

    // whether the next gadget has been triggered in process(...)
    bool triggered_in_process_;

    // whether the next gadget has been triggered becasue the triggerDim1 changes meet the required number of kspace
    // only used for triggerDim1!=DIM_NONE && triggerDim2!=DIM_NONE
    bool triggered_in_process_by_numOfKSpace_triggerDim1_;

    // whether the next gadget has been triggered in process(...) for the last acquisition
    // if so, extra triggering in close(...) is not needed
    bool triggered_in_process_last_acq_;

    size_t meas_max_ro_;
    ISMRMRD::EncodingCounters meas_max_idx_;
    size_t meas_max_channel_;

    // maximal idx for reference data
    ISMRMRD::EncodingCounters meas_max_idx_ref_;

    // track the trigger dim1 and dim2
    size_t prev_dim1_;
    size_t curr_dim1_;

    size_t prev_dim2_;
    size_t curr_dim2_;

    // store the previous acquisition head
    ISMRMRD::AcquisitionHeader prev_acq_header_;

    // for trigger dim1, need to count its times
    size_t count_dim1_;

    // a general workorder to store the buffered data
    WorkOrderType workOrder_;

    // indicator for the arrival of last acq
    bool last_acq_arrived_;

    // time stamp resolution (default, 0.0025s)
    float timeStampResolution_;

    // exporter
    Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

    // util for gtplus
    Gadgetron::gtPlus::gtPlusISMRMRDReconUtil< std::complex<float> > gtPlus_util_;

    // in verbose mode, more info is printed out
    bool verboseMode_;

private:

    // index for the time stamp
    std::vector<size_t> ind_time_stamp_;
};

}
