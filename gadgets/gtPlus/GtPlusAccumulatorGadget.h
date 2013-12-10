#pragma once

#include <complex>
#include "GtPlusExport.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd.h"
#include "GadgetIsmrmrdReadWrite.h"

// the buffered kspace is defined by the ISMRM 10 dimensions
// readout
// kspace_encode_step_1
// kspace_encode_step_2
// average
// slice
// contrast
// phase
// repetition
// set
// segment
// in the order of [RO E1 CHA AVE SLC E2 CON PHS REP SET SEG]

namespace Gadgetron
{

struct  EXPORTGTPLUS GadgetMessageImageExt : public ISMRMRD::ImageHeader
{
    // fields added to store the time_stamp and pmu_time_stamp for every incoming read-out line
    // if one line is not acquried, the corresponding time is -1
    std::vector<int>     time_stamps;
    std::vector<int>     pmu_time_stamps;

    GadgetMessageImageExt();
    ~GadgetMessageImageExt();

    void copy(GadgetMessageImageExt& aMessageImage);
    void set_matrix_size(unsigned int index, ACE_UINT16 size);
    void dump();
}; 

// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
//   0  1  2   3  4  5    6     7   8   9
// store a scan with 10 dimensions
struct  EXPORTGTPLUS GadgetMessageImageArray
{
    // size of the image array
    ACE_UINT16 matrix_size[10];

    // message information for every 2D image [Slice E2 Contrast Phase Rep Set Seg]
    GadgetMessageImageExt* imageArray_;

    GadgetMessageImageArray();
    GadgetMessageImageArray(int aSize[10]);
    ~GadgetMessageImageArray();

    void resize(int aSize[10]);
    void copy(GadgetMessageImageArray& imageArray);
    int get_offset(int slc, int par, int eco, int phs, int rep, int set, int seg);
    void extractMessageImageArrayForSLC(int slc, GadgetMessageImageArray& imageArray);
    void extractMessageImageArrayForREP(int rep, GadgetMessageImageArray& imageArray);

    void dump();
};

struct EXPORTGTPLUS KSpaceBuffer
{
    typedef hoNDArray< std::complex<float> > BufferType;
    typedef hoNDArray< unsigned short > ReflectBufferType;

    // reflect buffer shows whether a readouline is reflected or not

    // kspace data
    BufferType buffer_;
    ReflectBufferType reflect_;

    // reference ACS data
    BufferType ref_;
    ReflectBufferType refReflect_;

    // noise data
    BufferType noise_;

    // phase correction data
    BufferType phaseCorr_;
    ReflectBufferType phaseCorrReflect_;

    // other data, e.g. AIF data
    BufferType other_;

    // properties of kspace
    // kspace center readout number
    unsigned int kSpaceCentreRO_;
    // kspace center number for the first encoding dimension
    unsigned int kSpaceCentreEncode1_;
    // kspace center number for the second encoding dimension
    unsigned int kSpaceCentreEncode2_;

    // kspace max acquired readout number
    unsigned int kSpaceMaxRO_;
    // kspace max acquired number for the first encoding dimension
    unsigned int kSpaceMaxEncode1_;
    // kspace max acquired number for the second encoding dimension
    unsigned int kSpaceMaxEncode2_;

    // acceleration rate along the E1 and E2 dimensions
    unsigned int AccelFactE1_;
    unsigned int AccelFactE2_;

    // mode of calibration
    ISMRMRD::calibrationModeType::value CalibMode_;
    ISMRMRD::interleavingDimensionType::value InterleaveDim_;

    KSpaceBuffer();
    ~KSpaceBuffer();
};

// -----------------------------------------------------------------------------------------------------------

struct ReadOutBuffer
{
    ISMRMRD::AcquisitionHeader acqHead_;
    hoNDArray< std::complex<float> > data_;
    bool isReflect_;
};

class EXPORTGTPLUS GtPlusAccumulatorGadget : public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
{
public:
    GADGET_DECLARE(GtPlusAccumulatorGadget);

    typedef std::complex<float> ValueType;

    typedef Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< ValueType > > BaseClass;

    typedef std::vector< ReadOutBuffer > ReadOutBufferType;
    typedef hoNDArray< std::complex<float> > BufferType;
    typedef hoNDArray< unsigned short > ReflectBufferType;

    GtPlusAccumulatorGadget();
    ~GtPlusAccumulatorGadget();

    virtual int close(unsigned long flags);

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
    virtual bool checkStatus(uint64_t flag, int samples, bool& bIsKSpace, bool& bIsRef, bool& bIsNoise, bool& bIsPhaseCorr, bool& bIsReflect, bool& bIsOther);

    // store the image data
    virtual bool storeImageData(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2, bool isReflect);

    // fill the dynamically buffered data
    virtual bool fillBuffer(ReadOutBufferType& readOutBuffer, BufferType& buf, ReflectBufferType& reflectBuf);

    // fill the per 2D image info
    virtual bool fillImageInfo(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1, GadgetMessageImageArray* messageImage, const ISMRMRD::EncodingCounters& idx);

    // buffer for per 2D image information
    GadgetMessageImageArray* messageImage_;

    // buffer for image kspace data
    // if the partial fourier is used, the kspace center is put at the center of buffer
    // this means zeros will be padded accordingly
    KSpaceBuffer* kspaceBuffer_;

    // dynamic buffer for other kspace data
    ReadOutBufferType refBuffer_;
    ReadOutBufferType noiseBuffer_;
    ReadOutBufferType phaseCorrBuffer_;
    ReadOutBufferType otherBuffer_;

    // dimension for image kspace
    std::vector<unsigned int> dimensions_;

    // filed of view [mm]
    float field_of_view_[3];

    int image_counter_;
    int image_series_;

    // whether the next gadget has been triggered
    bool triggered_;

    int meas_max_ro_;
    ISMRMRD::EncodingCounters meas_max_idx_;
    int meas_max_channel_;
};

}
