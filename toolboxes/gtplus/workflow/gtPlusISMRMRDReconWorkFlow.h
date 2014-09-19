/** \file   gtPlusISMRMRDReconWorkFlow.h
    \brief  Define the base class for the GtPlus reconstruction workflow
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorker.h"

namespace Gadgetron { namespace gtPlus {

struct DimensionRecordCompare
{
    DimensionRecordCompare() {}
    ~DimensionRecordCompare() {}

    bool operator()(const std::pair<ISMRMRDDIM, size_t>& a, const std::pair<ISMRMRDDIM, size_t>& b) const
    {
        return (a.second > b.second);
    }
};

// [RO E1 CHA SLC E2 CON PHS REP SET SEG AVE]
#define GTPLUS_RECON_KSPACE_DIM_NUM 11
// [RO E1 CHA SLC E2 CON PHS REP SET AVE]
#define GTPLUS_RECON_IMAGE_DIM_NUM 10

template <typename T> 
class gtPlusISMRMRDReconWorkFlow
{
public:

    typedef std::pair<ISMRMRDDIM, size_t> DimensionRecordType;

    typedef typename realType<T>::Type real_value_type;

    gtPlusISMRMRDReconWorkFlow();
    gtPlusISMRMRDReconWorkFlow(gtPlusReconWorker<T>& worker, gtPlusReconWorkOrder<T>& workOrder);
    virtual ~gtPlusISMRMRDReconWorkFlow();

    void printInfo(std::ostream& os);

    virtual bool preProcessing() = 0;

    virtual bool recon() = 0;

    virtual bool postProcessing() = 0;

    // assemble the ISMRMRD dimension index
    // ind must have 10 elements
    bool ismrmrdDimIndex10D(std::vector<size_t>& ind, const ISMRMRDDIM& dim, size_t value);

    // find the permute order for ISMRMRD
    bool findISMRMRDPermuteOrder(const std::vector<ISMRMRDDIM>& dimsSrc, const std::vector<ISMRMRDDIM>& dimsDst, std::vector<size_t>& order);

    // print the dimension names
    std::string printISMRMRDDimensions(const std::vector<ISMRMRDDIM>& dims);

    // print the dimension size
    std::string printISMRMRDDimensionSize(const std::vector<size_t>& sizes);

    bool setDataArray(hoNDArray<T>& data);
    bool setDataArray(hoNDArray<T>& data, hoNDArray<real_value_type>& time_stamp, hoNDArray<real_value_type>& physio_time_stamp);
    bool setRefArray(hoNDArray<T>& ref);

    // -------- these member variables are made as public ------------- //

    // recon worker to do the computation
    gtPlusReconWorker<T>* worker_;

    // recon work order
    gtPlusReconWorkOrder<T>* workOrder_;

    // ----------------------------------
    // noise prewhitening
    // ----------------------------------
    // noise scan, 3D array [RO E1 CHA]
    hoNDArray<T>* noise_;

    // noise bandwidth (Hz/pixel)
    double noiseBW_;

    // noise equivalent bandwidth ratio for receiver
    double receriverBWRatio_;

    // ADC sampling time in second
    double ADCSamplingTimeinSecond_;

    // RO oversampling ratio
    double overSamplingRatioRO_;

    // ----------------------------------
    // final image sizes for RO/E1/E2
    // ----------------------------------
    size_t reconSizeRO_;
    size_t reconSizeE1_;
    size_t reconSizeE2_;

    float encodingFOV_RO_;
    float encodingFOV_E1_;
    float encodingFOV_E2_;

    float reconFOV_RO_;
    float reconFOV_E1_;
    float reconFOV_E2_;

    // ----------------------------------
    // dimension and starting indexes for this data_
    // in case this data_ is a portion of a larger dataset
    // ----------------------------------
    std::vector< DimensionRecordType > dataDimStartingIndexes_;

    // ----------------------------------
    // reconstruction results, complex images, 10D array [RO E1 CHA SLC E2 CON PHS REP SET AVE]
    // ----------------------------------
    hoNDArray<T> res_;
    // optional time stamps for the recon results, in the unit of seconds, 10D array [1 1 1 SLC E2 CON PHS REP SET AVE]
    // if not set, the stored image header will be used for time stamps
    hoNDArray<real_value_type> res_time_stamp_;
    hoNDArray<real_value_type> res_physio_time_stamp_;

    hoNDArray<T> res_second_;
    hoNDArray<real_value_type> res_time_stamp_second_;
    hoNDArray<real_value_type> res_physio_time_stamp_second_;

    // gfactor, not all reconstruction fills the gfactor
    // 10D array [RO E1 CHA SLC E2 CON PHS REP SET AVE]
    hoNDArray<T> gfactor_;

    // warp-around map, not all reconstruction fills the gfactor
    // 10D array [RO E1 2 SLC E2 CON PHS REP SET AVE]
    hoNDArray<T> wrap_around_map_;

    // ----------------------------------
    // debug and timing
    // ----------------------------------
    // clock for timing
    Gadgetron::GadgetronTimer gt_timer1_;
    Gadgetron::GadgetronTimer gt_timer2_;
    Gadgetron::GadgetronTimer gt_timer3_;

    bool performTiming_;

    // exporter
    Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

    // debug folder
    std::string debugFolder_;

    // ----------------------------------
    // input data array
    // ----------------------------------
    // image data, 11D [RO E1 CHA SLC E2 CON PHS REP SET SEG AVE]
    hoNDArray<T>* data_;
    // time stamp, 11D [1 E1 1 SLC E2 CON PHS REP SET SEG AVE]
    // these are set with data array
    hoNDArray<real_value_type>* time_stamp_;
    hoNDArray<real_value_type>* physio_time_stamp_;

    // reference calibration, 11D [RO E1 CHA SLC E2 CON PHS REP SET SEG AVE]
    hoNDArray<T>* ref_;

protected:

    // internal helper memory allocated for computation
    hoNDArray<T> dataCurr_;
    hoNDArray<T> refCurr_;
    hoNDArray<T> gfactorCurr_;
    hoNDArray<T> wrap_around_mapCurr_;

    // size of dimensions for image data
    DimensionRecordType RO_;
    DimensionRecordType E1_;
    DimensionRecordType CHA_;
    DimensionRecordType SLC_;
    DimensionRecordType E2_;
    DimensionRecordType CON_;
    DimensionRecordType PHS_;
    DimensionRecordType REP_;
    DimensionRecordType SET_;
    DimensionRecordType SEG_;
    DimensionRecordType AVE_;

    // size of dimensions for ref data
    DimensionRecordType RO_ref_;
    DimensionRecordType E1_ref_;
    DimensionRecordType CHA_ref_;
    DimensionRecordType SLC_ref_;
    DimensionRecordType E2_ref_;
    DimensionRecordType CON_ref_;
    DimensionRecordType PHS_ref_;
    DimensionRecordType REP_ref_;
    DimensionRecordType SET_ref_;
    DimensionRecordType SEG_ref_;
    DimensionRecordType AVE_ref_;

    // expected dimensions for results
    std::vector<ISMRMRDDIM> dimsRes_;

    // util
    gtPlusISMRMRDReconUtil<T> gtPlus_util_;
};

template <typename T> 
gtPlusISMRMRDReconWorkFlow<T>::gtPlusISMRMRDReconWorkFlow() 
: data_(NULL), time_stamp_(NULL), physio_time_stamp_(NULL), ref_(NULL), worker_(NULL), workOrder_(NULL), noise_(NULL), noiseBW_(1.0), receriverBWRatio_(1.0), overSamplingRatioRO_(1.0), ADCSamplingTimeinSecond_(1.0) , performTiming_(false)
{
    RO_.first = DIM_ReadOut;
    RO_.second = 1;

    E1_.first = DIM_Encoding1;
    E1_.second = 1;

    CHA_.first = DIM_Channel;
    CHA_.second = 1;

    SLC_.first = DIM_Slice;
    SLC_.second = 1;

    E2_.first = DIM_Encoding2;
    E2_.second = 1;

    CON_.first = DIM_Contrast;
    CON_.second = 1;

    PHS_.first = DIM_Phase;
    PHS_.second = 1;

    REP_.first = DIM_Repetition;
    REP_.second = 1;

    SET_.first = DIM_Set;
    SET_.second = 1;

    SEG_.first = DIM_Segment;
    SEG_.second = 1;

    AVE_.first = DIM_Average;
    AVE_.second = 1;

    RO_ref_.first = DIM_ReadOut;
    RO_ref_.second = 1;

    E1_ref_.first = DIM_Encoding1;
    E1_ref_.second = 1;

    CHA_ref_.first = DIM_Channel;
    CHA_ref_.second = 1;

    SLC_ref_.first = DIM_Slice;
    SLC_ref_.second = 1;

    E2_ref_.first = DIM_Encoding2;
    E2_ref_.second = 1;

    CON_ref_.first = DIM_Contrast;
    CON_ref_.second = 1;

    PHS_ref_.first = DIM_Phase;
    PHS_ref_.second = 1;

    REP_ref_.first = DIM_Repetition;
    REP_ref_.second = 1;

    SET_ref_.first = DIM_Set;
    SET_ref_.second = 1;

    SEG_ref_.first = DIM_Segment;
    SEG_ref_.second = 1;

    AVE_ref_.first = DIM_Average;
    AVE_ref_.second = 1;

    dimsRes_.resize(GTPLUS_RECON_IMAGE_DIM_NUM);
    dimsRes_[0] = DIM_ReadOut;
    dimsRes_[1] = DIM_Encoding1;
    dimsRes_[2] = DIM_Channel;
    dimsRes_[3] = DIM_Slice;
    dimsRes_[4] = DIM_Encoding2;
    dimsRes_[5] = DIM_Contrast;
    dimsRes_[6] = DIM_Phase;
    dimsRes_[7] = DIM_Repetition;
    dimsRes_[8] = DIM_Set;
    dimsRes_[9] = DIM_Average;

    dataDimStartingIndexes_.resize(GTPLUS_RECON_KSPACE_DIM_NUM);
    dataDimStartingIndexes_[0] = DimensionRecordType(DIM_ReadOut, 0);
    dataDimStartingIndexes_[1] = DimensionRecordType(DIM_Encoding1, 0);
    dataDimStartingIndexes_[2] = DimensionRecordType(DIM_Channel, 0);
    dataDimStartingIndexes_[3] = DimensionRecordType(DIM_Slice, 0);
    dataDimStartingIndexes_[4] = DimensionRecordType(DIM_Encoding2, 0);
    dataDimStartingIndexes_[5] = DimensionRecordType(DIM_Contrast, 0);
    dataDimStartingIndexes_[6] = DimensionRecordType(DIM_Phase, 0);
    dataDimStartingIndexes_[7] = DimensionRecordType(DIM_Repetition, 0);
    dataDimStartingIndexes_[8] = DimensionRecordType(DIM_Set, 0);
    dataDimStartingIndexes_[9] = DimensionRecordType(DIM_Segment, 0);
    dataDimStartingIndexes_[10] = DimensionRecordType(DIM_Average, 0);

    gt_timer1_.set_timing_in_destruction(false);
    gt_timer2_.set_timing_in_destruction(false);
    gt_timer3_.set_timing_in_destruction(false);
}

template <typename T> 
gtPlusISMRMRDReconWorkFlow<T>::gtPlusISMRMRDReconWorkFlow(gtPlusReconWorker<T>& worker, gtPlusReconWorkOrder<T>& workOrder) 
: data_(NULL), ref_(NULL), worker_(&worker), workOrder_(&workOrder), noise_(NULL),
  noiseBW_(1.0), receriverBWRatio_(1.0), overSamplingRatioRO_(1.0), ADCSamplingTimeinSecond_(1.0) , performTiming_(false)
{
    RO_.second = 1;
    E1_.second = 1;
    CHA_.second = 1;
    SLC_.second = 1;
    E2_.second = 1;
    CON_.second = 1;
    PHS_.second = 1;
    REP_.second = 1;
    SET_.second = 1;
    SEG_.second = 1;
    AVE_.second = 1;

    RO_ref_.second = 1;
    E1_ref_.second = 1;
    CHA_ref_.second = 1;
    SLC_ref_.second = 1;
    E2_ref_.second = 1;
    CON_ref_.second = 1;
    PHS_ref_.second = 1;
    REP_ref_.second = 1;
    SET_ref_.second = 1;
    SEG_ref_.second = 1;
    AVE_ref_.second = 1;

    dimsRes_.resize(GTPLUS_RECON_IMAGE_DIM_NUM);
    dimsRes_[0] = DIM_ReadOut;
    dimsRes_[1] = DIM_Encoding1;
    dimsRes_[2] = DIM_Channel;
    dimsRes_[3] = DIM_Slice;
    dimsRes_[4] = DIM_Encoding2;
    dimsRes_[5] = DIM_Contrast;
    dimsRes_[6] = DIM_Phase;
    dimsRes_[7] = DIM_Repetition;
    dimsRes_[8] = DIM_Set;
    dimsRes_[9] = DIM_Average;

    dataDimStartingIndexes_.resize(GTPLUS_RECON_KSPACE_DIM_NUM);
    dataDimStartingIndexes_[0] = DimensionRecordType(DIM_ReadOut, 0);
    dataDimStartingIndexes_[1] = DimensionRecordType(DIM_Encoding1, 0);
    dataDimStartingIndexes_[2] = DimensionRecordType(DIM_Channel, 0);
    dataDimStartingIndexes_[3] = DimensionRecordType(DIM_Slice, 0);
    dataDimStartingIndexes_[4] = DimensionRecordType(DIM_Encoding2, 0);
    dataDimStartingIndexes_[5] = DimensionRecordType(DIM_Contrast, 0);
    dataDimStartingIndexes_[6] = DimensionRecordType(DIM_Phase, 0);
    dataDimStartingIndexes_[7] = DimensionRecordType(DIM_Repetition, 0);
    dataDimStartingIndexes_[8] = DimensionRecordType(DIM_Set, 0);
    dataDimStartingIndexes_[9] = DimensionRecordType(DIM_Segment, 0);
    dataDimStartingIndexes_[10] = DimensionRecordType(DIM_Average, 0);

    gt_timer1_.set_timing_in_destruction(false);
    gt_timer2_.set_timing_in_destruction(false);
    gt_timer3_.set_timing_in_destruction(false);
}

template <typename T> 
gtPlusISMRMRDReconWorkFlow<T>::~gtPlusISMRMRDReconWorkFlow() 
{
}

template <typename T> 
void gtPlusISMRMRDReconWorkFlow<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD Recon workflow -------------" << endl;
    os << "Implementation of general reconstruction workflow for ISMRMRD convention" << endl;
    os << "the gtPlusISMRMRDReconWorkFlow defines and implements the reconstruction workflow for the ISMRMRD definition" << endl;
    os << "the reconstruction is split into three stages:" << endl;
    os << "1) PreProcessing" << endl;
    os << "2) Reconstruction" << endl;
    os << "3) PostProcessing" << endl;
    os << endl;
    os << "These three steps can have different operations for different sampling patterns or imaging applications" << endl;
    os << "----------------------------------------------------------" << endl;
}

template <typename T> 
inline bool gtPlusISMRMRDReconWorkFlow<T>::
ismrmrdDimIndex10D(std::vector<size_t>& ind, const ISMRMRDDIM& dim, size_t value)
{
    GADGET_CHECK_RETURN_FALSE(ind.size()>(size_t)(dim-DIM_ReadOut));
    ind[dim-DIM_ReadOut] = value;
    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlow<T>::
findISMRMRDPermuteOrder(const std::vector<ISMRMRDDIM>& dimsSrc, const std::vector<ISMRMRDDIM>& dimsDst, std::vector<size_t>& order)
{
    GADGET_CHECK_RETURN_FALSE(dimsSrc.size()==dimsDst.size());
    order.resize(dimsSrc.size());

    size_t NDim = dimsSrc.size();
    size_t src, dst;

    for ( dst=0; dst<NDim; dst++ )
    {
        for ( src=0; src<NDim; src++ )
        {
            if ( dimsSrc[src] == dimsDst[dst] )
                break;
        }

        order[dst] = src;
    }

    return true;
}

template <typename T> 
std::string gtPlusISMRMRDReconWorkFlow<T>::
printISMRMRDDimensions(const std::vector<ISMRMRDDIM>& dims)
{
    using namespace std;

    if ( dims.empty() ) return std::string("[ ]");

    size_t NDim = dims.size();

    size_t ii;

    std::ostringstream os;

    os << "[ ";
    for ( ii=0; ii<NDim; ii++ )
    {
        ISMRMRDDIM dim = dims[ii];
        switch (dim)
        {
            case DIM_ReadOut:
                os << "DIM_ReadOut ";
            break;

            case DIM_Encoding1:
                os << "Encoding1 ";
            break;

            case DIM_Channel:
                os << "Channel ";
            break;

            case DIM_Slice:
                os << "Slice ";
            break;

            case DIM_Encoding2:
                os << "Encoding2 ";
            break;

            case DIM_Contrast:
                os << "Contrast ";
            break;

            case DIM_Phase:
                os << "Phase ";
            break;

            case DIM_Repetition:
                os << "Repitition ";
            break;

            case DIM_Set:
                os << "Set ";
            break;

            case DIM_Segment:
                os << "Segment ";
            break;

            case DIM_Average:
                os << "Average ";
            break;

            default:
                os << " Other";
        }
    }
    os << "]" << endl;

    std::string dimStr(os.str());
    return dimStr;
}

template <typename T> 
std::string gtPlusISMRMRDReconWorkFlow<T>::
printISMRMRDDimensionSize(const std::vector<size_t>& sizes)
{
    using namespace std;

    if ( sizes.empty() ) return std::string("[ ]");

    size_t NDim = sizes.size();

    size_t ii;

    std::ostringstream os;

    os << "[ ";
    for ( ii=0; ii<NDim; ii++ )
    {
        os << sizes[ii] << " ";
    }
    os << "]" << endl;

    std::string sizeStr(os.str());
    return sizeStr;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlow<T>::setDataArray(hoNDArray<T>& data)
{
    try
    {
        data_ = &data;

        RO_.second = data.get_size(0);
        E1_.second = data.get_size(1);
        CHA_.second = data.get_size(2);
        SLC_.second = data.get_size(3);
        E2_.second = data.get_size(4);
        CON_.second = data.get_size(5);
        PHS_.second = data.get_size(6);
        REP_.second = data.get_size(7);
        SET_.second = data.get_size(8);
        SEG_.second = data.get_size(9);
        AVE_.second = data.get_size(10);
    }
    catch(...)
    {
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlow<T>::setDataArray(hoNDArray<T>& data, hoNDArray<real_value_type>& time_stamp, hoNDArray<real_value_type>& physio_time_stamp)
{
    try
    {
        data_ = &data;
        time_stamp_ = &time_stamp;
        physio_time_stamp_ = &physio_time_stamp;

        RO_.second = data.get_size(0);
        E1_.second = data.get_size(1);
        CHA_.second = data.get_size(2);
        SLC_.second = data.get_size(3);
        E2_.second = data.get_size(4);
        CON_.second = data.get_size(5);
        PHS_.second = data.get_size(6);
        REP_.second = data.get_size(7);
        SET_.second = data.get_size(8);
        SEG_.second = data.get_size(9);
        AVE_.second = data.get_size(10);
    }
    catch(...)
    {
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlow<T>::setRefArray(hoNDArray<T>& ref)
{
    try
    {
        ref_ = &ref;

        RO_ref_.second     = ref.get_size(0);
        E1_ref_.second     = ref.get_size(1);
        CHA_ref_.second    = ref.get_size(2);
        SLC_ref_.second    = ref.get_size(3);
        E2_ref_.second     = ref.get_size(4);
        CON_ref_.second    = ref.get_size(5);
        PHS_ref_.second    = ref.get_size(6);
        REP_ref_.second    = ref.get_size(7);
        SET_ref_.second    = ref.get_size(8);
        SEG_ref_.second    = ref.get_size(9);
        AVE_ref_.second    = ref.get_size(10);
    }
    catch(...)
    {
        return false;
    }

    return true;
}

}}
