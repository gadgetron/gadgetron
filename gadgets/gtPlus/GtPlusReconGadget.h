/** \file   GtPlusReconGadget.h
    \brief  This is the base class gadget for both 2DT and 3DT reconstruction.
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "hoNDMetaAttributes.h"
#include "ismrmrd.h"
#include "ismrmrd_xml.h"
#include "GadgetronTimer.h"

#include "hoNDArray_utils.h"

#include "GtPlusGadgetImageArray.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorkOrder.h"

#include "GadgetStreamController.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

#define SNR_NOISEFLOOR_SCALEFACTOR 8

namespace Gadgetron
{

class EXPORTGTPLUSGADGET GtPlusReconGadget : public Gadgetron::Gadget2< GtPlusGadgetImageArray, Gadgetron::gtPlus::gtPlusReconWorkOrder<std::complex<float> > >
{
public:
    GADGET_DECLARE(GtPlusReconGadget);

    typedef float real_value_type;
    typedef std::complex<real_value_type> ValueType;

    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder<ValueType> WorkOrderType;

    typedef Gadget2< GtPlusGadgetImageArray, WorkOrderType > BaseClass;

    typedef std::pair<Gadgetron::gtPlus::ISMRMRDDIM, size_t> DimensionRecordType;

    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder<ValueType>::CloudNodeType CloudNodeType;
    typedef std::vector<CloudNodeType> CloudType;

    GtPlusReconGadget();
    ~GtPlusReconGadget();

    // image series number
    int image_series_;

    // the min/max dynamic range of magnitude images
    size_t min_intensity_value_;
    size_t max_intensity_value_;

    // maximal intensity value when converted to unsigned short
    size_t max_intensity_value_US_;

    // scaling factor for recon results
    double scalingFactor_;

    // scaling factor for gfactor images
    double scalingFactor_gfactor_;

    // scaling factor for snr images
    double scalingFactor_snr_image_;

    // scaling factor for std map
    double scalingFactor_std_map_;

    // start frame to compute std map, to avoid transitional signal
    unsigned int start_frame_for_std_map_;

    // whether to use the fixed intensity scaling factor
    bool use_constant_scalingFactor_;

    // time stamp resolution (default, 2.5ms)
    float timeStampResolution_;

    // pixel spacing when exporting the images
    double aSpacing_[6];

    // field of view in mm
    double FOV_RO_;
    double FOV_E1_;
    double FOV_E2_;

    // debug folder
    std::string debugFolder_;
    std::string debugFolder_fullPath_;

    // debug folder 2
    std::string debugFolder2_;
    std::string debugFolder2_fullPath_;

    // whether to perform timing
    bool performTiming_;

    // whether to recon kspace
    bool recon_kspace_needed_;

    // parameters for gt-plus recon
    Gadgetron::gtPlus::gtPlusReconWorkOrderPara workOrderPara_;

    // --------------------------------------------------
    // utility functions
    // --------------------------------------------------

    // generate the debug folder path
    // debugFolderPath = ${GADGETRON_HOME}/debugFolder
    virtual bool generateDebugFolderPath(const std::string& debugFolder, std::string& debugFolderPath);

    // get the current moment
    void getCurrentMoment(std::string& procTime);

    // compute image number using ICE way
    size_t computeSeriesImageNumber (ISMRMRD::ImageHeader& imheader, size_t nCHA=1, size_t cha=0, size_t nE2=1, size_t e2=0);

    // to handle partial fourier, add pre or post zeros
    // PrePostZeros: 0 no zeros; 1 pre zeros; 2 post zeros
    bool addPrePostZeros(int centreNo, int sampleNo, int& PrePostZeros);

    // find the dimension index
    bool findStartingDimIndex(const std::vector<DimensionRecordType>& dimStartingIndexes, Gadgetron::gtPlus::ISMRMRDDIM& dim, size_t ind);

    // compute SNR image and std map
    bool computeSNRImage(const hoNDArray<ValueType>& res, const hoNDArray<ValueType>& gfactor, unsigned int startInd, bool withAcceleration, hoNDArray<ValueType>& snrImage, hoNDArray<ValueType>& stdMap);

    // scale the recon images
    bool scalingImages(hoNDArray<ValueType>& res);

    // scale the magnitude images
    bool scalingMagnitude(hoNDArray<float>& mag);

    // recompute the image geometry parameters if the recon FOV is different from encoding FOV
    bool recomputeImageGeometry(GtPlusGadgetImageArray* images, GtPlusGadgetImageExt& imageHeader, size_t slc, size_t e2, size_t con, size_t phs, size_t rep, size_t set, size_t seg, size_t maxE2);

    // send out the recon results
    virtual bool sendOutRecon(GtPlusGadgetImageArray* images, const hoNDArray<ValueType>& res, int seriesNum, const std::vector<DimensionRecordType>& dimStartingIndexes, const std::string& prefix, const std::string& dataRole);

    // special sending function for the interactive cases
    virtual bool sendOutRecon2D(GtPlusGadgetImageArray* images, const hoNDArray<float>& res, int seriesNum, int imageNum);
    virtual bool sendOutRecon2D(GtPlusGadgetImageArray* images, const hoNDArray<ValueType>& res, int seriesNum, int imageNum);

    // compute the kspace filter
    bool generateKSpaceFilter(WorkOrderType& workOrder);

protected:

    // --------------------------------------------------
    // functional functions
    // --------------------------------------------------

    // default interface function
    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(Gadgetron::GadgetContainerMessage< GtPlusGadgetImageArray >* m1, Gadgetron::GadgetContainerMessage< WorkOrderType > * m2);

    // read in parameters
    virtual bool readParameters();

    // parse the cloud file if any
    virtual bool parseGTCloudNodeFile(const std::string& filename, CloudType& gtCloud);

public:

    // --------------------------------------------------
    // variables used for data buffer and processing
    // --------------------------------------------------

    // dimension of incoming array
    std::vector<size_t> dimensions_;

    // number of acquisition channels
    size_t num_acq_channels_;

    // encoding matrix size (the real sampled size)
    size_t matrix_size_encoding_[3];

    // encoding filed of view [mm]
    float field_of_view_encoding_[3];

    // recon matrix size (the final image size)
    size_t matrix_size_recon_[3];

    // recon filed of view [mm]
    float field_of_view_recon_[3];

    // number of E1/E2 after zero-filling resize
    size_t reconE1_;
    size_t reconE2_;

    // acceleration factor
    double acceFactorE1_;
    double acceFactorE2_;

    // calibration mode
    Gadgetron::gtPlus::ISMRMRDCALIBMODE CalibMode_;
    Gadgetron::gtPlus::ISMRMRDDIM InterleaveDim_;

    // acquired max indexes
    size_t kSpaceMaxAcqE1No_;
    size_t kSpaceMaxAcqE2No_;

    // number of times the process function is called
    unsigned int processed_called_times_;

    // kspace filter for RO/E1/E2
    // for the partial fourier, zero-padding resize or asymmetric echo
    // if the kspace filter is not selected, the default filter will be used anyway

    // kspace filter
    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterRO_type_;
    double filterRO_sigma_;
    double filterRO_width_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE1_type_;
    double filterE1_sigma_;
    double filterE1_width_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE2_type_;
    double filterE2_sigma_;
    double filterE2_width_;

    // ref data filter
    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterRO_ref_type_;
    double filterRO_ref_sigma_;
    double filterRO_ref_width_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE1_ref_type_;
    double filterE1_ref_sigma_;
    double filterE1_ref_width_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE2_ref_type_;
    double filterE2_ref_sigma_;
    double filterE2_ref_width_;

    // partial fourier filter
    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterRO_pf_type_;
    double filterRO_pf_sigma_;
    double filterRO_pf_width_;
    bool filterRO_pf_densityComp_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE1_pf_type_;
    double filterE1_pf_sigma_;
    double filterE1_pf_width_;
    bool filterE1_pf_densityComp_;

    Gadgetron::gtPlus::ISMRMRDKSPACEFILTER filterE2_pf_type_;
    double filterE2_pf_sigma_;
    double filterE2_pf_width_;
    bool filterE2_pf_densityComp_;

    /// cloud related definition
    bool CloudComputing_;
    unsigned int CloudSize_;

    CloudType gt_cloud_;

    // cloud node file
    std::string cloud_node_file_;

    // encoding space size
    ISMRMRD::EncodingCounters meas_max_idx_;

    Gadgetron::gtPlus::gtPlusISMRMRDReconUtil<ValueType> gtPlus_util_;
    Gadgetron::gtPlus::gtPlusISMRMRDReconUtilComplex<ValueType> gtPlus_util_complex_;

    // clock for timing
    Gadgetron::GadgetronTimer gt_timer1_;
    Gadgetron::GadgetronTimer gt_timer2_;
    Gadgetron::GadgetronTimer gt_timer3_;

    // exporter
    Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

    // in verbose mode, more info is printed out
    bool verboseMode_;

    // memory manager
    boost::shared_ptr<Gadgetron::gtPlus::gtPlusMemoryManager> mem_manager_;
};

}
