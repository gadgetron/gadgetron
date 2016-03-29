/** \file   GtPlusReconGadget.h
    \brief  This is the base class gadget for both 2DT and 3DT reconstruction.
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "ismrmrd/meta.h"
#include "GadgetronTimer.h"

#include "hoNDArray_utils.h"

#include "GtPlusGadgetImageArray.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorkOrder.h"

#include "GadgetStreamController.h"

#include "GtPlusReconGadgetUtil.h"

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

    typedef std::pair<Gadgetron::ISMRMRDDIM, size_t> DimensionRecordType;

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

    // scaling factor for wrap around map
    double scalingFactor_wrap_around_map_;

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

    // whether the second set of recon results is required
    bool recon_res_second_required_;

    // whether to send out recon results
    bool send_out_recon_;
    bool send_out_recon_second_;

    // parameters for gt-plus recon
    Gadgetron::gtPlus::gtPlusReconWorkOrderPara workOrderPara_;

    // --------------------------------------------------
    // utility functions
    // --------------------------------------------------

    // compute image number using ICE way
    size_t computeSeriesImageNumber (ISMRMRD::ImageHeader& imheader, size_t nCHA=1, size_t cha=0, size_t nE2=1, size_t e2=0);

    // to handle partial fourier, add pre or post zeros
    // PrePostZeros: 0 no zeros; 1 pre zeros; 2 post zeros
    bool addPrePostZeros(int centreNo, int sampleNo, int& PrePostZeros);

    // find the dimension index
    bool findStartingDimIndex(const std::vector<DimensionRecordType>& dimStartingIndexes, Gadgetron::ISMRMRDDIM& dim, size_t ind);

    // compute SNR image and std map
    bool computeSNRImage(const hoNDArray<ValueType>& res, const hoNDArray<ValueType>& gfactor, unsigned int startInd, bool withAcceleration, hoNDArray<ValueType>& snrImage, hoNDArray<ValueType>& stdMap);

    // scale the recon images
    bool scalingImages(hoNDArray<ValueType>& res);

    // scale the magnitude images
    bool scalingMagnitude(hoNDArray<float>& mag);

    // recompute the image geometry parameters if the recon FOV is different from encoding FOV
    bool recomputeImageGeometry(GtPlusGadgetImageArray* images, GtPlusGadgetImageExt& imageHeader, size_t slc, size_t e2, size_t con, size_t phs, size_t rep, size_t set, size_t seg, size_t ave, size_t maxE2);

    // get the acquisition and PMU time stamps
    bool getTimeStamp(GtPlusGadgetImageArray* images, WorkOrderType& workOrder, hoNDArray<real_value_type>& timeStamp,  hoNDArray<real_value_type>& pmuTimeStamp);

    // send out the recon results
    virtual bool sendOutRecon(GtPlusGadgetImageArray* images, const hoNDArray<ValueType>& res, int seriesNum, const std::vector<DimensionRecordType>& dimStartingIndexes, const std::string& prefix, const std::string& dataRole);
    virtual bool sendOutRecon(GtPlusGadgetImageArray* images, const hoNDArray<ValueType>& res, const hoNDArray<real_value_type>& timeStamp, const hoNDArray<real_value_type>& physioTimeStamp, int seriesNum, const std::vector<DimensionRecordType>& dimStartingIndexes, const std::string& prefix, const std::string& dataRole);

    // special sending function for the interactive cases
    virtual bool sendOutRecon2D(GtPlusGadgetImageArray* images, const hoNDArray<float>& res, int seriesNum, int imageNum);
    virtual bool sendOutRecon2D(GtPlusGadgetImageArray* images, const hoNDArray<ValueType>& res, int seriesNum, int imageNum);

    // compute the kspace filter
    bool generateKSpaceFilter(WorkOrderType& workOrder);
    //void GDEBUG_CONDITION_STREAM(bool verboseMode_, const char* arg2);

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

    // close call
    int close(unsigned long flags);

public:

    /// ------------------------------------------------------------------------------------
    /// image scaling
    GADGET_PROPERTY(min_intensity_value, int, "Minimal intensity value for auto image scaling", 64);
    GADGET_PROPERTY(max_intensity_value, int, "Maxmimal intensity value for auto image scaling", 4095);
    GADGET_PROPERTY(scalingFactor, float, "Default scaling ratio", 4.0);
    GADGET_PROPERTY(scalingFactor_gfactor, float, "Default scaling ratio for gfactor", 100.0);
    GADGET_PROPERTY(scalingFactor_wrap_around_map, float, "Default scaling ratio for wrap around map", 100.0);
    GADGET_PROPERTY(scalingFactor_snr_image, float, "Default scaling ratio for snr map", 10.0);
    GADGET_PROPERTY(scalingFactor_std_map, float, "Default scaling ratio for std map", 1000.0);
    GADGET_PROPERTY(start_frame_for_std_map, int, "Start frame to compute std map", 5);
    GADGET_PROPERTY(use_constant_scalingFactor, bool, "Whether to use constrant scaling", true);

    /// ------------------------------------------------------------------------------------
    /// debug and image sending
    GADGET_PROPERTY(image_series, int, "Image series number", 0);
    GADGET_PROPERTY(verboseMode, bool, "Whether to print more information", false);
    GADGET_PROPERTY(debugFolder, std::string, "If set, the debug output will be written out", "");
    GADGET_PROPERTY(debugFolder2, std::string, "If set, the debug output will be written out", "");
    GADGET_PROPERTY(timeStampResolution, float, "Time tick resolution in second", 0.0025f);
    GADGET_PROPERTY(send_out_recon, bool, "Whether to send out recon images", true);
    GADGET_PROPERTY(send_out_recon_second, bool, "Whether to send out extra recon images", true);
    GADGET_PROPERTY(performTiming, bool, "Whether to perform timing on some computational steps", false);

    /// ------------------------------------------------------------------------------------
    /// kspace filter parameters
    GADGET_PROPERTY_LIMITS(filterRO, std::string, "Kspace filter for RO dimension", "Gaussian",
        GadgetPropertyLimitsEnumeration, "Gaussian", "Hanning", "TaperedHanning", "None");

    GADGET_PROPERTY(filterRO_sigma, double, "Filter sigma for gaussian for RO dimension", 1.0);
    GADGET_PROPERTY(filterRO_width, double, "Filter width for tapered hanning for RO dimension", 0.15);

    // ------------------------------------------------------------------------------------

    GADGET_PROPERTY_LIMITS(filterE1, std::string, "Kspace filter for E1 dimension", "Gaussian",
        GadgetPropertyLimitsEnumeration, "Gaussian", "Hanning", "TaperedHanning", "None");

    GADGET_PROPERTY(filterE1_sigma, double, "Filter sigma for gaussian for E1 dimension", 1.0);
    GADGET_PROPERTY(filterE1_width, double, "Filter width for tapered hanning for E1 dimension", 0.15);

    // ------------------------------------------------------------------------------------

    GADGET_PROPERTY_LIMITS(filterE2, std::string, "Kspace filter for E2 dimension", "Gaussian",
        GadgetPropertyLimitsEnumeration, "Gaussian", "Hanning", "TaperedHanning", "None");

    GADGET_PROPERTY(filterE2_sigma, double, "Filter sigma for gaussian for E2 dimension", 1.0);
    GADGET_PROPERTY(filterE2_width, double, "Filter width for tapered hanning for E2 dimension", 0.15);

    // ------------------------------------------------------------------------------------

    GADGET_PROPERTY_LIMITS(filterRefRO, std::string, "Kspace filter for ref data for RO dimension", "Hanning",
        GadgetPropertyLimitsEnumeration, "Hanning", "None");

    GADGET_PROPERTY(filterRefRO_sigma, double, "Filter sigma for gaussian for RO dimension", 1.5);
    GADGET_PROPERTY(filterRefRO_width, double, "Filter width for tapered hanning for RO dimension", 0.15);

    // ------------------------------------------------------------------------------------

    GADGET_PROPERTY_LIMITS(filterRefE1, std::string, "Kspace filter for ref data for E1 dimension", "Hanning",
        GadgetPropertyLimitsEnumeration, "Hanning", "None");

    GADGET_PROPERTY(filterRefE1_sigma, double, "Filter sigma for gaussian for E1 dimension", 1.5);
    GADGET_PROPERTY(filterRefE1_width, double, "Filter width for tapered hanning for E1 dimension", 0.15);

    // ------------------------------------------------------------------------------------

    GADGET_PROPERTY_LIMITS(filterRefE2, std::string, "Kspace filter for ref data for E2 dimension", "Hanning",
        GadgetPropertyLimitsEnumeration, "Hanning", "None");

    GADGET_PROPERTY(filterRefE2_sigma, double, "Filter sigma for gaussian for E2 dimension", 1.5);
    GADGET_PROPERTY(filterRefE2_width, double, "Filter width for tapered hanning for E2 dimension", 0.15);

    // ------------------------------------------------------------------------------------

    GADGET_PROPERTY_LIMITS(filterPartialFourierRO, std::string, "Kspace filter for partial fourier for RO dimension", "TaperedHanning",
        GadgetPropertyLimitsEnumeration, "TaperedHanning", "None");

    GADGET_PROPERTY(filterPartialFourierRO_sigma, double, "Partial fourier filter sigma for gaussian for RO dimension", 1.5);
    GADGET_PROPERTY(filterPartialFourierRO_width, double, "Partial fourier filter width for tapered hanning for RO dimension", 0.15);
    GADGET_PROPERTY(filterPartialFourierRO_densityComp, bool, "Whether to apply density compensation for RO dimension", false);

    // ------------------------------------------------------------------------------------

    GADGET_PROPERTY_LIMITS(filterPartialFourierE1, std::string, "Kspace filter for partial fourier for E1 dimension", "TaperedHanning",
        GadgetPropertyLimitsEnumeration, "TaperedHanning", "None");

    GADGET_PROPERTY(filterPartialFourierE1_sigma, double, "Partial fourier filter sigma for gaussian for E1 dimension", 1.5);
    GADGET_PROPERTY(filterPartialFourierE1_width, double, "Partial fourier filter width for tapered hanning for E1 dimension", 0.15);
    GADGET_PROPERTY(filterPartialFourierE1_densityComp, bool, "Whether to apply density compensation for E1 dimension", false);

    // ------------------------------------------------------------------------------------

    GADGET_PROPERTY_LIMITS(filterPartialFourierE2, std::string, "Kspace filter for partial fourier for E2 dimension", "TaperedHanning",
        GadgetPropertyLimitsEnumeration, "TaperedHanning", "None");

    GADGET_PROPERTY(filterPartialFourierE2_sigma, double, "Partial fourier filter sigma for gaussian for E2 dimension", 1.5);
    GADGET_PROPERTY(filterPartialFourierE2_width, double, "Partial fourier filter width for tapered hanning for E2 dimension", 0.15);
    GADGET_PROPERTY(filterPartialFourierE2_densityComp, bool, "Whether to apply density compensation for E2 dimension", false);

    /// ------------------------------------------------------------------------------------
    /// cloud computing
    GADGET_PROPERTY(CloudComputing, bool, "Whether to use cloud", false);
    GADGET_PROPERTY(cloudNodeFile, std::string, "Cloud node file", "my_Cloud.txt");
    GADGET_PROPERTY(CloudNodeXMLConfiguration, std::string, "Cloud node xml configuration file when using cloud bus", "GT_Cartesian_CloudNode.xml");

    /// ------------------------------------------------------------------------------------
    /// coil compression
    GADGET_PROPERTY(upstream_coil_compression, bool, "Whether to perform upstream coil compression", false);
    GADGET_PROPERTY(upstream_coil_compression_thres, double, "Threadhold for upstream coil compression", -1);
    GADGET_PROPERTY(upstream_coil_compression_num_modesKept, int, "Number of modes to keep for upstream coil compression", -1);

    GADGET_PROPERTY(downstream_coil_compression, bool, "Whether to perform downstream coil compression", true);
    GADGET_PROPERTY(coil_compression_thres, double, "Threadhold for downstream coil compression", 0.002);
    GADGET_PROPERTY(coil_compression_num_modesKept, int, "Number of modes to keep for downstream coil compression", -1);

    /// ------------------------------------------------------------------------------------
    /// coil map estimation
    GADGET_PROPERTY_LIMITS(coil_map_algorithm, std::string, "Coil map estimation method", "ISMRMRD_SOUHEIL",
        GadgetPropertyLimitsEnumeration, "ISMRMRD_SOUHEIL", "ISMRMRD_SOUHEIL_ITER");

    GADGET_PROPERTY(csm_kSize, int, "For ISMRMRD_SOUHEIL, kernel size", 7);
    GADGET_PROPERTY(csm_powermethod_num, int, "For ISMRMRD_SOUHEIL, number to apply power method", 3);
    GADGET_PROPERTY(csm_true_3D, bool, "Whether to use 3D kernel for 3D coil map estimation", true);
    GADGET_PROPERTY(csm_iter_num, int, "For ISMRMRD_SOUHEIL_ITER, number of iterations", 5);
    GADGET_PROPERTY(csm_iter_thres, double, "For ISMRMRD_SOUHEIL_ITER, iteration threshold", 1e-5);

    /// ------------------------------------------------------------------------------------
    /// recon parameters
    GADGET_PROPERTY_LIMITS(recon_algorithm, std::string, "Reconstruction algorithm", "ISMRMRD_GRAPPA",
        GadgetPropertyLimitsEnumeration, "ISMRMRD_GRAPPA", "ISMRMRD_SENSE", "ISMRMRD_SPIRIT", 
                            "ISMRMRD_SOFTSENSE", "ISMRMRD_L1SOFTSENSE", "ISMRMRD_2DTBINNING", "ISMRMRD_2DTBINNING_FLOW", 
                            "ISMRMRD_L1SPIRIT_SLEP", "ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP", "ISMRMRD_NONE");

    GADGET_PROPERTY(recon_auto_parameters, bool, "Whether to estimate recon algorithm parametes automatically", true);
    GADGET_PROPERTY(gfactor_needed, bool, "Whether to compute gfactor map", false);
    GADGET_PROPERTY(wrap_around_map_needed, bool, "Whether to compute wrap-around map", false);
    GADGET_PROPERTY(recon_kspace_needed, bool, "Whether to compute multi-channel full kspace", false);

    /// ------------------------------------------------------------------------------------
    /// grappa parameters
    GADGET_PROPERTY(grappa_kSize_RO, int, "Grappa kernel size RO", 5);
    GADGET_PROPERTY(grappa_kSize_E1, int, "Grappa kernel size E1", 4);
    GADGET_PROPERTY(grappa_kSize_E2, int, "Grappa kernel size E2", 4);
    GADGET_PROPERTY(grappa_reg_lamda, double, "Grappa regularization threshold", 0.0005);
    GADGET_PROPERTY(grappa_calib_over_determine_ratio, double, "Grappa calibration overdermination ratio", 0);

    /// ------------------------------------------------------------------------------------
    /// spirit parameters
    GADGET_PROPERTY(spirit_kSize_RO, int, "Spirit kernel size RO", 7);
    GADGET_PROPERTY(spirit_kSize_E1, int, "Spirit kernel size E1", 7);
    GADGET_PROPERTY(spirit_kSize_E2, int, "Spirit kernel size E2", 5);
    GADGET_PROPERTY(spirit_reg_lamda, double, "Spirit regularization threshold", 0.005);
    GADGET_PROPERTY(spirit_calib_over_determine_ratio, double, "Spirit calibration overdermination ratio", 0);
    GADGET_PROPERTY(spirit_solve_symmetric, bool, "Whether to solve symmetric spirit kernel size E2", false);
    GADGET_PROPERTY(spirit_iter_max, int, "Spirit number of iterations", 150);
    GADGET_PROPERTY(spirit_iter_thres, double, "Spirit threshold for iteration", 0.0015);
    GADGET_PROPERTY(spirit_print_iter, bool, "Spirit print iteration info", false);

    GADGET_PROPERTY(spirit_perform_linear, bool, "Whether to perform linear spirit", true);
    GADGET_PROPERTY(spirit_perform_nonlinear, bool, "Whether to perform non-linear spirit", true);
    GADGET_PROPERTY(spirit_parallel_imaging_lamda, double, "Parallel imaging term lamda", 1.0);
    GADGET_PROPERTY(spirit_image_reg_lamda, double, "Imaging regularization term lamda", 0.00005);
    GADGET_PROPERTY(spirit_data_fidelity_lamda, double, "Data fedility term lamda", 1.0);

    GADGET_PROPERTY(spirit_ncg_iter_max, int, "Nonlinear cg solver numbre of iterations", 10);
    GADGET_PROPERTY(spirit_ncg_iter_thres, double, "Nonlinear cg solver iteration threshold", 0.0001);
    GADGET_PROPERTY(spirit_ncg_print_iter, bool, "Spirit print nonlinear cg iteration info", true);

    GADGET_PROPERTY(spirit_use_coil_sen_map, bool, "Spirit using coil map in image regularization", false);
    GADGET_PROPERTY(spirit_use_moco_enhancement, bool, "Spirit using motion correction enhancement in image regularization", false);
    GADGET_PROPERTY(spirit_recon_moco_images, bool, "Spirit recon motion corrected image series", false);
    GADGET_PROPERTY(spirit_RO_enhancement_ratio, double, "RO enhancement ratio for image regularization", 1.0);
    GADGET_PROPERTY(spirit_E1_enhancement_ratio, double, "E1 enhancement ratio for image regularization", 1.0);
    GADGET_PROPERTY(spirit_E2_enhancement_ratio, double, "E2 enhancement ratio for image regularization", 1.0);
    GADGET_PROPERTY(spirit_temporal_enhancement_ratio, double, "Temporal enhancement ratio for image regularization", 5.0);

    GADGET_PROPERTY(spirit_2D_scale_per_chunk, bool, "Spirit scaling per chunk for 2DT", false);
    GADGET_PROPERTY(spirit_3D_scale_per_chunk, bool, "Spirit scalng per chunk for 3D", false);

    /// ------------------------------------------------------------------------------------
    /// retro gating parameters
    GADGET_PROPERTY(retro_gated_interp_method, std::string, "Interpolation method used for retro-gating", "ISMRMRD_INTERP_RETRO_GATING_LINEAR");

    /// ------------------------------------------------------------------------------------
    /// recon job parameters
    GADGET_PROPERTY(job_split_by_S, bool, "Every S leads to a recon job", false);
    GADGET_PROPERTY(job_num_of_N, int, "Recon job size along N", 32);
    GADGET_PROPERTY(job_max_Megabytes, int, "Maximal recon job size in MegaBytes", 2048);
    GADGET_PROPERTY(job_overlap, int, "Recon job overlap size", 2);
    GADGET_PROPERTY(job_perform_on_control_node, bool, "Whether to perform recon job on control node", false);

    /// ------------------------------------------------------------------------------------
    /// partial fourier handling parameters
    GADGET_PROPERTY_LIMITS(partialFourier_algo, std::string, "Partial fourier handling method", "ISMRMRD_PF_ZEROFILLING_FILTER",
        GadgetPropertyLimitsEnumeration, "ISMRMRD_PF_HOMODYNE", "ISMRMRD_PF_FENGHUANG", "ISMRMRD_PF_POCS", "ISMRMRD_PF_ZEROFILLING_FILTER", "ISMRMRD_PF_ZEROFILLING", "ISMRMRD_PF_NONE");

    GADGET_PROPERTY(partialFourier_homodyne_iters, int, "Number of iterations for homodyne PF handling", 6);
    GADGET_PROPERTY(partialFourier_homodyne_thres, double, "Threshold for homodyne PF handling", 0.01);
    GADGET_PROPERTY(partialFourier_homodyne_densityComp, bool, "Whether to perform density compensation for homodyne PF handling", false);

    GADGET_PROPERTY(partialFourier_POCS_iters, int, "Number of iterations for POCS PF handling", 6);
    GADGET_PROPERTY(partialFourier_POCS_thres, double, "Threshold for POSC PF handling", 0.01);
    GADGET_PROPERTY(partialFourier_POCS_transitBand, int, "Transition band width for POCS PF handling", 24);
    GADGET_PROPERTY(partialFourier_POCS_transitBand_E2, int, "Transition band width for POCS PF handling for E2 dimension", 16);

    GADGET_PROPERTY(partialFourier_FengHuang_kSize_RO, int, "RO kernel size for FengHuang PF handling", 5);
    GADGET_PROPERTY(partialFourier_FengHuang_kSize_E1, int, "E1 kernel size for FengHuang PF handling", 5);
    GADGET_PROPERTY(partialFourier_FengHuang_kSize_E2, int, "E2 kernel size for FengHuang PF handling", 5);
    GADGET_PROPERTY(partialFourier_FengHuang_thresReg, double, "Threshold for FengHuang PF handling", 0.01);
    GADGET_PROPERTY(partialFourier_FengHuang_sameKernel_allN, bool, "Whether all N have the same kernel for FengHuang PF handling", false);
    GADGET_PROPERTY(partialFourier_FengHuang_transitBand, int, "Transition band width for FengHuang PF handling", 24);
    GADGET_PROPERTY(partialFourier_FengHuang_transitBand_E2, int, "Transition band width for FengHuang PF handling for E2 dimension", 16);

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

    // encoding limits for E1 and E2
    // used to set up kspace filter
    size_t min_E1_;
    size_t max_E1_;
    size_t center_E1_;

    size_t min_E2_;
    size_t max_E2_;
    size_t center_E2_;

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
    Gadgetron::ISMRMRDCALIBMODE CalibMode_;
    Gadgetron::ISMRMRDDIM InterleaveDim_;

    // acquired max indexes
    size_t kSpaceMaxAcqE1No_;
    size_t kSpaceMaxAcqE2No_;

    // number of times the process function is called
    unsigned int processed_called_times_;

    // kspace filter for RO/E1/E2
    // for the partial fourier, zero-padding resize or asymmetric echo
    // if the kspace filter is not selected, the default filter will be used anyway

    // kspace filter
    Gadgetron::ISMRMRDKSPACEFILTER filterRO_type_;
    double filterRO_sigma_;
    double filterRO_width_;

    Gadgetron::ISMRMRDKSPACEFILTER filterE1_type_;
    double filterE1_sigma_;
    double filterE1_width_;

    Gadgetron::ISMRMRDKSPACEFILTER filterE2_type_;
    double filterE2_sigma_;
    double filterE2_width_;

    // ref data filter
    Gadgetron::ISMRMRDKSPACEFILTER filterRO_ref_type_;
    double filterRO_ref_sigma_;
    double filterRO_ref_width_;

    Gadgetron::ISMRMRDKSPACEFILTER filterE1_ref_type_;
    double filterE1_ref_sigma_;
    double filterE1_ref_width_;

    Gadgetron::ISMRMRDKSPACEFILTER filterE2_ref_type_;
    double filterE2_ref_sigma_;
    double filterE2_ref_width_;

    // partial fourier filter
    Gadgetron::ISMRMRDKSPACEFILTER filterRO_pf_type_;
    double filterRO_pf_sigma_;
    double filterRO_pf_width_;
    bool filterRO_pf_densityComp_;

    Gadgetron::ISMRMRDKSPACEFILTER filterE1_pf_type_;
    double filterE1_pf_sigma_;
    double filterE1_pf_width_;
    bool filterE1_pf_densityComp_;

    Gadgetron::ISMRMRDKSPACEFILTER filterE2_pf_type_;
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

    // define the maximal number of threads used
    // number_of_used_threads = thread_number_ratio_ * max_available_threads_number
    // 0 means all threads are used
    float thread_number_ratio_;

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
};

}
