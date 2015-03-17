/** \file   gtPlusISMRMRDReconWorkOrder.h
    \brief  Define the GtPlus reconstruction workorder and parameters
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"
#include "gtPlusISMRMRDReconUtil.h"

// MACROS FOR PRINTING
#define GADGET_PARA_PRINT(content) { GDEBUG_STREAM(#content << " is " << content); }

namespace Gadgetron { namespace gtPlus {

#define MAX_MOCO_LEVEL 16

struct gtPlusReconWorkOrderPara
{
    ISMRMRDCALIBMODE CalibMode_;
    ISMRMRDDIM InterleaveDim_;

    // acceleration factor along E1 and E2
    double acceFactorE1_;
    double acceFactorE2_;

    // kspace center for RO/E1/E2
    size_t kSpaceCenterRO_;
    size_t kSpaceCenterEncode1_;
    size_t kSpaceCenterEncode2_;

    // kspace max acquired number for RO/E1/E2
    size_t kSpaceMaxRO_;
    size_t kSpaceMaxEncode1_;
    size_t kSpaceMaxEncode2_;

    // for asymmetric echo
    // the sampled range for RO
    // if <0, all RO ranges are used
    int start_RO_;
    int end_RO_;

    // sampled range for E1
    int start_E1_;
    int end_E1_;

    // sampled range for E2
    int start_E2_;
    int end_E2_;

    // work order has to have some interaction with work flow
    // if true work flow will buffer kernel computed from this work order
    bool workFlow_BufferKernel_;
    // if true, work flow will use its buffered kernel for this work order
    bool workFlow_use_BufferedKernel_;

    // number of channels for the reconstruction results
    // most cases, it is 1
    size_t num_channels_res_;

    // -----------------------------------------------------------------------
    // parameters
    // -----------------------------------------------------------------------

    // -------------------------------
    // coil compression
    // -------------------------------
    bool upstream_coil_compression_;
    double upstream_coil_compression_thres_;
    int upstream_coil_compression_num_modesKept_;

    bool downstream_coil_compression_;
    double coil_compression_thres_;
    int coil_compression_num_modesKept_;

    // -------------------------------
    // coil sensitivity estimation
    // -------------------------------
    Gadgetron::ISMRMRDCOILMAPALGO coil_map_algorithm_;

    // for ISMRMRD_SOUHEIL
    size_t csm_kSize_;
    size_t csm_powermethod_num_;
    // for 3D acquisition, whether to use the true 3D data correlation matrix
    bool csm_true_3D_;

    // for ISMRMRD_SOUHEIL_ITER
    size_t csm_iter_num_;
    double csm_iter_thres_;

    // whether to use gpu for csm estimation
    bool csm_use_gpu_;

    // -------------------------------
    // parameters for variant reconstruction algorithms
    // -------------------------------
    Gadgetron::ISMRMRDALGO recon_algorithm_;
    bool recon_auto_parameters_;

    bool gfactor_needed_;

    bool wrap_around_map_needed_;

    /// --------------
    // grappa
    /// --------------

    size_t grappa_kSize_RO_;
    size_t grappa_kSize_E1_;
    size_t grappa_kSize_E2_;
    double grappa_reg_lamda_;
    double grappa_calib_over_determine_ratio_;
    bool grappa_use_gpu_;

    /// --------------
    // SPIRiT
    /// --------------
    size_t spirit_kSize_RO_;
    size_t spirit_kSize_E1_;
    size_t spirit_kSize_E2_;

    size_t spirit_oSize_RO_;
    size_t spirit_oSize_E1_;
    size_t spirit_oSize_E2_;

    double spirit_reg_lamda_;
    double spirit_calib_over_determine_ratio_;

    bool spirit_solve_symmetric_;

    size_t spirit_iter_max_;
    double spirit_iter_thres_;
    bool spirit_print_iter_;

    bool spirit_use_gpu_;

    /// --------------
    // L1 SPIRiT
    /// --------------
    bool spirit_perform_linear_;
    bool spirit_perform_grappa_linear_;
    bool spirit_perform_nonlinear_;

    double spirit_parallel_imaging_lamda_;
    double spirit_image_reg_lamda_;
    double spirit_data_fidelity_lamda_;

    size_t spirit_ncg_iter_max_;
    double spirit_ncg_iter_thres_;
    bool spirit_ncg_print_iter_;
    double spirit_ncg_scale_factor_;

    size_t spirit_slep_iter_max_;
    double spirit_slep_iter_thres_;
    bool spirit_slep_print_iter_;
    bool spirit_slep_keep_third_dimension_coeff_;
    bool spirit_slep_keep_approx_coeff_;
    double spirit_slep_scale_factor_;

    bool spirit_use_coil_sen_map_;
    bool spirit_use_moco_enhancement_;
    bool spirit_recon_moco_images_;

    bool spirit_2D_scale_per_chunk_;
    bool spirit_3D_scale_per_chunk_;

    double spirit_RO_enhancement_ratio_;
    double spirit_E1_enhancement_ratio_;
    double spirit_E2_enhancement_ratio_;
    double spirit_temporal_enhancement_ratio_;

    /// --------------
    /// parameters for retro-gating
    /// --------------
    // number of retro-gated phases
    // if 0, retro-gating is not prescribed
    size_t retro_gated_images_;

    // how many readout lines in each segment for retro-gating
    size_t retro_gated_segment_size_;

    // which method used for retro-gating
    ISMRMRDINTERPRETROGATING retro_gated_interp_method_;

    /// --------------
    /// parameters for binning
    /// --------------

    // number of target cardiac phases
    size_t kspace_binning_number_of_cardiac_phases_;

    // minimal allowed cardiac phase width used for binning, in ms
    // if the binned temporal window is smaller than this threshold,
    // the binned window will be increased
    // if <=0, then this value will not take effect
    double kspace_binning_minimal_cardiac_phase_width_;

    // whether to perform binning recon with multiple channel complex data
    bool kspace_binning_multiple_channel_recon_;

    // whether to perform non-linear recon
    bool kspace_binning_iterative_non_linear_recon_;

    // non-linear recon using slep optimizer
    bool kspace_binning_iterative_non_linear_recon_slep_;

    // whether to use coil map when warpping multiple channel images
    bool kspace_binning_multiple_channel_recon_with_coil_map_;

    // whether to compute navigator signal
    bool kspace_binning_compute_navigator_signal_;

    // for navigator detection
    size_t kspace_binning_navigator_moco_level_;
    size_t kspace_binning_navigator_moco_iter_[MAX_MOCO_LEVEL];
    double kspace_binning_navigator_hilbert_strength_;
    double kspace_binning_navigator_dissimilarity_sigma_;
    bool  kspace_binning_navigator_bidirectional_moco_;

    // parameters for the moco in kspace binning
    size_t kspace_binning_moco_level_;
    size_t kspace_binning_moco_iter_[MAX_MOCO_LEVEL];
    double kspace_binning_moco_hilbert_strength_;
    double kspace_binning_moco_dissimilarity_sigma_;
    bool  kspace_binning_bidirectional_moco_;

    // whether to perform soft combination
    bool kspace_binning_soft_combination_;

    // navigator signal acceptance window
    double kspace_binning_navigator_window_wide_;
    double kspace_binning_navigator_window_narrow_;

    // method for warpping the complex images ("BSpline", "Linear")
    ISMRMRDINTERP kspace_binning_method_warpping_;

    // whether to exclude the last cardiac cycle for binning
    bool kspace_binning_exclude_last_cardiac_cycle_;

    // some blocks around central kspace must be filled
    size_t kspace_binning_number_of_central_kspace_blocks_;

    // maximal allowed temporal ratio window
    double kspace_binning_max_temporal_window_;

    // temporal ratio window used for binning
    double kspace_binning_temporal_window_;

    // interpolation method to generate best cardiac cycle ('Linear', 'Spline')
    ISMRMRDINTERP kspace_binning_best_cardiac_cycle_interpolator_;

    // recon using certain length of data (if <=0, use the whole data), in the unit of seconds
    double kspace_binning_data_length_used_for_recon_;

    // fill hole with nearest neighbor
    bool kspace_binning_fill_kspace_with_neighbors_;

    // for the flow binning, whether the flow encoding is performed insided every e1
    bool kspace_binning_flow_in_e1_;

    // whether to jointly recon all flow encoding directions
    // if false, every flow encoding direction will be reconed seperately
    bool kspace_binning_flow_recon_jointly_;

    /// --------------
    /// parameters for motion compensated recon
    /// --------------
    size_t motion_comp_num_of_PD_images_;

    // -------------------------------
    // job split
    // -------------------------------
    bool job_split_by_S_;
    size_t job_num_of_N_;
    size_t job_max_Megabytes_;
    size_t job_overlap_;
    // whether to perform computation on the control node
    bool job_perform_on_control_node_;

    // -------------------------------
    // partial fourier handling
    // -------------------------------
    // partial fourier handling algorithms
    ISMRMRDPFALGO partialFourier_algo_;

    // homodyne filter
    // number of iterations
    size_t partialFourier_homodyne_iters_;
    // threshold to stop the iteration
    double partialFourier_homodyne_thres_;
    // density compensation for homodyne filter results
    bool partialFourier_homodyne_densityComp_;

    // POCS
    // number of iterations
    size_t partialFourier_POCS_iters_;
    // threshold to stop the iteration
    double partialFourier_POCS_thres_;
    // transit band width
    size_t partialFourier_POCS_transitBand_;
    // transit band width for E2
    size_t partialFourier_POCS_transitBand_E2_;

    // Feng Huang method
    // kernel size
    size_t partialFourier_FengHuang_kSize_RO_;
    size_t partialFourier_FengHuang_kSize_E1_;
    size_t partialFourier_FengHuang_kSize_E2_;
    // threshold for kernel estimation
    double partialFourier_FengHuang_thresReg_;
    // same kernel for all N
    bool partialFourier_FengHuang_sameKernel_allN_;
    // transit band width
    size_t partialFourier_FengHuang_transitBand_;
    // transit band width for E2
    size_t partialFourier_FengHuang_transitBand_E2_;

    gtPlusReconWorkOrderPara()
    {
        CalibMode_ = ISMRMRD_noacceleration;
        InterleaveDim_ = DIM_NONE;

        acceFactorE1_ = 1;
        acceFactorE2_ = 1;

        kSpaceCenterRO_ = 0;
        kSpaceCenterEncode1_ = 0;
        kSpaceCenterEncode2_ = 0;

        kSpaceMaxRO_ = 1;
        kSpaceMaxEncode1_ = 1;
        kSpaceMaxEncode2_ = 1;

        start_RO_ = -1;
        end_RO_ = -1;

        start_E1_ = -1;
        end_E1_ = -1;

        start_E2_ = -1;
        end_E2_ = -1;

        workFlow_BufferKernel_ = false;
        workFlow_use_BufferedKernel_ = false;

        num_channels_res_ = 1;

        // ----------------------------------------------

        upstream_coil_compression_ = false;
        upstream_coil_compression_thres_ = 1e-3;
        upstream_coil_compression_num_modesKept_ = -1;

        downstream_coil_compression_ = true;
        coil_compression_thres_ = 1e-3;
        coil_compression_num_modesKept_ = -1;

        coil_map_algorithm_ = ISMRMRD_SOUHEIL;
        csm_kSize_ = 7;
        csm_powermethod_num_ = 3;
        csm_true_3D_ = false;
        csm_iter_num_ = 5;
        csm_iter_thres_ = 1e-3;
        csm_use_gpu_ = false;

        // ----------------------------------------------

        recon_algorithm_ = ISMRMRD_GRAPPA;
        recon_auto_parameters_ = true;
        gfactor_needed_ = false;
        wrap_around_map_needed_ = false;

        // ----------------------------------------------

        grappa_kSize_RO_ = 5;
        grappa_kSize_E1_ = 4;
        grappa_kSize_E2_ = 4;
        grappa_reg_lamda_ = 0.0005;
        grappa_calib_over_determine_ratio_ = 0;
        grappa_use_gpu_ = false;

        // ----------------------------------------------

        spirit_kSize_RO_ = 7;
        spirit_kSize_E1_ = 7;
        spirit_kSize_E2_ = 7;

        spirit_oSize_RO_ = 1;
        spirit_oSize_E1_ = 1;
        spirit_oSize_E2_ = 1;

        spirit_reg_lamda_ = 0.005;
        spirit_calib_over_determine_ratio_ = 0;

        spirit_use_gpu_ = false;

        spirit_solve_symmetric_ = false;

        spirit_iter_max_ = 70;
        spirit_iter_thres_ = 1e-5;
        spirit_print_iter_ = false;

        // ----------------------------------------------

        spirit_perform_linear_ = true;
        spirit_perform_grappa_linear_ = false;
        spirit_perform_nonlinear_ = true;

        spirit_parallel_imaging_lamda_ = 1.0;
        spirit_image_reg_lamda_ = 1e-3;
        spirit_data_fidelity_lamda_ = 0;

        spirit_ncg_iter_max_ = 10;
        spirit_ncg_iter_thres_ = 1e-3;
        spirit_ncg_print_iter_ = false;
        spirit_ncg_scale_factor_ = -1.0;

        spirit_slep_iter_max_ = 5;
        spirit_slep_iter_thres_ = 1e-5;
        spirit_slep_print_iter_ = false;
        spirit_slep_keep_third_dimension_coeff_ = false;
        spirit_slep_keep_approx_coeff_ = true;
        spirit_slep_scale_factor_ = -1.0;

        spirit_use_coil_sen_map_ = true;
        spirit_use_moco_enhancement_ = false;
        spirit_recon_moco_images_ = false;

        spirit_RO_enhancement_ratio_ = 1;
        spirit_E1_enhancement_ratio_ = 1;
        spirit_E2_enhancement_ratio_ = 1;
        spirit_temporal_enhancement_ratio_ = 1;

        spirit_2D_scale_per_chunk_ = false;
        spirit_3D_scale_per_chunk_ = true;

        // ----------------------------------------------

        retro_gated_images_ = 0;
        retro_gated_segment_size_ = 0;
        retro_gated_interp_method_ = ISMRMRD_INTERP_RETRO_GATING_BSPLINE;

        // ----------------------------------------------

        kspace_binning_number_of_cardiac_phases_ = 30;
        kspace_binning_minimal_cardiac_phase_width_ = 33; // 33ms, 30 phases for the heart rate of 60

        kspace_binning_multiple_channel_recon_ = true;
        kspace_binning_iterative_non_linear_recon_ = true;
        kspace_binning_iterative_non_linear_recon_slep_ = true;
        kspace_binning_multiple_channel_recon_with_coil_map_ = false;
        kspace_binning_compute_navigator_signal_ = true;

        kspace_binning_navigator_moco_level_ = 4;

        size_t ii;
        for ( ii=0; ii<MAX_MOCO_LEVEL; ii++ ) kspace_binning_navigator_moco_iter_[ii] = 0;
        kspace_binning_navigator_moco_iter_[0] = 1;
        kspace_binning_navigator_moco_iter_[1] = 100;
        kspace_binning_navigator_moco_iter_[2] = 100;
        kspace_binning_navigator_moco_iter_[3] = 100;

        kspace_binning_navigator_hilbert_strength_ = 6.0;
        kspace_binning_navigator_dissimilarity_sigma_ = 2.0;
        kspace_binning_navigator_bidirectional_moco_ = false;

        kspace_binning_moco_level_ = 5;
        for ( ii=0; ii<MAX_MOCO_LEVEL; ii++ ) kspace_binning_moco_iter_[ii] = 0;
        kspace_binning_moco_iter_[0] = 100;
        kspace_binning_moco_iter_[1] = 100;
        kspace_binning_moco_iter_[2] = 100;
        kspace_binning_moco_iter_[3] = 100;
        kspace_binning_moco_iter_[4] = 100;

        kspace_binning_moco_hilbert_strength_ = 12.0;
        kspace_binning_moco_dissimilarity_sigma_ = 2.0;
        kspace_binning_bidirectional_moco_ = false;
        kspace_binning_soft_combination_ = true;
        kspace_binning_navigator_window_wide_ = 0.75;
        kspace_binning_navigator_window_narrow_ = 0.5;
        kspace_binning_method_warpping_ = ISMRMRD_INTERP_BSPLINE;
        kspace_binning_exclude_last_cardiac_cycle_ = false;
        kspace_binning_number_of_central_kspace_blocks_ = 0;
        kspace_binning_max_temporal_window_ = 1.0;
        kspace_binning_temporal_window_ = 4.0;
        kspace_binning_best_cardiac_cycle_interpolator_= ISMRMRD_INTERP_SPLINE;
        kspace_binning_data_length_used_for_recon_ = 0;
        kspace_binning_fill_kspace_with_neighbors_ = false;
        kspace_binning_flow_in_e1_ = true;
        kspace_binning_flow_recon_jointly_ = true;

        // ----------------------------------------------

        motion_comp_num_of_PD_images_ = 0;

        // ----------------------------------------------

        job_split_by_S_ = false;
        job_num_of_N_ = 0;
        job_max_Megabytes_ = 20*1024;
        job_overlap_ = 2;
        job_perform_on_control_node_ = true;

        // ----------------------------------------------

        partialFourier_algo_ = ISMRMRD_PF_ZEROFILLING_FILTER;

        partialFourier_homodyne_iters_ = 6;
        partialFourier_homodyne_thres_ = 1e-2;
        partialFourier_homodyne_densityComp_ = false;

        partialFourier_POCS_iters_ = 6;
        partialFourier_POCS_thres_ = 1e-2;
        partialFourier_POCS_transitBand_ = 16;
        partialFourier_POCS_transitBand_E2_ = 16;

        partialFourier_FengHuang_kSize_RO_ = 5;
        partialFourier_FengHuang_kSize_E1_ = 5;
        partialFourier_FengHuang_kSize_E2_ = 5;
        partialFourier_FengHuang_thresReg_ = 0.005;
        partialFourier_FengHuang_sameKernel_allN_ = false;
        partialFourier_FengHuang_transitBand_ = 16;
        partialFourier_FengHuang_transitBand_E2_ = 16;
    }

    ~gtPlusReconWorkOrderPara() {}
};

template <typename T> 
class gtPlusReconWorkOrder : public gtPlusReconWorkOrderPara
{
public:

    typedef typename realType<T>::Type real_value_type;

    gtPlusReconWorkOrder();
    virtual ~gtPlusReconWorkOrder();

    // reset the status of work order
    // all computed calibration/coil sensitivity results
    // are deleted
    virtual bool reset();

    // check and modify inconsistency in the work order
    virtual bool enforceConsistency(ISMRMRDDIM& /*lastDim*/);

    typedef std::pair<ISMRMRDDIM, size_t> DimensionRecordType;

    // duplicate a workorder without copying the data arrays
    virtual void duplicatePara(gtPlusReconWorkOrderPara& worder) const;
    virtual void duplicate(gtPlusReconWorkOrder<T>& worder) const;

    virtual void copyFromPara(const gtPlusReconWorkOrderPara& worder);

    virtual void printInfo(std::ostream& os) const;
    virtual void print(std::ostream& os) const;

    // -------------------------------
    // input
    // -------------------------------
    // kspace data
    hoNDArray<T> data_;
    // ref data
    hoNDArray<T> ref_;

    // noise data
    hoNDArray<T> noise_;

    // phase correction data
    hoNDArray<T> phaseCorr_;

    // other data
    hoNDArray<T> other_;

    // sometime, the initial kspace can be provided
    hoNDArray<T> kspace_initial_;

    // acqusition time stamp in the unit of second for kspace data lines
    // for the embedded mode, the time stamps of ref lines are also stored
    hoNDArray<real_value_type> time_stamp_;

    // physio time stamp in the unit of second for kspace data lines
    // for the embedded mode, the physio time stamps of ref lines are also stored
    hoNDArray<real_value_type> physio_time_stamp_;

    // dimension starting indexes for the data_
    std::vector< DimensionRecordType > dataDimStartingIndexes_;

    // to support EPI and other trajectories
    // if 1, the readout line is acquired inversely, otherwise, 0
    hoNDArray<unsigned short> reflect_;
    hoNDArray<unsigned short> reflect_ref_;
    hoNDArray<unsigned short> reflect_phaseCorr_;
    hoNDArray<unsigned short> reflect_other_;

    // -------------------------------
    // output
    // -------------------------------
    // reconstructed kspace
    hoNDArray<T> fullkspace_;

    // reconstructed images
    hoNDArray<T> complexIm_;

    // time stamp and physio stamp for reconed images, in the unit of seconds
    // if these fields are not set, the buffered image header will be used
    hoNDArray<real_value_type> recon_time_stamp_;
    hoNDArray<real_value_type> recon_physio_time_stamp_;

    // extra reconstructed results
    // some methods can generate more than one set of reconstruction results
    hoNDArray<T> fullkspace_second_;
    hoNDArray<T> complexIm_second_;
    hoNDArray<real_value_type> recon_time_stamp_second_;
    hoNDArray<real_value_type> recon_physio_time_stamp_second_;

    // gfactor
    hoNDArray<T> gfactor_;

    // wrap-around eig map
    hoNDArray<T> wrap_around_map_;

    // -------------------------------
    // buffers for computation
    // -------------------------------
    // ref for recon
    hoNDArray<T> ref_recon_;
    // ref for coil map
    hoNDArray<T> ref_coil_map_;

    // store the estimated kernel, kernel in image domain
    // if these fields are set before recon, they will be used
    boost::shared_ptr< hoNDArray<T> > kernel_; // [RO E1 srcCHA dstCHA dstE1 1 or N S]
    boost::shared_ptr< hoNDArray<T> > kernelIm_; // [RO E1 srcCHA dstCHA 1 or N S]
    boost::shared_ptr< hoNDArray<T> > unmixingCoeffIm_; // [RO E1 srcCHA 1 or N S]
    boost::shared_ptr< std::vector<hoMatrix<T> > > coilCompressionCoef_; // [dstCHA srcCHA] matrices
    boost::shared_ptr< hoNDArray<T> > coilMap_; // [RO E1 dstCHA 1 or N S]

    // -------------------------------
    // kspace filter for RO/E1/E2 dimension, applied to the reconstruction results
    // -------------------------------
    // 1D filter for kspace data
    hoNDArray<T> filterRO_;
    hoNDArray<T> filterE1_;
    hoNDArray<T> filterE2_;
    // 2D and 3D filter, overwrite the 1D filters
    hoNDArray<T> filterROE1_;
    hoNDArray<T> filterROE1E2_;

    // -------------------------------
    // kspace filter for RO/E1/E2 dimension, applied to the ref data for coil map estimation
    // -------------------------------
    // filter for ref data
    hoNDArray<T> filterRO_ref_;
    hoNDArray<T> filterE1_ref_;
    hoNDArray<T> filterE2_ref_;

    hoNDArray<T> filterROE1_ref_;
    hoNDArray<T> filterROE1E2_ref_;

    // -------------------------------
    // kspace filter for RO/E1/E2 dimension, applied to the data edge in case of partial fourier or asymmetric echo
    // -------------------------------
    // filter for partial fourier/asymmetric echo
    hoNDArray<T> filterRO_partialfourier_;
    hoNDArray<T> filterE1_partialfourier_;
    hoNDArray<T> filterE2_partialfourier_;

    hoNDArray<T> filterROE1_partialfourier_;
    hoNDArray<T> filterROE1E2_partialfourier_;

    // -------------------------------
    // parameters for cloud computing
    // -------------------------------
    bool CloudComputing_;
    unsigned int CloudSize_;

    typedef boost::tuple<std::string, std::string, std::string, unsigned int> CloudNodeType;
    typedef std::vector<CloudNodeType> CloudType;

    CloudType gt_cloud_;
};

template <typename T> 
gtPlusReconWorkOrder<T>::gtPlusReconWorkOrder() : gtPlusReconWorkOrderPara()
{
    hoNDArray<T>* tmp = new hoNDArray<T>();
    kernel_ = boost::shared_ptr< hoNDArray<T> >(tmp);

    tmp = new hoNDArray<T>();
    kernelIm_ = boost::shared_ptr< hoNDArray<T> >(tmp);

    tmp = new hoNDArray<T>();
    unmixingCoeffIm_ = boost::shared_ptr< hoNDArray<T> >(tmp);

    std::vector<hoMatrix<T> >* tmpCoilCoef = new std::vector<hoMatrix<T> >();
    coilCompressionCoef_ = boost::shared_ptr< std::vector<hoMatrix<T> > >(tmpCoilCoef);

    tmp = new hoNDArray<T>();
    coilMap_ = boost::shared_ptr< hoNDArray<T> >(tmp);

    CloudComputing_ = false;
    CloudSize_ = 0;
}

template <typename T> 
gtPlusReconWorkOrder<T>::~gtPlusReconWorkOrder()
{
}

template <typename T> 
bool gtPlusReconWorkOrder<T>::reset()
{
    return true;
}

template <typename T> 
bool gtPlusReconWorkOrder<T>::enforceConsistency(ISMRMRDDIM& /*lastDim*/)
{
    return true;
}

template <typename T> 
void gtPlusReconWorkOrder<T>::duplicatePara(gtPlusReconWorkOrderPara& worder) const
{
    worder.CalibMode_                                  = CalibMode_;
    worder.InterleaveDim_                              = InterleaveDim_;

    worder.acceFactorE1_                               = acceFactorE1_;
    worder.acceFactorE2_                               = acceFactorE2_;

    worder.kSpaceCenterRO_                             = kSpaceCenterRO_;
    worder.kSpaceCenterEncode1_                        = kSpaceCenterEncode1_;
    worder.kSpaceCenterEncode2_                        = kSpaceCenterEncode2_;

    worder.kSpaceMaxRO_                                = kSpaceMaxRO_;
    worder.kSpaceMaxEncode1_                           = kSpaceMaxEncode1_;
    worder.kSpaceMaxEncode2_                           = kSpaceMaxEncode2_;

    worder.workFlow_BufferKernel_                      = workFlow_BufferKernel_;
    worder.workFlow_use_BufferedKernel_                = workFlow_use_BufferedKernel_;
    worder.num_channels_res_                           = num_channels_res_;

    worder.upstream_coil_compression_                  = upstream_coil_compression_;
    worder.upstream_coil_compression_thres_            = upstream_coil_compression_thres_;
    worder.upstream_coil_compression_num_modesKept_    = upstream_coil_compression_num_modesKept_;

    worder.downstream_coil_compression_                = downstream_coil_compression_;
    worder.coil_compression_thres_                     = coil_compression_thres_;
    worder.coil_compression_num_modesKept_             = coil_compression_num_modesKept_;

    worder.coil_map_algorithm_                         = coil_map_algorithm_;
    worder.csm_kSize_                                  = csm_kSize_;
    worder.csm_powermethod_num_                        = csm_powermethod_num_;
    worder.csm_true_3D_                                = csm_true_3D_;
    worder.csm_iter_num_                               = csm_iter_num_;
    worder.csm_iter_thres_                             = csm_iter_thres_;
    worder.csm_use_gpu_                                = csm_use_gpu_;

    worder.start_RO_                                   = start_RO_;
    worder.end_RO_                                     = end_RO_;

    worder.start_E1_                                   = start_E1_;
    worder.end_E1_                                     = end_E1_;

    worder.start_E2_                                   = start_E2_;
    worder.end_E2_                                     = end_E2_;

    worder.recon_algorithm_                            = recon_algorithm_;
    worder.recon_auto_parameters_                      = recon_auto_parameters_;
    worder.gfactor_needed_                             = gfactor_needed_;
    worder.wrap_around_map_needed_                     = wrap_around_map_needed_;

    worder.grappa_kSize_RO_                            = grappa_kSize_RO_;
    worder.grappa_kSize_RO_                            = grappa_kSize_RO_;
    worder.grappa_kSize_E1_                            = grappa_kSize_E1_;
    worder.grappa_kSize_E2_                            = grappa_kSize_E2_;
    worder.grappa_reg_lamda_                           = grappa_reg_lamda_;
    worder.grappa_calib_over_determine_ratio_          = grappa_calib_over_determine_ratio_;
    worder.grappa_use_gpu_                             = grappa_use_gpu_;

    worder.spirit_kSize_RO_                            = spirit_kSize_RO_;
    worder.spirit_kSize_E1_                            = spirit_kSize_E1_;
    worder.spirit_kSize_E2_                            = spirit_kSize_E2_;
    worder.spirit_oSize_RO_                            = spirit_oSize_RO_;
    worder.spirit_oSize_E1_                            = spirit_oSize_E1_;
    worder.spirit_oSize_E2_                            = spirit_oSize_E2_;
    worder.spirit_reg_lamda_                           = spirit_reg_lamda_;
    worder.spirit_use_gpu_                             = spirit_use_gpu_;
    worder.spirit_calib_over_determine_ratio_          = spirit_calib_over_determine_ratio_;
    worder.spirit_solve_symmetric_                     = spirit_solve_symmetric_;
    worder.spirit_iter_max_                            = spirit_iter_max_;
    worder.spirit_iter_thres_                          = spirit_iter_thres_;
    worder.spirit_print_iter_                          = spirit_print_iter_;

    worder.spirit_perform_linear_                      = spirit_perform_linear_;
    worder.spirit_perform_grappa_linear_               = spirit_perform_grappa_linear_;
    worder.spirit_perform_nonlinear_                   = spirit_perform_nonlinear_;
    worder.spirit_parallel_imaging_lamda_              = spirit_parallel_imaging_lamda_;
    worder.spirit_image_reg_lamda_                     = spirit_image_reg_lamda_;
    worder.spirit_data_fidelity_lamda_                 = spirit_data_fidelity_lamda_;
    worder.spirit_ncg_iter_max_                        = spirit_ncg_iter_max_;
    worder.spirit_ncg_iter_thres_                      = spirit_ncg_iter_thres_;
    worder.spirit_ncg_scale_factor_                    = spirit_ncg_scale_factor_;
    worder.spirit_ncg_print_iter_                      = spirit_ncg_print_iter_;
    worder.spirit_slep_iter_max_                       = spirit_slep_iter_max_;
    worder.spirit_slep_iter_thres_                     = spirit_slep_iter_thres_;
    worder.spirit_slep_print_iter_                     = spirit_slep_print_iter_;
    worder.spirit_slep_keep_third_dimension_coeff_     = spirit_slep_keep_third_dimension_coeff_;
    worder.spirit_slep_keep_approx_coeff_              = spirit_slep_keep_approx_coeff_;
    worder.spirit_slep_scale_factor_                   = spirit_slep_scale_factor_;
    worder.spirit_use_coil_sen_map_                    = spirit_use_coil_sen_map_;
    worder.spirit_use_moco_enhancement_                = spirit_use_moco_enhancement_;
    worder.spirit_recon_moco_images_                   = spirit_recon_moco_images_;
    worder.spirit_RO_enhancement_ratio_                = spirit_RO_enhancement_ratio_;
    worder.spirit_E1_enhancement_ratio_                = spirit_E1_enhancement_ratio_;
    worder.spirit_E2_enhancement_ratio_                = spirit_E2_enhancement_ratio_;
    worder.spirit_temporal_enhancement_ratio_          = spirit_temporal_enhancement_ratio_;
    worder.spirit_2D_scale_per_chunk_                  = spirit_2D_scale_per_chunk_;
    worder.spirit_3D_scale_per_chunk_                  = spirit_3D_scale_per_chunk_;

    worder.retro_gated_images_                         = retro_gated_images_;
    worder.retro_gated_segment_size_                   = retro_gated_segment_size_;
    worder.retro_gated_interp_method_                  = retro_gated_interp_method_;

    worder.kspace_binning_number_of_cardiac_phases_                 = kspace_binning_number_of_cardiac_phases_;
    worder.kspace_binning_minimal_cardiac_phase_width_              = kspace_binning_minimal_cardiac_phase_width_;
    worder.kspace_binning_multiple_channel_recon_                   = kspace_binning_multiple_channel_recon_;
    worder.kspace_binning_iterative_non_linear_recon_               = kspace_binning_iterative_non_linear_recon_;
    worder.kspace_binning_iterative_non_linear_recon_slep_          = kspace_binning_iterative_non_linear_recon_slep_;
    worder.kspace_binning_multiple_channel_recon_with_coil_map_     = kspace_binning_multiple_channel_recon_with_coil_map_;
    worder.kspace_binning_compute_navigator_signal_                 = kspace_binning_compute_navigator_signal_;
    worder.kspace_binning_navigator_moco_level_                     = kspace_binning_navigator_moco_level_;
    memcpy(worder.kspace_binning_navigator_moco_iter_, kspace_binning_navigator_moco_iter_, sizeof(size_t)*MAX_MOCO_LEVEL);
    worder.kspace_binning_navigator_hilbert_strength_               = kspace_binning_navigator_hilbert_strength_;
    worder.kspace_binning_navigator_dissimilarity_sigma_            = kspace_binning_navigator_dissimilarity_sigma_;
    worder.kspace_binning_navigator_bidirectional_moco_             = kspace_binning_navigator_bidirectional_moco_;
    worder.kspace_binning_moco_level_                               = kspace_binning_moco_level_;
    memcpy(worder.kspace_binning_moco_iter_, kspace_binning_moco_iter_, sizeof(size_t)*MAX_MOCO_LEVEL);
    worder.kspace_binning_moco_hilbert_strength_                    = kspace_binning_moco_hilbert_strength_;
    worder.kspace_binning_moco_dissimilarity_sigma_                 = kspace_binning_moco_dissimilarity_sigma_;
    worder.kspace_binning_bidirectional_moco_                       = kspace_binning_bidirectional_moco_;
    worder.kspace_binning_soft_combination_                         = kspace_binning_soft_combination_;
    worder.kspace_binning_navigator_window_wide_                    = kspace_binning_navigator_window_wide_;
    worder.kspace_binning_navigator_window_narrow_                  = kspace_binning_navigator_window_narrow_;
    worder.kspace_binning_method_warpping_                          = kspace_binning_method_warpping_;
    worder.kspace_binning_exclude_last_cardiac_cycle_               = kspace_binning_exclude_last_cardiac_cycle_;
    worder.kspace_binning_number_of_central_kspace_blocks_          = kspace_binning_number_of_central_kspace_blocks_;
    worder.kspace_binning_max_temporal_window_                      = kspace_binning_max_temporal_window_;
    worder.kspace_binning_temporal_window_                          = kspace_binning_temporal_window_;
    worder.kspace_binning_best_cardiac_cycle_interpolator_          = kspace_binning_best_cardiac_cycle_interpolator_;
    worder.kspace_binning_data_length_used_for_recon_               = kspace_binning_data_length_used_for_recon_;
    worder.kspace_binning_fill_kspace_with_neighbors_               = kspace_binning_fill_kspace_with_neighbors_;
    worder.kspace_binning_flow_in_e1_                               = kspace_binning_flow_in_e1_;
    worder.kspace_binning_flow_recon_jointly_                       = kspace_binning_flow_recon_jointly_;

    worder.motion_comp_num_of_PD_images_                            = motion_comp_num_of_PD_images_;

    worder.job_split_by_S_                             = job_split_by_S_;
    worder.job_num_of_N_                               = job_num_of_N_;
    worder.job_max_Megabytes_                          = job_max_Megabytes_;
    worder.job_overlap_                                = job_overlap_;
    worder.job_perform_on_control_node_                = job_perform_on_control_node_;

    worder.partialFourier_algo_                        = partialFourier_algo_;

    worder.partialFourier_homodyne_iters_              = partialFourier_homodyne_iters_;
    worder.partialFourier_homodyne_thres_              = partialFourier_homodyne_thres_;
    worder.partialFourier_homodyne_densityComp_        = partialFourier_homodyne_densityComp_;

    worder.partialFourier_POCS_iters_                  = partialFourier_POCS_iters_;
    worder.partialFourier_POCS_thres_                  = partialFourier_POCS_thres_;
    worder.partialFourier_POCS_transitBand_            = partialFourier_POCS_transitBand_;
    worder.partialFourier_POCS_transitBand_E2_         = partialFourier_POCS_transitBand_E2_;

    worder.partialFourier_FengHuang_kSize_RO_          = partialFourier_FengHuang_kSize_RO_;
    worder.partialFourier_FengHuang_kSize_E1_          = partialFourier_FengHuang_kSize_E1_;
    worder.partialFourier_FengHuang_kSize_E2_          = partialFourier_FengHuang_kSize_E2_;
    worder.partialFourier_FengHuang_thresReg_          = partialFourier_FengHuang_thresReg_;
    worder.partialFourier_FengHuang_sameKernel_allN_   = partialFourier_FengHuang_sameKernel_allN_;
    worder.partialFourier_FengHuang_transitBand_       = partialFourier_FengHuang_transitBand_;
    worder.partialFourier_FengHuang_transitBand_E2_    = partialFourier_FengHuang_transitBand_E2_;
}

template <typename T> 
void gtPlusReconWorkOrder<T>::duplicate(gtPlusReconWorkOrder<T>& worder) const
{
    this->duplicatePara(worder);

    worder.dataDimStartingIndexes_      = dataDimStartingIndexes_;

    worder.filterRO_                    = filterRO_;
    worder.filterE1_                    = filterE1_;
    worder.filterE2_                    = filterE2_;
    worder.filterROE1_                  = filterROE1_;
    worder.filterROE1E2_                = filterROE1E2_;

    worder.filterRO_ref_                = filterRO_ref_;
    worder.filterE1_ref_                = filterE1_ref_;
    worder.filterE2_ref_                = filterE2_ref_;
    worder.filterROE1_ref_              = filterROE1_ref_;
    worder.filterROE1E2_ref_            = filterROE1E2_ref_;

    worder.filterRO_partialfourier_     = filterRO_partialfourier_;
    worder.filterE1_partialfourier_     = filterE1_partialfourier_;
    worder.filterE2_partialfourier_     = filterE2_partialfourier_;
    worder.filterROE1_partialfourier_   = filterROE1_partialfourier_;
    worder.filterROE1E2_partialfourier_ = filterROE1E2_partialfourier_;

    worder.CloudComputing_              = CloudComputing_;
    worder.CloudSize_                   = CloudSize_;
    worder.gt_cloud_                    = gt_cloud_;
}

template <typename T> 
void gtPlusReconWorkOrder<T>::copyFromPara(const gtPlusReconWorkOrderPara& worder)
{
    CalibMode_                                  = worder.CalibMode_;
    InterleaveDim_                              = worder.InterleaveDim_;

    acceFactorE1_                               = worder.acceFactorE1_;
    acceFactorE2_                               = worder.acceFactorE2_;

    kSpaceCenterRO_                             = worder.kSpaceCenterRO_;
    kSpaceCenterEncode1_                        = worder.kSpaceCenterEncode1_;
    kSpaceCenterEncode2_                        = worder.kSpaceCenterEncode2_;

    kSpaceMaxRO_                                = worder.kSpaceMaxRO_;
    kSpaceMaxEncode1_                           = worder.kSpaceMaxEncode1_;
    kSpaceMaxEncode2_                           = worder.kSpaceMaxEncode2_;

    workFlow_BufferKernel_                      = worder.workFlow_BufferKernel_;
    workFlow_use_BufferedKernel_                = worder.workFlow_use_BufferedKernel_;
    num_channels_res_                           = worder.num_channels_res_;

    upstream_coil_compression_                  = worder.upstream_coil_compression_;
    upstream_coil_compression_thres_            = worder.upstream_coil_compression_thres_;
    upstream_coil_compression_num_modesKept_    = worder.upstream_coil_compression_num_modesKept_;

    downstream_coil_compression_                = worder.downstream_coil_compression_;
    coil_compression_thres_                     = worder.coil_compression_thres_;
    coil_compression_num_modesKept_             = worder.coil_compression_num_modesKept_;

    coil_map_algorithm_                         = worder.coil_map_algorithm_;
    csm_kSize_                                  = worder.csm_kSize_;
    csm_powermethod_num_                        = worder.csm_powermethod_num_;
    csm_true_3D_                                = worder.csm_true_3D_;
    csm_iter_num_                               = worder.csm_iter_num_;
    csm_iter_thres_                             = worder.csm_iter_thres_;
    csm_use_gpu_                                = worder.csm_use_gpu_;

    start_RO_                                   = worder.start_RO_;
    end_RO_                                     = worder.end_RO_;

    start_E1_                                   = worder.start_E1_;
    end_E1_                                     = worder.end_E1_;

    start_E2_                                   = worder.start_E2_;
    end_E2_                                     = worder.end_E2_;

    recon_algorithm_                            = worder.recon_algorithm_;
    recon_auto_parameters_                      = worder.recon_auto_parameters_;
    gfactor_needed_                             = worder.gfactor_needed_;
    wrap_around_map_needed_                     = worder.wrap_around_map_needed_;

    grappa_kSize_RO_                            = worder.grappa_kSize_RO_;
    grappa_kSize_RO_                            = worder.grappa_kSize_RO_;
    grappa_kSize_E1_                            = worder.grappa_kSize_E1_;
    grappa_kSize_E2_                            = worder.grappa_kSize_E2_;
    grappa_reg_lamda_                           = worder.grappa_reg_lamda_;
    grappa_calib_over_determine_ratio_          = worder.grappa_calib_over_determine_ratio_;
    grappa_use_gpu_                             = worder.grappa_use_gpu_;

    spirit_kSize_RO_                            = worder.spirit_kSize_RO_;
    spirit_kSize_E1_                            = worder.spirit_kSize_E1_;
    spirit_kSize_E2_                            = worder.spirit_kSize_E2_;
    spirit_oSize_RO_                            = worder.spirit_oSize_RO_;
    spirit_oSize_E1_                            = worder.spirit_oSize_E1_;
    spirit_oSize_E2_                            = worder.spirit_oSize_E2_;
    spirit_reg_lamda_                           = worder.spirit_reg_lamda_;
    spirit_use_gpu_                             = worder.spirit_use_gpu_;
    spirit_calib_over_determine_ratio_          = worder.spirit_calib_over_determine_ratio_;
    spirit_solve_symmetric_                     = worder.spirit_solve_symmetric_;
    spirit_iter_max_                            = worder.spirit_iter_max_;
    spirit_iter_thres_                          = worder.spirit_iter_thres_;
    spirit_print_iter_                          = worder.spirit_print_iter_;

    spirit_perform_linear_                      = worder.spirit_perform_linear_;
    spirit_perform_grappa_linear_               = worder.spirit_perform_grappa_linear_;
    spirit_perform_nonlinear_                   = worder.spirit_perform_nonlinear_;
    spirit_parallel_imaging_lamda_              = worder.spirit_parallel_imaging_lamda_;
    spirit_image_reg_lamda_                     = worder.spirit_image_reg_lamda_;
    spirit_data_fidelity_lamda_                 = worder.spirit_data_fidelity_lamda_;
    spirit_ncg_iter_max_                        = worder.spirit_ncg_iter_max_;
    spirit_ncg_iter_thres_                      = worder.spirit_ncg_iter_thres_;
    spirit_ncg_scale_factor_                    = worder.spirit_ncg_scale_factor_;
    spirit_ncg_print_iter_                      = worder.spirit_ncg_print_iter_;
    spirit_slep_iter_max_                       = worder.spirit_slep_iter_max_;
    spirit_slep_iter_thres_                     = worder.spirit_slep_iter_thres_;
    spirit_slep_print_iter_                     = worder.spirit_slep_print_iter_;
    spirit_slep_keep_third_dimension_coeff_     = worder.spirit_slep_keep_third_dimension_coeff_;
    spirit_slep_keep_approx_coeff_              = worder.spirit_slep_keep_approx_coeff_;
    spirit_slep_scale_factor_                   = worder.spirit_slep_scale_factor_;
    spirit_use_coil_sen_map_                    = worder.spirit_use_coil_sen_map_;
    spirit_use_moco_enhancement_                = worder.spirit_use_moco_enhancement_;
    spirit_recon_moco_images_                   = worder.spirit_recon_moco_images_;
    spirit_RO_enhancement_ratio_                = worder.spirit_RO_enhancement_ratio_;
    spirit_E1_enhancement_ratio_                = worder.spirit_E1_enhancement_ratio_;
    spirit_E2_enhancement_ratio_                = worder.spirit_E2_enhancement_ratio_;
    spirit_temporal_enhancement_ratio_          = worder.spirit_temporal_enhancement_ratio_;
    spirit_2D_scale_per_chunk_                  = worder.spirit_2D_scale_per_chunk_;
    spirit_3D_scale_per_chunk_                  = worder.spirit_3D_scale_per_chunk_;

    retro_gated_images_                         = worder.retro_gated_images_;
    retro_gated_segment_size_                   = worder.retro_gated_segment_size_;
    retro_gated_interp_method_                  = worder.retro_gated_interp_method_;

    kspace_binning_number_of_cardiac_phases_          = worder.kspace_binning_number_of_cardiac_phases_;
    kspace_binning_minimal_cardiac_phase_width_          = worder.kspace_binning_minimal_cardiac_phase_width_;
    kspace_binning_multiple_channel_recon_         = worder.kspace_binning_multiple_channel_recon_;
    kspace_binning_iterative_non_linear_recon_              = worder.kspace_binning_iterative_non_linear_recon_;
    kspace_binning_iterative_non_linear_recon_slep_              = worder.kspace_binning_iterative_non_linear_recon_slep_;
    kspace_binning_multiple_channel_recon_with_coil_map_ = worder.kspace_binning_multiple_channel_recon_with_coil_map_;
    kspace_binning_compute_navigator_signal_       = worder.kspace_binning_compute_navigator_signal_;
    kspace_binning_navigator_moco_level_                = worder.kspace_binning_navigator_moco_level_;
    memcpy(kspace_binning_navigator_moco_iter_, worder.kspace_binning_navigator_moco_iter_, sizeof(size_t)*MAX_MOCO_LEVEL);
    kspace_binning_navigator_hilbert_strength_                    = worder.kspace_binning_navigator_hilbert_strength_;
    kspace_binning_navigator_dissimilarity_sigma_                 = worder.kspace_binning_navigator_dissimilarity_sigma_;
    kspace_binning_navigator_bidirectional_moco_         = worder.kspace_binning_navigator_bidirectional_moco_;
    kspace_binning_moco_level_                   = worder.kspace_binning_moco_level_;
    memcpy(kspace_binning_moco_iter_, worder.kspace_binning_moco_iter_, sizeof(size_t)*MAX_MOCO_LEVEL);
    kspace_binning_moco_hilbert_strength_                       = worder.kspace_binning_moco_hilbert_strength_;
    kspace_binning_moco_dissimilarity_sigma_                    = worder.kspace_binning_moco_dissimilarity_sigma_;
    kspace_binning_bidirectional_moco_            = worder.kspace_binning_bidirectional_moco_;
    kspace_binning_soft_combination_            = worder.kspace_binning_soft_combination_;
    kspace_binning_navigator_window_wide_                  = worder.kspace_binning_navigator_window_wide_;
    kspace_binning_navigator_window_narrow_                = worder.kspace_binning_navigator_window_narrow_;
    kspace_binning_method_warpping_              = worder.kspace_binning_method_warpping_;
    kspace_binning_exclude_last_cardiac_cycle_            = worder.kspace_binning_exclude_last_cardiac_cycle_;
    kspace_binning_number_of_central_kspace_blocks_         = worder.kspace_binning_number_of_central_kspace_blocks_;
    kspace_binning_max_temporal_window_    = worder.kspace_binning_max_temporal_window_;
    kspace_binning_temporal_window_    = worder.kspace_binning_temporal_window_;
    kspace_binning_best_cardiac_cycle_interpolator_        = worder.kspace_binning_best_cardiac_cycle_interpolator_;
    kspace_binning_data_length_used_for_recon_            = worder.kspace_binning_data_length_used_for_recon_;
    kspace_binning_fill_kspace_with_neighbors_ = worder.kspace_binning_fill_kspace_with_neighbors_;
    kspace_binning_flow_in_e1_ = worder.kspace_binning_flow_in_e1_;
    kspace_binning_flow_recon_jointly_ = worder.kspace_binning_flow_recon_jointly_;

    motion_comp_num_of_PD_images_ = worder.motion_comp_num_of_PD_images_;

    job_split_by_S_                             = worder.job_split_by_S_;
    job_num_of_N_                               = worder.job_num_of_N_;
    job_max_Megabytes_                          = worder.job_max_Megabytes_;
    job_overlap_                                = worder.job_overlap_;
    job_perform_on_control_node_                = worder.job_perform_on_control_node_;

    partialFourier_algo_                        = worder.partialFourier_algo_;

    partialFourier_homodyne_iters_              = worder.partialFourier_homodyne_iters_;
    partialFourier_homodyne_thres_              = worder.partialFourier_homodyne_thres_;
    partialFourier_homodyne_densityComp_        = worder.partialFourier_homodyne_densityComp_;

    partialFourier_POCS_iters_                  = worder.partialFourier_POCS_iters_;
    partialFourier_POCS_thres_                  = worder.partialFourier_POCS_thres_;
    partialFourier_POCS_transitBand_            = worder.partialFourier_POCS_transitBand_;
    partialFourier_POCS_transitBand_E2_         = worder.partialFourier_POCS_transitBand_E2_;

    partialFourier_FengHuang_kSize_RO_          = worder.partialFourier_FengHuang_kSize_RO_;
    partialFourier_FengHuang_kSize_E1_          = worder.partialFourier_FengHuang_kSize_E1_;
    partialFourier_FengHuang_kSize_E2_          = worder.partialFourier_FengHuang_kSize_E2_;
    partialFourier_FengHuang_thresReg_          = worder.partialFourier_FengHuang_thresReg_;
    partialFourier_FengHuang_sameKernel_allN_   = worder.partialFourier_FengHuang_sameKernel_allN_;
    partialFourier_FengHuang_transitBand_       = worder.partialFourier_FengHuang_transitBand_;
    partialFourier_FengHuang_transitBand_E2_    = worder.partialFourier_FengHuang_transitBand_E2_;
}

template <typename T> 
void gtPlusReconWorkOrder<T>::printInfo(std::ostream& os) const
{
    using namespace std;
    GADGET_PARA_PRINT(CalibMode_);
    GADGET_PARA_PRINT(InterleaveDim_);
    GADGET_PARA_PRINT(acceFactorE1_);
    GADGET_PARA_PRINT(acceFactorE2_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(kSpaceCenterRO_);
    GADGET_PARA_PRINT(kSpaceCenterEncode1_);
    GADGET_PARA_PRINT(kSpaceCenterEncode2_);
    GADGET_PARA_PRINT(kSpaceMaxRO_);
    GADGET_PARA_PRINT(kSpaceMaxEncode1_);
    GADGET_PARA_PRINT(kSpaceMaxEncode2_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(workFlow_BufferKernel_);
    GADGET_PARA_PRINT(workFlow_use_BufferedKernel_);
    GADGET_PARA_PRINT(num_channels_res_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(upstream_coil_compression_);
    GADGET_PARA_PRINT(upstream_coil_compression_thres_);
    GADGET_PARA_PRINT(upstream_coil_compression_num_modesKept_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(downstream_coil_compression_);
    GADGET_PARA_PRINT(coil_compression_thres_);
    GADGET_PARA_PRINT(coil_compression_num_modesKept_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(coil_map_algorithm_);
    GADGET_PARA_PRINT(csm_kSize_);
    GADGET_PARA_PRINT(csm_powermethod_num_);
    GADGET_PARA_PRINT(csm_true_3D_);
    GADGET_PARA_PRINT(csm_iter_num_);
    GADGET_PARA_PRINT(csm_iter_thres_);
    GADGET_PARA_PRINT(csm_use_gpu_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(start_RO_);
    GADGET_PARA_PRINT(end_RO_);
    GADGET_PARA_PRINT(start_E1_);
    GADGET_PARA_PRINT(end_E1_);
    GADGET_PARA_PRINT(start_E2_);
    GADGET_PARA_PRINT(end_E2_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(recon_algorithm_);
    GADGET_PARA_PRINT(recon_auto_parameters_);
    GADGET_PARA_PRINT(gfactor_needed_);
    GADGET_PARA_PRINT(wrap_around_map_needed_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(grappa_kSize_RO_);
    GADGET_PARA_PRINT(grappa_kSize_E1_);
    GADGET_PARA_PRINT(grappa_kSize_E2_);
    GADGET_PARA_PRINT(grappa_reg_lamda_);
    GADGET_PARA_PRINT(grappa_calib_over_determine_ratio_);
    GADGET_PARA_PRINT(grappa_use_gpu_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(spirit_kSize_RO_);
    GADGET_PARA_PRINT(spirit_kSize_E1_);
    GADGET_PARA_PRINT(spirit_kSize_E2_);
    GADGET_PARA_PRINT(spirit_oSize_RO_);
    GADGET_PARA_PRINT(spirit_oSize_E1_);
    GADGET_PARA_PRINT(spirit_oSize_E2_);
    GADGET_PARA_PRINT(spirit_reg_lamda_);
    GADGET_PARA_PRINT(spirit_use_gpu_);
    GADGET_PARA_PRINT(spirit_calib_over_determine_ratio_);
    GADGET_PARA_PRINT(spirit_solve_symmetric_);
    GADGET_PARA_PRINT(spirit_iter_max_);
    GADGET_PARA_PRINT(spirit_iter_thres_);
    GADGET_PARA_PRINT(spirit_print_iter_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(spirit_perform_linear_);
    GADGET_PARA_PRINT(spirit_perform_grappa_linear_);
    GADGET_PARA_PRINT(spirit_perform_nonlinear_);
    GADGET_PARA_PRINT(spirit_parallel_imaging_lamda_);
    GADGET_PARA_PRINT(spirit_image_reg_lamda_);
    GADGET_PARA_PRINT(spirit_data_fidelity_lamda_);
    GADGET_PARA_PRINT(spirit_ncg_iter_max_);
    GADGET_PARA_PRINT(spirit_ncg_iter_thres_);
    GADGET_PARA_PRINT(spirit_ncg_scale_factor_);
    GADGET_PARA_PRINT(spirit_ncg_print_iter_);
    GADGET_PARA_PRINT(spirit_slep_iter_max_);
    GADGET_PARA_PRINT(spirit_slep_iter_thres_);
    GADGET_PARA_PRINT(spirit_slep_print_iter_);
    GADGET_PARA_PRINT(spirit_slep_keep_third_dimension_coeff_);
    GADGET_PARA_PRINT(spirit_slep_keep_approx_coeff_);
    GADGET_PARA_PRINT(spirit_slep_scale_factor_);
    GADGET_PARA_PRINT(spirit_use_coil_sen_map_);
    GADGET_PARA_PRINT(spirit_use_moco_enhancement_);
    GADGET_PARA_PRINT(spirit_recon_moco_images_);
    GADGET_PARA_PRINT(spirit_RO_enhancement_ratio_);
    GADGET_PARA_PRINT(spirit_E1_enhancement_ratio_);
    GADGET_PARA_PRINT(spirit_E2_enhancement_ratio_);
    GADGET_PARA_PRINT(spirit_temporal_enhancement_ratio_);
    GADGET_PARA_PRINT(spirit_2D_scale_per_chunk_);
    GADGET_PARA_PRINT(spirit_3D_scale_per_chunk_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(retro_gated_images_);
    GADGET_PARA_PRINT(retro_gated_segment_size_);
    GADGET_PARA_PRINT(retro_gated_interp_method_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(kspace_binning_number_of_cardiac_phases_);
    GADGET_PARA_PRINT(kspace_binning_minimal_cardiac_phase_width_);
    GADGET_PARA_PRINT(kspace_binning_multiple_channel_recon_);
    GADGET_PARA_PRINT(kspace_binning_iterative_non_linear_recon_);
    GADGET_PARA_PRINT(kspace_binning_iterative_non_linear_recon_slep_);
    GADGET_PARA_PRINT(kspace_binning_multiple_channel_recon_with_coil_map_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(kspace_binning_compute_navigator_signal_);

    GADGET_PARA_PRINT(kspace_binning_navigator_moco_level_);
    std::stringstream ostr;
    ostr << " [ ";
    size_t ii;
    for ( ii=0; ii<kspace_binning_navigator_moco_level_; ii++ )
    {
        ostr << kspace_binning_navigator_moco_iter_[ii] << " ";
    }
    ostr << " ] " << std::endl;
    GDEBUG_STREAM(ostr.str());

    GADGET_PARA_PRINT(kspace_binning_navigator_hilbert_strength_);
    GADGET_PARA_PRINT(kspace_binning_navigator_dissimilarity_sigma_);
    GADGET_PARA_PRINT(kspace_binning_navigator_bidirectional_moco_);
    GDEBUG_STREAM("---------------------");

    GADGET_PARA_PRINT(kspace_binning_moco_level_);
    std::stringstream ostr_moco_level;
    ostr_moco_level << " [ ";
    for ( ii=0; ii<kspace_binning_moco_level_; ii++ )
    {
        ostr_moco_level << kspace_binning_moco_iter_[ii] << " ";
    }
    ostr_moco_level << " ] " << std::endl;
    GDEBUG_STREAM(ostr_moco_level.str());

    GADGET_PARA_PRINT(kspace_binning_moco_hilbert_strength_);
    GADGET_PARA_PRINT(kspace_binning_moco_dissimilarity_sigma_);
    GADGET_PARA_PRINT(kspace_binning_bidirectional_moco_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(kspace_binning_soft_combination_);
    GADGET_PARA_PRINT(kspace_binning_navigator_window_wide_);
    GADGET_PARA_PRINT(kspace_binning_navigator_window_narrow_);
    GADGET_PARA_PRINT(kspace_binning_method_warpping_);
    GADGET_PARA_PRINT(kspace_binning_exclude_last_cardiac_cycle_);
    GADGET_PARA_PRINT(kspace_binning_number_of_central_kspace_blocks_);
    GADGET_PARA_PRINT(kspace_binning_max_temporal_window_);
    GADGET_PARA_PRINT(kspace_binning_temporal_window_);
    GADGET_PARA_PRINT(kspace_binning_best_cardiac_cycle_interpolator_);
    GADGET_PARA_PRINT(kspace_binning_data_length_used_for_recon_);
    GADGET_PARA_PRINT(kspace_binning_fill_kspace_with_neighbors_);
    GADGET_PARA_PRINT(kspace_binning_flow_in_e1_);
    GADGET_PARA_PRINT(kspace_binning_flow_recon_jointly_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(motion_comp_num_of_PD_images_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(job_split_by_S_);
    GADGET_PARA_PRINT(job_num_of_N_);
    GADGET_PARA_PRINT(job_max_Megabytes_);
    GADGET_PARA_PRINT(job_overlap_);
    GADGET_PARA_PRINT(job_perform_on_control_node_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(partialFourier_algo_);
    GADGET_PARA_PRINT(partialFourier_homodyne_iters_);
    GADGET_PARA_PRINT(partialFourier_homodyne_thres_);
    GADGET_PARA_PRINT(partialFourier_homodyne_densityComp_);
    GADGET_PARA_PRINT(partialFourier_POCS_iters_);
    GADGET_PARA_PRINT(partialFourier_POCS_thres_);
    GADGET_PARA_PRINT(partialFourier_POCS_transitBand_);
    GADGET_PARA_PRINT(partialFourier_POCS_transitBand_E2_);
    GADGET_PARA_PRINT(partialFourier_FengHuang_kSize_RO_);
    GADGET_PARA_PRINT(partialFourier_FengHuang_kSize_E1_);
    GADGET_PARA_PRINT(partialFourier_FengHuang_kSize_E2_);
    GADGET_PARA_PRINT(partialFourier_FengHuang_thresReg_);
    GADGET_PARA_PRINT(partialFourier_FengHuang_sameKernel_allN_);
    GADGET_PARA_PRINT(partialFourier_FengHuang_transitBand_);
    GADGET_PARA_PRINT(partialFourier_FengHuang_transitBand_E2_);
    GDEBUG_STREAM("---------------------");
    GADGET_PARA_PRINT(CloudComputing_);
    GADGET_PARA_PRINT(CloudSize_);
    for ( unsigned int nn=0; nn<gt_cloud_.size(); nn++ )
    {
        GADGET_PARA_PRINT(gt_cloud_[nn]);
    }
}

template <typename T> 
void gtPlusReconWorkOrder<T>::print(std::ostream& os) const
{
    using namespace std;
    os << "-------------- gtPlusReconWorkOrder ---------------" << endl;
    printInfo(os);
    os << "---------------------------------------------------" << endl;
}

}}

#include "gtPlusISMRMRDReconWorkOrder2DT.h"
#include "gtPlusISMRMRDReconWorkOrder3DT.h"
