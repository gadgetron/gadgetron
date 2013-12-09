/** \file   gtPlusISMRMRDReconWorkOrder.h
    \brief  Define the GtPlus reconstruction workorder and parameters
    \author Hui Xue
*/

#pragma once

#include "ismrmrd.h"
#include "gtPlusISMRMRDReconUtil.h"

namespace Gadgetron { namespace gtPlus {

struct gtPlusReconWorkOrderPara
{
    ISMRMRDCALIBMODE CalibMode_;
    ISMRMRDDIM InterleaveDim_;

    // acceleration factor along E1 and E2
    double acceFactorE1_;
    double acceFactorE2_;

    // kspace center for RO/E1/E2
    unsigned long long kSpaceCenterRO_;
    unsigned long long kSpaceCenterEncode1_;
    unsigned long long kSpaceCenterEncode2_;

    // kspace max acquired number for RO/E1/E2
    unsigned long long kSpaceMaxRO_;
    unsigned long long kSpaceMaxEncode1_;
    unsigned long long kSpaceMaxEncode2_;

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
    unsigned long long num_channels_res_;

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
    Gadgetron::gtPlus::ISMRMRDCOILMAPALGO coil_map_algorithm_;

    // for ISMRMRD_SOUHEIL
    unsigned long long csm_kSize_;
    unsigned long long csm_powermethod_num_;

    // for ISMRMRD_SOUHEIL_ITER
    unsigned long long csm_iter_num_;
    double csm_iter_thres_;

    // -------------------------------
    // parameters for variant reconstruction algorithms
    // -------------------------------
    Gadgetron::gtPlus::ISMRMRDALGO recon_algorithm_;
    bool recon_auto_parameters_;

    // grappa
    unsigned long long grappa_kSize_RO_;
    unsigned long long grappa_kSize_E1_;
    unsigned long long grappa_kSize_E2_;
    double grappa_reg_lamda_;

    // sense

    // soft sense

    // SPIRiT
    unsigned long long spirit_kSize_RO_;
    unsigned long long spirit_kSize_E1_;
    unsigned long long spirit_kSize_E2_;

    double spirit_reg_lamda_;
    double spirit_calib_over_determine_ratio_;

    bool spirit_solve_symmetric_;

    unsigned long long spirit_iter_max_;
    double spirit_iter_thres_;
    bool spirit_print_iter_;

    // L1 SPIRiT
    bool spirit_perform_linear_;
    bool spirit_perform_nonlinear_;

    double spirit_parallel_imaging_lamda_;
    double spirit_image_reg_lamda_;
    double spirit_data_fidelity_lamda_;

    unsigned long long spirit_ncg_iter_max_;
    double spirit_ncg_iter_thres_;
    bool spirit_ncg_print_iter_;
    double spirit_ncg_scale_factor_;

    bool spirit_use_coil_sen_map_;
    bool spirit_use_moco_enhancement_;
    bool spirit_recon_moco_images_;

    bool spirit_2D_scale_per_chunk_;
    bool spirit_3D_scale_per_chunk_;

    double spirit_RO_enhancement_ratio_;
    double spirit_E1_enhancement_ratio_;
    double spirit_E2_enhancement_ratio_;
    double spirit_temporal_enhancement_ratio_;

    // L1 soft sense

    // -------------------------------
    // job split
    // -------------------------------
    bool job_split_by_S_;
    unsigned long long job_num_of_N_;
    unsigned long long job_max_Megabytes_;
    unsigned long long job_overlap_;
    // whether to perform computation on the control node
    bool job_perform_on_control_node_;

    // -------------------------------
    // partial fourier handling
    // -------------------------------
    // partial fourier handling algorithms
    ISMRMRDPFALGO partialFourier_algo_;

    // homodyne filter
    // number of iterations
    unsigned long long partialFourier_homodyne_iters_;
    // threshold to stop the iteration
    double partialFourier_homodyne_thres_;
    // density compensation for homodyne filter results
    bool partialFourier_homodyne_densityComp_;

    // POCS
    // number of iterations
    unsigned long long partialFourier_POCS_iters_;
    // threshold to stop the iteration
    double partialFourier_POCS_thres_;
    // transit band width
    unsigned long long partialFourier_POCS_transitBand_;
    // transit band width for E2
    unsigned long long partialFourier_POCS_transitBand_E2_;

    // Feng Huang method
    // kernel size
    unsigned long long partialFourier_FengHuang_kSize_RO_;
    unsigned long long partialFourier_FengHuang_kSize_E1_;
    unsigned long long partialFourier_FengHuang_kSize_E2_;
    // threshold for kernel estimation
    double partialFourier_FengHuang_thresReg_;
    // same kernel for all N
    bool partialFourier_FengHuang_sameKernel_allN_;
    // transit band width
    unsigned long long partialFourier_FengHuang_transitBand_;
    // transit band width for E2
    unsigned long long partialFourier_FengHuang_transitBand_E2_;

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

        upstream_coil_compression_ = false;
        upstream_coil_compression_thres_ = 1e-3;
        upstream_coil_compression_num_modesKept_ = -1;

        downstream_coil_compression_ = true;
        coil_compression_thres_ = 1e-3;
        coil_compression_num_modesKept_ = -1;

        coil_map_algorithm_ = ISMRMRD_SOUHEIL;
        csm_kSize_ = 7;
        csm_powermethod_num_ = 3;
        csm_iter_num_ = 5;
        csm_iter_thres_ = 1e-3;

        recon_algorithm_ = ISMRMRD_GRAPPA;
        recon_auto_parameters_ = true;

        grappa_kSize_RO_ = 5;
        grappa_kSize_E1_ = 4;
        grappa_kSize_E2_ = 4;
        grappa_reg_lamda_ = 0.0005;

        spirit_kSize_RO_ = 7;
        spirit_kSize_E1_ = 7;
        spirit_kSize_E2_ = 7;

        spirit_reg_lamda_ = 0.005;
        spirit_calib_over_determine_ratio_ = 0;

        spirit_solve_symmetric_ = false;

        spirit_iter_max_ = 70;
        spirit_iter_thres_ = 1e-5;
        spirit_print_iter_ = false;

        spirit_perform_linear_ = true;
        spirit_perform_nonlinear_ = true;

        spirit_parallel_imaging_lamda_ = 1.0;
        spirit_image_reg_lamda_ = 1e-3;
        spirit_data_fidelity_lamda_ = 0;

        spirit_ncg_iter_max_ = 10;
        spirit_ncg_iter_thres_ = 1e-3;
        spirit_ncg_print_iter_ = false;
        spirit_ncg_scale_factor_ = 1.0;

        spirit_use_coil_sen_map_ = true;
        spirit_use_moco_enhancement_ = false;
        spirit_recon_moco_images_ = false;

        spirit_RO_enhancement_ratio_ = 1;
        spirit_E1_enhancement_ratio_ = 1;
        spirit_E2_enhancement_ratio_ = 1;
        spirit_temporal_enhancement_ratio_ = 1;

        spirit_2D_scale_per_chunk_ = false;
        spirit_3D_scale_per_chunk_ = true;

        job_split_by_S_ = false;
        job_num_of_N_ = 0;
        job_max_Megabytes_ = 20*1024;
        job_overlap_ = 2;
        job_perform_on_control_node_ = true;

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

    gtPlusReconWorkOrder();
    virtual ~gtPlusReconWorkOrder();

    // reset the status of work order
    // all computed calibration/coil sensitivity results
    // are deleted
    virtual bool reset();

    // check and modify inconsistency in the work order
    virtual bool enforceConsistency(ISMRMRDDIM& /*lastDim*/);

    typedef std::pair<ISMRMRDDIM, unsigned long long> DimensionRecordType;

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

    // gfactor
    hoNDArray<T> gfactor_;

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
    worder.CalibMode_ = CalibMode_;
    worder.InterleaveDim_ = InterleaveDim_;

    worder.acceFactorE1_ = acceFactorE1_;
    worder.acceFactorE2_ = acceFactorE2_;

    worder.kSpaceCenterRO_ = kSpaceCenterRO_;
    worder.kSpaceCenterEncode1_ = kSpaceCenterEncode1_;
    worder.kSpaceCenterEncode2_ = kSpaceCenterEncode2_;

    worder.kSpaceMaxRO_ = kSpaceMaxRO_;
    worder.kSpaceMaxEncode1_ = kSpaceMaxEncode1_;
    worder.kSpaceMaxEncode2_ = kSpaceMaxEncode2_;

    worder.workFlow_BufferKernel_ = workFlow_BufferKernel_;
    worder.workFlow_use_BufferedKernel_ = workFlow_use_BufferedKernel_;
    worder.num_channels_res_ = num_channels_res_;

    worder.upstream_coil_compression_ = upstream_coil_compression_;
    worder.upstream_coil_compression_thres_ = upstream_coil_compression_thres_;
    worder.upstream_coil_compression_num_modesKept_ = upstream_coil_compression_num_modesKept_;

    worder.downstream_coil_compression_ = downstream_coil_compression_;
    worder.coil_compression_thres_ = coil_compression_thres_;
    worder.coil_compression_num_modesKept_ = coil_compression_num_modesKept_;

    worder.coil_map_algorithm_ = coil_map_algorithm_;
    worder.csm_kSize_ = csm_kSize_;
    worder.csm_powermethod_num_ = csm_powermethod_num_;
    worder.csm_iter_num_ = csm_iter_num_;
    worder.csm_iter_thres_ = csm_iter_thres_;

    worder.start_RO_ = start_RO_;
    worder.end_RO_ = end_RO_;

    worder.start_E1_ = start_E1_;
    worder.end_E1_ = end_E1_;

    worder.start_E2_ = start_E2_;
    worder.end_E2_ = end_E2_;

    worder.recon_algorithm_ = recon_algorithm_;
    worder.recon_auto_parameters_ = recon_auto_parameters_;

    worder.grappa_kSize_RO_ = grappa_kSize_RO_;
    worder.grappa_kSize_RO_ = grappa_kSize_RO_;
    worder.grappa_kSize_E1_ = grappa_kSize_E1_;
    worder.grappa_kSize_E2_ = grappa_kSize_E2_;
    worder.grappa_reg_lamda_ = grappa_reg_lamda_;

    worder.spirit_kSize_RO_ = spirit_kSize_RO_;
    worder.spirit_kSize_E1_ = spirit_kSize_E1_;
    worder.spirit_kSize_E2_ = spirit_kSize_E2_;
    worder.spirit_reg_lamda_ = spirit_reg_lamda_;
    worder.spirit_calib_over_determine_ratio_ = spirit_calib_over_determine_ratio_;
    worder.spirit_solve_symmetric_ = spirit_solve_symmetric_;
    worder.spirit_iter_max_ = spirit_iter_max_;
    worder.spirit_iter_thres_ = spirit_iter_thres_;
    worder.spirit_print_iter_ = spirit_print_iter_;

    worder.spirit_perform_linear_ = spirit_perform_linear_;
    worder.spirit_perform_nonlinear_ = spirit_perform_nonlinear_;
    worder.spirit_parallel_imaging_lamda_ = spirit_parallel_imaging_lamda_;
    worder.spirit_image_reg_lamda_ = spirit_image_reg_lamda_;
    worder.spirit_data_fidelity_lamda_ = spirit_data_fidelity_lamda_;
    worder.spirit_ncg_iter_max_ = spirit_ncg_iter_max_;
    worder.spirit_ncg_iter_thres_ = spirit_ncg_iter_thres_;
    worder.spirit_ncg_scale_factor_ = spirit_ncg_scale_factor_;
    worder.spirit_ncg_print_iter_ = spirit_ncg_print_iter_;
    worder.spirit_use_coil_sen_map_ = spirit_use_coil_sen_map_;
    worder.spirit_use_moco_enhancement_ = spirit_use_moco_enhancement_;
    worder.spirit_recon_moco_images_ = spirit_recon_moco_images_;
    worder.spirit_RO_enhancement_ratio_ = spirit_RO_enhancement_ratio_;
    worder.spirit_E1_enhancement_ratio_ = spirit_E1_enhancement_ratio_;
    worder.spirit_E2_enhancement_ratio_ = spirit_E2_enhancement_ratio_;
    worder.spirit_temporal_enhancement_ratio_ = spirit_temporal_enhancement_ratio_;
    worder.spirit_2D_scale_per_chunk_ = spirit_2D_scale_per_chunk_;
    worder.spirit_3D_scale_per_chunk_ = spirit_3D_scale_per_chunk_;

    worder.job_split_by_S_ = job_split_by_S_;
    worder.job_num_of_N_ = job_num_of_N_;
    worder.job_max_Megabytes_ = job_max_Megabytes_;
    worder.job_overlap_ = job_overlap_;
    worder.job_perform_on_control_node_ = job_perform_on_control_node_;

    worder.partialFourier_algo_ = partialFourier_algo_;

    worder.partialFourier_homodyne_iters_ = partialFourier_homodyne_iters_;
    worder.partialFourier_homodyne_thres_ = partialFourier_homodyne_thres_;
    worder.partialFourier_homodyne_densityComp_ = partialFourier_homodyne_densityComp_;

    worder.partialFourier_POCS_iters_ = partialFourier_POCS_iters_;
    worder.partialFourier_POCS_thres_ = partialFourier_POCS_thres_;
    worder.partialFourier_POCS_transitBand_ = partialFourier_POCS_transitBand_;
    worder.partialFourier_POCS_transitBand_E2_ = partialFourier_POCS_transitBand_E2_;

    worder.partialFourier_FengHuang_kSize_RO_ = partialFourier_FengHuang_kSize_RO_;
    worder.partialFourier_FengHuang_kSize_E1_ = partialFourier_FengHuang_kSize_E1_;
    worder.partialFourier_FengHuang_kSize_E2_ = partialFourier_FengHuang_kSize_E2_;
    worder.partialFourier_FengHuang_thresReg_ = partialFourier_FengHuang_thresReg_;
    worder.partialFourier_FengHuang_sameKernel_allN_ = partialFourier_FengHuang_sameKernel_allN_;
    worder.partialFourier_FengHuang_transitBand_ = partialFourier_FengHuang_transitBand_;
    worder.partialFourier_FengHuang_transitBand_E2_ = partialFourier_FengHuang_transitBand_E2_;
}

template <typename T> 
void gtPlusReconWorkOrder<T>::duplicate(gtPlusReconWorkOrder<T>& worder) const
{
    this->duplicatePara(worder);

    worder.dataDimStartingIndexes_ = dataDimStartingIndexes_;

    worder.filterRO_ = filterRO_;
    worder.filterE1_ = filterE1_;
    worder.filterE2_ = filterE2_;
    worder.filterROE1_ = filterROE1_;
    worder.filterROE1E2_ = filterROE1E2_;

    worder.filterRO_ref_ = filterRO_ref_;
    worder.filterE1_ref_ = filterE1_ref_;
    worder.filterE2_ref_ = filterE2_ref_;
    worder.filterROE1_ref_ = filterROE1_ref_;
    worder.filterROE1E2_ref_ = filterROE1E2_ref_;

    worder.filterRO_partialfourier_ = filterRO_partialfourier_;
    worder.filterE1_partialfourier_ = filterE1_partialfourier_;
    worder.filterE2_partialfourier_ = filterE2_partialfourier_;
    worder.filterROE1_partialfourier_ = filterROE1_partialfourier_;
    worder.filterROE1E2_partialfourier_ = filterROE1E2_partialfourier_;

    worder.CloudComputing_ = CloudComputing_;
    worder.CloudSize_ = CloudSize_;
    worder.gt_cloud_ = gt_cloud_;
}

template <typename T> 
void gtPlusReconWorkOrder<T>::copyFromPara(const gtPlusReconWorkOrderPara& worder)
{
    CalibMode_ = worder.CalibMode_;
    InterleaveDim_ = worder.InterleaveDim_;

    acceFactorE1_ = worder.acceFactorE1_;
    acceFactorE2_ = worder.acceFactorE2_;

    kSpaceCenterRO_ = worder.kSpaceCenterRO_;
    kSpaceCenterEncode1_ = worder.kSpaceCenterEncode1_;
    kSpaceCenterEncode2_ = worder.kSpaceCenterEncode2_;

    kSpaceMaxRO_ = worder.kSpaceMaxRO_;
    kSpaceMaxEncode1_ = worder.kSpaceMaxEncode1_;
    kSpaceMaxEncode2_ = worder.kSpaceMaxEncode2_;

    workFlow_BufferKernel_ = worder.workFlow_BufferKernel_;
    workFlow_use_BufferedKernel_ = worder.workFlow_use_BufferedKernel_;
    num_channels_res_ = worder.num_channels_res_;

    upstream_coil_compression_ = worder.upstream_coil_compression_;
    upstream_coil_compression_thres_ = worder.upstream_coil_compression_thres_;
    upstream_coil_compression_num_modesKept_ = worder.upstream_coil_compression_num_modesKept_;

    downstream_coil_compression_ = worder.downstream_coil_compression_;
    coil_compression_thres_ = worder.coil_compression_thres_;
    coil_compression_num_modesKept_ = worder.coil_compression_num_modesKept_;

    coil_map_algorithm_ = worder.coil_map_algorithm_;
    csm_kSize_ = worder.csm_kSize_;
    csm_powermethod_num_ = worder.csm_powermethod_num_;
    csm_iter_num_ = worder.csm_iter_num_;
    csm_iter_thres_ = worder.csm_iter_thres_;

    start_RO_ = worder.start_RO_;
    end_RO_ = worder.end_RO_;

    start_E1_ = worder.start_E1_;
    end_E1_ = worder.end_E1_;

    start_E2_ = worder.start_E2_;
    end_E2_ = worder.end_E2_;

    recon_algorithm_ = worder.recon_algorithm_;
    recon_auto_parameters_ = worder.recon_auto_parameters_;

    grappa_kSize_RO_ = worder.grappa_kSize_RO_;
    grappa_kSize_RO_ = worder.grappa_kSize_RO_;
    grappa_kSize_E1_ = worder.grappa_kSize_E1_;
    grappa_kSize_E2_ = worder.grappa_kSize_E2_;
    grappa_reg_lamda_ = worder.grappa_reg_lamda_;

    spirit_kSize_RO_ = worder.spirit_kSize_RO_;
    spirit_kSize_E1_ = worder.spirit_kSize_E1_;
    spirit_kSize_E2_ = worder.spirit_kSize_E2_;
    spirit_reg_lamda_ = worder.spirit_reg_lamda_;
    spirit_calib_over_determine_ratio_ = worder.spirit_calib_over_determine_ratio_;
    spirit_solve_symmetric_ = worder.spirit_solve_symmetric_;
    spirit_iter_max_ = worder.spirit_iter_max_;
    spirit_iter_thres_ = worder.spirit_iter_thres_;
    spirit_print_iter_ = worder.spirit_print_iter_;

    spirit_perform_linear_ = worder.spirit_perform_linear_;
    spirit_perform_nonlinear_ = worder.spirit_perform_nonlinear_;
    spirit_parallel_imaging_lamda_ = worder.spirit_parallel_imaging_lamda_;
    spirit_image_reg_lamda_ = worder.spirit_image_reg_lamda_;
    spirit_data_fidelity_lamda_ = worder.spirit_data_fidelity_lamda_;
    spirit_ncg_iter_max_ = worder.spirit_ncg_iter_max_;
    spirit_ncg_iter_thres_ = worder.spirit_ncg_iter_thres_;
    spirit_ncg_scale_factor_ = worder.spirit_ncg_scale_factor_;
    spirit_ncg_print_iter_ = worder.spirit_ncg_print_iter_;
    spirit_use_coil_sen_map_ = worder.spirit_use_coil_sen_map_;
    spirit_use_moco_enhancement_ = worder.spirit_use_moco_enhancement_;
    spirit_recon_moco_images_ = worder.spirit_recon_moco_images_;
    spirit_RO_enhancement_ratio_ = worder.spirit_RO_enhancement_ratio_;
    spirit_E1_enhancement_ratio_ = worder.spirit_E1_enhancement_ratio_;
    spirit_E2_enhancement_ratio_ = worder.spirit_E2_enhancement_ratio_;
    spirit_temporal_enhancement_ratio_ = worder.spirit_temporal_enhancement_ratio_;
    spirit_2D_scale_per_chunk_ = worder.spirit_2D_scale_per_chunk_;
    spirit_3D_scale_per_chunk_ = worder.spirit_3D_scale_per_chunk_;

    job_split_by_S_ = worder.job_split_by_S_;
    job_num_of_N_ = worder.job_num_of_N_;
    job_max_Megabytes_ = worder.job_max_Megabytes_;
    job_overlap_ = worder.job_overlap_;
    job_perform_on_control_node_ = worder.job_perform_on_control_node_;

    partialFourier_algo_ = worder.partialFourier_algo_;

    partialFourier_homodyne_iters_ = worder.partialFourier_homodyne_iters_;
    partialFourier_homodyne_thres_ = worder.partialFourier_homodyne_thres_;
    partialFourier_homodyne_densityComp_ = worder.partialFourier_homodyne_densityComp_;

    partialFourier_POCS_iters_ = worder.partialFourier_POCS_iters_;
    partialFourier_POCS_thres_ = worder.partialFourier_POCS_thres_;
    partialFourier_POCS_transitBand_ = worder.partialFourier_POCS_transitBand_;
    partialFourier_POCS_transitBand_E2_ = worder.partialFourier_POCS_transitBand_E2_;

    partialFourier_FengHuang_kSize_RO_ = worder.partialFourier_FengHuang_kSize_RO_;
    partialFourier_FengHuang_kSize_E1_ = worder.partialFourier_FengHuang_kSize_E1_;
    partialFourier_FengHuang_kSize_E2_ = worder.partialFourier_FengHuang_kSize_E2_;
    partialFourier_FengHuang_thresReg_ = worder.partialFourier_FengHuang_thresReg_;
    partialFourier_FengHuang_sameKernel_allN_ = worder.partialFourier_FengHuang_sameKernel_allN_;
    partialFourier_FengHuang_transitBand_ = worder.partialFourier_FengHuang_transitBand_;
    partialFourier_FengHuang_transitBand_E2_ = worder.partialFourier_FengHuang_transitBand_E2_;
}

template <typename T> 
void gtPlusReconWorkOrder<T>::printInfo(std::ostream& os) const
{
    using namespace std;
    GADGET_OSTREAM_PRINT(os, CalibMode_);
    GADGET_OSTREAM_PRINT(os, InterleaveDim_);
    GADGET_OSTREAM_PRINT(os, acceFactorE1_);
    GADGET_OSTREAM_PRINT(os, acceFactorE2_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, kSpaceCenterRO_);
    GADGET_OSTREAM_PRINT(os, kSpaceCenterEncode1_);
    GADGET_OSTREAM_PRINT(os, kSpaceCenterEncode2_);
    GADGET_OSTREAM_PRINT(os, kSpaceMaxRO_);
    GADGET_OSTREAM_PRINT(os, kSpaceMaxEncode1_);
    GADGET_OSTREAM_PRINT(os, kSpaceMaxEncode2_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, workFlow_BufferKernel_);
    GADGET_OSTREAM_PRINT(os, workFlow_use_BufferedKernel_);
    GADGET_OSTREAM_PRINT(os, num_channels_res_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, upstream_coil_compression_);
    GADGET_OSTREAM_PRINT(os, upstream_coil_compression_thres_);
    GADGET_OSTREAM_PRINT(os, upstream_coil_compression_num_modesKept_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, downstream_coil_compression_);
    GADGET_OSTREAM_PRINT(os, coil_compression_thres_);
    GADGET_OSTREAM_PRINT(os, coil_compression_num_modesKept_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, coil_map_algorithm_);
    GADGET_OSTREAM_PRINT(os, csm_kSize_);
    GADGET_OSTREAM_PRINT(os, csm_powermethod_num_);
    GADGET_OSTREAM_PRINT(os, csm_iter_num_);
    GADGET_OSTREAM_PRINT(os, csm_iter_thres_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, start_RO_);
    GADGET_OSTREAM_PRINT(os, end_RO_);
    GADGET_OSTREAM_PRINT(os, start_E1_);
    GADGET_OSTREAM_PRINT(os, end_E1_);
    GADGET_OSTREAM_PRINT(os, start_E2_);
    GADGET_OSTREAM_PRINT(os, end_E2_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, recon_algorithm_);
    GADGET_OSTREAM_PRINT(os, recon_auto_parameters_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, grappa_kSize_RO_);
    GADGET_OSTREAM_PRINT(os, grappa_kSize_E1_);
    GADGET_OSTREAM_PRINT(os, grappa_kSize_E2_);
    GADGET_OSTREAM_PRINT(os, grappa_reg_lamda_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, spirit_kSize_RO_);
    GADGET_OSTREAM_PRINT(os, spirit_kSize_E1_);
    GADGET_OSTREAM_PRINT(os, spirit_kSize_E2_);
    GADGET_OSTREAM_PRINT(os, spirit_reg_lamda_);
    GADGET_OSTREAM_PRINT(os, spirit_calib_over_determine_ratio_);
    GADGET_OSTREAM_PRINT(os, spirit_solve_symmetric_);
    GADGET_OSTREAM_PRINT(os, spirit_iter_max_);
    GADGET_OSTREAM_PRINT(os, spirit_iter_thres_);
    GADGET_OSTREAM_PRINT(os, spirit_print_iter_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, spirit_perform_linear_);
    GADGET_OSTREAM_PRINT(os, spirit_perform_nonlinear_);
    GADGET_OSTREAM_PRINT(os, spirit_parallel_imaging_lamda_);
    GADGET_OSTREAM_PRINT(os, spirit_image_reg_lamda_);
    GADGET_OSTREAM_PRINT(os, spirit_data_fidelity_lamda_);
    GADGET_OSTREAM_PRINT(os, spirit_ncg_iter_max_);
    GADGET_OSTREAM_PRINT(os, spirit_ncg_iter_thres_);
    GADGET_OSTREAM_PRINT(os, spirit_ncg_scale_factor_);
    GADGET_OSTREAM_PRINT(os, spirit_ncg_print_iter_);
    GADGET_OSTREAM_PRINT(os, spirit_use_coil_sen_map_);
    GADGET_OSTREAM_PRINT(os, spirit_use_moco_enhancement_);
    GADGET_OSTREAM_PRINT(os, spirit_recon_moco_images_);
    GADGET_OSTREAM_PRINT(os, spirit_RO_enhancement_ratio_);
    GADGET_OSTREAM_PRINT(os, spirit_E1_enhancement_ratio_);
    GADGET_OSTREAM_PRINT(os, spirit_E2_enhancement_ratio_);
    GADGET_OSTREAM_PRINT(os, spirit_temporal_enhancement_ratio_);
    GADGET_OSTREAM_PRINT(os, spirit_2D_scale_per_chunk_);
    GADGET_OSTREAM_PRINT(os, spirit_3D_scale_per_chunk_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, job_split_by_S_);
    GADGET_OSTREAM_PRINT(os, job_num_of_N_);
    GADGET_OSTREAM_PRINT(os, job_max_Megabytes_);
    GADGET_OSTREAM_PRINT(os, job_overlap_);
    GADGET_OSTREAM_PRINT(os, job_perform_on_control_node_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, partialFourier_algo_);
    GADGET_OSTREAM_PRINT(os, partialFourier_homodyne_iters_);
    GADGET_OSTREAM_PRINT(os, partialFourier_homodyne_thres_);
    GADGET_OSTREAM_PRINT(os, partialFourier_homodyne_densityComp_);
    GADGET_OSTREAM_PRINT(os, partialFourier_POCS_iters_);
    GADGET_OSTREAM_PRINT(os, partialFourier_POCS_thres_);
    GADGET_OSTREAM_PRINT(os, partialFourier_POCS_transitBand_);
    GADGET_OSTREAM_PRINT(os, partialFourier_POCS_transitBand_E2_);
    GADGET_OSTREAM_PRINT(os, partialFourier_FengHuang_kSize_RO_);
    GADGET_OSTREAM_PRINT(os, partialFourier_FengHuang_kSize_E1_);
    GADGET_OSTREAM_PRINT(os, partialFourier_FengHuang_kSize_E2_);
    GADGET_OSTREAM_PRINT(os, partialFourier_FengHuang_thresReg_);
    GADGET_OSTREAM_PRINT(os, partialFourier_FengHuang_sameKernel_allN_);
    GADGET_OSTREAM_PRINT(os, partialFourier_FengHuang_transitBand_);
    GADGET_OSTREAM_PRINT(os, partialFourier_FengHuang_transitBand_E2_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, CloudComputing_);
    GADGET_OSTREAM_PRINT(os, CloudSize_);
    for ( unsigned int nn=0; nn<gt_cloud_.size(); nn++ )
    {
        GADGET_OSTREAM_PRINT(os, gt_cloud_[nn]);
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
