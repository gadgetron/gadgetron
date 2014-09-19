/** \file   gtPlusISMRMRDReconWorkOrder3DT.h
    \brief  Define the GtPlus reconstruction workorder and parameters for 3DT reconstruction
    \author Hui Xue
*/

#pragma once

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorkOrder3DT : public gtPlusReconWorkOrder<T>
{
public:

    typedef gtPlusReconWorkOrder<T> BaseClass;

    gtPlusReconWorkOrder3DT();
    virtual ~gtPlusReconWorkOrder3DT();

    // reset the status of work order
    // all computed calibration/coil sensitivity results are deleted
    virtual bool reset();

    // check and modify inconsistency in the work order
    virtual bool enforceConsistency(ISMRMRDDIM& lastDim);

    virtual void duplicate(gtPlusReconWorkOrder3DT<T>& worder);

    virtual void printInfo(std::ostream& os) const;
    virtual void print(std::ostream& os) const;

    // kspace_: [RO E1 E2 CHA N], for 3D recon, N can be 1
    // ref_: [RO E1 E2 CHA M], M can equal to N or 1 or others
    // fullkspace_: [RO E1 E2 CHA N]
    // complexIm_: [RO E1 E2 N], after coil combination
    // coilMap_: [RO E1 E2 CHA 1 or N]
    // gfactor_: [RO E1 E2 CHA 1 or N]

    // the fifth dimension can be the temporal dimension or others

    // default behavior
    // a) the coil compression coefficients are computed once for all N

    // embedded mode
    // a) perform recon and estimate kernel/coil sensitivity for every 3D kspace [RO E1 E2 CHA]
    // b) coil combination uses different coil maps for every N
    // c) if the kspace recon is performed, the coil combination map is reestimated on the fullkspace for every 3D images
    // d) the ref lines are filled back to fullkspace_

    // separate mode
    // a) perform recon and estimate kernel/coil sensitivity for every 3D kspace [RO E1 E2 CHA] if M==N
    // b) if M==1, the kernel is only estimated once for all N
    // c) coil combination uses different coil maps for every N
    // d) if the kspace recon is performed, the coil combination map is reestimated on the fullkspace for every 3D images

    // interleave
    // a) the average-all ref is used
    // b) kernel/coil sensitivity is estimated once for all N

    using BaseClass::data_;
    using BaseClass::ref_;
    using BaseClass::ref_recon_;
    using BaseClass::ref_coil_map_;
    using BaseClass::CalibMode_;
    using BaseClass::InterleaveDim_;
    using BaseClass::acceFactorE1_;
    using BaseClass::acceFactorE2_;
    using BaseClass::num_channels_res_;

    using BaseClass::coilMap_;
    using BaseClass::fullkspace_;
    using BaseClass::complexIm_;
    using BaseClass::recon_time_stamp_;
    using BaseClass::recon_physio_time_stamp_;

    using BaseClass::fullkspace_second_;
    using BaseClass::complexIm_second_;
    using BaseClass::recon_time_stamp_second_;
    using BaseClass::recon_physio_time_stamp_second_;

    using BaseClass::gfactor_;
    using BaseClass::wrap_around_map_;

    using BaseClass::upstream_coil_compression_;
    using BaseClass::upstream_coil_compression_thres_;
    using BaseClass::upstream_coil_compression_num_modesKept_;

    using BaseClass::downstream_coil_compression_;
    using BaseClass::coil_compression_thres_;
    using BaseClass::coil_compression_num_modesKept_;
    using BaseClass::csm_kSize_;
    using BaseClass::csm_powermethod_num_;
    using BaseClass::csm_true_3D_;
    using BaseClass::csm_iter_num_;
    using BaseClass::csm_iter_thres_;
    using BaseClass::csm_use_gpu_;
    using BaseClass::start_RO_;
    using BaseClass::end_RO_;
    using BaseClass::start_E1_;
    using BaseClass::end_E1_;
    using BaseClass::start_E2_;
    using BaseClass::end_E2_;

    using BaseClass::filterRO_;
    using BaseClass::filterE1_;
    using BaseClass::filterE2_;
    using BaseClass::filterROE1_;
    using BaseClass::filterROE1E2_;

    using BaseClass::filterRO_ref_;
    using BaseClass::filterE1_ref_;
    using BaseClass::filterE2_ref_;
    using BaseClass::filterROE1_ref_;
    using BaseClass::filterROE1E2_ref_;

    using BaseClass::filterRO_partialfourier_;
    using BaseClass::filterE1_partialfourier_;
    using BaseClass::filterE2_partialfourier_;
    using BaseClass::filterROE1_partialfourier_;
    using BaseClass::filterROE1E2_partialfourier_;

    using BaseClass::recon_algorithm_;

    using BaseClass::grappa_kSize_RO_;
    using BaseClass::grappa_kSize_E1_;
    using BaseClass::grappa_kSize_E2_;
    using BaseClass::grappa_reg_lamda_;
    using BaseClass::grappa_calib_over_determine_ratio_;
    using BaseClass::grappa_use_gpu_;

    using BaseClass::spirit_kSize_RO_;
    using BaseClass::spirit_kSize_E1_;
    using BaseClass::spirit_kSize_E2_;
    using BaseClass::spirit_oSize_RO_;
    using BaseClass::spirit_oSize_E1_;
    using BaseClass::spirit_oSize_E2_;
    using BaseClass::spirit_reg_lamda_;
    using BaseClass::spirit_use_gpu_;
    using BaseClass::spirit_iter_max_;
    using BaseClass::spirit_iter_thres_;
    using BaseClass::spirit_print_iter_;

    using BaseClass::spirit_perform_linear_;
    using BaseClass::spirit_perform_nonlinear_;
    using BaseClass::spirit_parallel_imaging_lamda_;
    using BaseClass::spirit_image_reg_lamda_;
    using BaseClass::spirit_data_fidelity_lamda_;
    using BaseClass::spirit_ncg_iter_max_;
    using BaseClass::spirit_ncg_iter_thres_;
    using BaseClass::spirit_ncg_scale_factor_;
    using BaseClass::spirit_ncg_print_iter_;
    using BaseClass::spirit_slep_iter_max_;
    using BaseClass::spirit_slep_iter_thres_;
    using BaseClass::spirit_slep_print_iter_;
    using BaseClass::spirit_slep_keep_third_dimension_coeff_;
    using BaseClass::spirit_slep_scale_factor_;
    using BaseClass::spirit_use_coil_sen_map_;
    using BaseClass::spirit_use_moco_enhancement_;
    using BaseClass::spirit_recon_moco_images_;
    using BaseClass::spirit_RO_enhancement_ratio_;
    using BaseClass::spirit_E1_enhancement_ratio_;
    using BaseClass::spirit_E2_enhancement_ratio_;
    using BaseClass::spirit_temporal_enhancement_ratio_;

    using BaseClass::job_split_by_S_;
    using BaseClass::job_num_of_N_;
    using BaseClass::job_max_Megabytes_;
    using BaseClass::job_overlap_;

    using BaseClass::partialFourier_algo_;
    using BaseClass::partialFourier_homodyne_iters_;
    using BaseClass::partialFourier_homodyne_thres_;
    using BaseClass::partialFourier_homodyne_densityComp_;
    using BaseClass::partialFourier_POCS_iters_;
    using BaseClass::partialFourier_POCS_thres_;
    using BaseClass::partialFourier_POCS_transitBand_;
    using BaseClass::partialFourier_FengHuang_kSize_RO_;
    using BaseClass::partialFourier_FengHuang_kSize_E1_;
    using BaseClass::partialFourier_FengHuang_kSize_E2_;
    using BaseClass::partialFourier_FengHuang_thresReg_;
    using BaseClass::partialFourier_FengHuang_sameKernel_allN_;
    using BaseClass::partialFourier_FengHuang_transitBand_;
    using BaseClass::partialFourier_FengHuang_transitBand_E2_;

    using BaseClass::CloudComputing_;
    using BaseClass::CloudSize_;
    using BaseClass::gt_cloud_;

    using BaseClass::kernel_; // [RO E1 E2 srcCHA dstCHA dstRO dstE1 dstE2]
    using BaseClass::kernelIm_; // [RO E1 E2 srcCHA dstCHA]
    using BaseClass::unmixingCoeffIm_; // [RO E1 E2 srcCHA 1 or N]
    using BaseClass::coilCompressionCoef_;

    // parameters to change the default behavior

    // if true, the actual full kspace is computed, not only the coil combined complex images
    bool recon_kspace_needed_;

    // if true, no coil compression will be performed
    bool coil_compression_;
    // if true, the same coil compression coefficient is computed for all N
    bool same_coil_compression_coeff_allN_;

    // no acceleration
    // if true, the average of all M ref will be used
    // the coil sensitivity will be only estimated once for all N
    bool no_acceleration_averageall_ref_;
    // if true, the same coil combination coefficients will be used for all N
    bool no_acceleration_same_combinationcoeff_allN_;
    // if no_acceleration_same_combinationcoeff_allN_==true, select the N for coil combination coefficient estimation
    size_t no_acceleration_whichN_combinationcoeff_;

    // embedded mode
    // if true, the average of all M ref will be used
    // the kernel/sensitivity will be only estimated once for all N
    bool embedded_averageall_ref_;
    // if true, the coil map will be estimated from the fullkspace_
    bool embedded_fullres_coilmap_;
    // if true, the same coil combination coefficients will be used for all N
    bool embedded_same_combinationcoeff_allN_;
    // if embedded_same_combinationcoeff_allN_==true, select the N for coil combination coefficient estimation
    // if -1, the average-all N is used for coil combination
    int embedded_whichN_combinationcoeff_;
    // if true, the ref lines will be filled back to fullkspace
    bool embedded_ref_fillback_;

    // separate mode
    // if true, the average of all M ref will be used
    // the kernel/sensitivity will be only estimated once for all N
    bool separate_averageall_ref_;
    // if true, the coil map will be estimated from the fullkspace_
    bool separate_fullres_coilmap_;
    // if true, the same coil combination coefficients will be used for all N
    bool separate_same_combinationcoeff_allN_;
    // if separate_same_combinationcoeff_allN_==true, select the 3D kspace used for coil combination coefficient estimation
    // if -1, the average-all N is used for coil combination
    int separate_whichN_combinationcoeff_;

    // interleaved mode
};

template <typename T> 
gtPlusReconWorkOrder3DT<T>::gtPlusReconWorkOrder3DT() : BaseClass()
{
    recon_kspace_needed_ = false;
    coil_compression_ = true;
    same_coil_compression_coeff_allN_ = false;

    no_acceleration_averageall_ref_ = false;
    no_acceleration_same_combinationcoeff_allN_ = false;

    embedded_averageall_ref_ = false;
    embedded_fullres_coilmap_ = true;
    embedded_same_combinationcoeff_allN_ = false;
    embedded_whichN_combinationcoeff_ = false;
    embedded_ref_fillback_ = true;

    separate_averageall_ref_ = false;
    separate_fullres_coilmap_ = true;
    separate_same_combinationcoeff_allN_ = false;
    separate_whichN_combinationcoeff_ = false;
}

template <typename T> 
gtPlusReconWorkOrder3DT<T>::~gtPlusReconWorkOrder3DT()
{
}

template <typename T> 
bool gtPlusReconWorkOrder3DT<T>::reset()
{
    try
    {
        kernel_->clear();
        kernelIm_->clear();
        unmixingCoeffIm_->clear();
        coilCompressionCoef_->clear();
        coilMap_->clear();

        fullkspace_.clear();
        complexIm_.clear();
        recon_time_stamp_.clear();
        recon_physio_time_stamp_.clear();

        fullkspace_second_.clear();
        complexIm_second_.clear();
        recon_time_stamp_second_.clear();
        recon_physio_time_stamp_second_.clear();

        gfactor_.clear();
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusReconWorkOrder3DT<T>::reset() ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorkOrder3DT<T>::enforceConsistency(ISMRMRDDIM& lastDim)
{
    if ( lastDim == DIM_Slice )
    {
        no_acceleration_averageall_ref_ = false;
        no_acceleration_same_combinationcoeff_allN_ = false;

        embedded_averageall_ref_ = false;
        embedded_same_combinationcoeff_allN_ = false;

        separate_averageall_ref_ = false;
        separate_same_combinationcoeff_allN_ = false;
    }

    return true;
}

template <typename T> 
void gtPlusReconWorkOrder3DT<T>::duplicate(gtPlusReconWorkOrder3DT<T>& worder)
{
    BaseClass::duplicate(worder);

    worder.recon_kspace_needed_ = recon_kspace_needed_;
    worder.coil_compression_ = coil_compression_;
    worder.same_coil_compression_coeff_allN_ = same_coil_compression_coeff_allN_;

    worder.no_acceleration_averageall_ref_ = no_acceleration_averageall_ref_;
    worder.no_acceleration_same_combinationcoeff_allN_ = no_acceleration_same_combinationcoeff_allN_;
    worder.no_acceleration_whichN_combinationcoeff_ = no_acceleration_whichN_combinationcoeff_;

    worder.embedded_averageall_ref_ = embedded_averageall_ref_;
    worder.embedded_fullres_coilmap_ = embedded_fullres_coilmap_;
    worder.embedded_same_combinationcoeff_allN_ = embedded_same_combinationcoeff_allN_;
    worder.embedded_whichN_combinationcoeff_ = embedded_whichN_combinationcoeff_;
    worder.embedded_ref_fillback_ = embedded_ref_fillback_;

    worder.separate_averageall_ref_ = separate_averageall_ref_;
    worder.separate_fullres_coilmap_ = separate_fullres_coilmap_;
    worder.separate_same_combinationcoeff_allN_ = separate_same_combinationcoeff_allN_;
    worder.separate_whichN_combinationcoeff_ = separate_whichN_combinationcoeff_;
}

template <typename T> 
void gtPlusReconWorkOrder3DT<T>::printInfo(std::ostream& os) const
{
    using namespace std;
    BaseClass::printInfo(os);

    GADGET_OSTREAM_PRINT(os, recon_kspace_needed_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, coil_compression_);
    GADGET_OSTREAM_PRINT(os, same_coil_compression_coeff_allN_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, no_acceleration_averageall_ref_);
    GADGET_OSTREAM_PRINT(os, no_acceleration_same_combinationcoeff_allN_);
    GADGET_OSTREAM_PRINT(os, no_acceleration_whichN_combinationcoeff_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, embedded_averageall_ref_);
    GADGET_OSTREAM_PRINT(os, embedded_fullres_coilmap_);
    GADGET_OSTREAM_PRINT(os, embedded_same_combinationcoeff_allN_);
    GADGET_OSTREAM_PRINT(os, embedded_whichN_combinationcoeff_);
    GADGET_OSTREAM_PRINT(os, embedded_ref_fillback_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, separate_averageall_ref_);
    GADGET_OSTREAM_PRINT(os, separate_fullres_coilmap_);
    GADGET_OSTREAM_PRINT(os, separate_same_combinationcoeff_allN_);
    GADGET_OSTREAM_PRINT(os, separate_whichN_combinationcoeff_);
}

template <typename T> 
void gtPlusReconWorkOrder3DT<T>::print(std::ostream& os) const
{
    using namespace std;
    os << "-------------- gtPlusReconWorkOrder3DT ---------------" << endl;
    printInfo(os);
    os << "------------------------------------------------------" << endl;
}

}}
