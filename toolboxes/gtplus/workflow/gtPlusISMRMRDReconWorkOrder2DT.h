/** \file   gtPlusISMRMRDReconWorkOrder2DT.h
    \brief  Define the GtPlus reconstruction workorder and parameters for 2DT reconstruction
    \author Hui Xue
*/

#pragma once

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorkOrder2DT : public gtPlusReconWorkOrder<T>
{
public:

    typedef gtPlusReconWorkOrder<T> BaseClass;

    gtPlusReconWorkOrder2DT();
    virtual ~gtPlusReconWorkOrder2DT();

    virtual bool reset();

    virtual bool enforceConsistency(ISMRMRDDIM& lastDim);
    virtual void duplicate(gtPlusReconWorkOrder2DT<T>& worder);

    virtual void printInfo(std::ostream& os) const;
    virtual void print(std::ostream& os) const;

    // kspace_: [RO E1 CHA N S], for 2D recon, N can be 1
    // ref_: [RO E1 CHA M S], M can equal to N or 1 or others
    // fullkspace_: [RO E1 CHA N S]
    // complexIm_: [RO E1 N S], after coil combination or [RO E1 num_channels_res_ N S] if num_channels_res_ > 1
    // coilMap_: [RO E1 CHA 1 or N S]
    // gfactor_: [RO E1 CHA 1 or N S]

    // the fifth dimension can be SLC or SET or others

    // default behavior
    // a) the coil compression coefficients are computed once across all S
    // b) the kernel or coil sensitivity are estimated for every S

    // embedded mode
    // a) perform recon and estimate kernel/coil sensitivity for every 2D kspace [RO E1 CHA]
    // b) coil combination uses different coil maps for every S
    // c) if the kspace recon is performed, the coil combination map is reestimated on the fullkspace for every 2D images
    // d) the ref lines are filled back to fullkspace_

    // separate mode
    // a) perform recon and estimate kernel/coil sensitivity for every 2D kspace [RO E1 CHA] if M==N
    // b) if M==1, the kernel is only estimated once for every S
    // c) coil combination uses different coil maps for every S
    // d) if the kspace recon is performed, the coil combination map is reestimated on the fullkspace for every 2D images

    // interleave
    // a) the average-all ref is used
    // b) kernel/coil sensitivity is estimated once for every S

    using BaseClass::data_;
    using BaseClass::ref_;
    using BaseClass::ref_recon_;
    using BaseClass::ref_coil_map_;
    using BaseClass::kspace_initial_;
    using BaseClass::CalibMode_;
    using BaseClass::InterleaveDim_;
    using BaseClass::acceFactorE1_;
    using BaseClass::acceFactorE2_;
    using BaseClass::num_channels_res_;
    using BaseClass::coilMap_; // [RO E1 dstCHA 1 or N S]

    using BaseClass::fullkspace_; // [RO E1 dstCHA N S]
    using BaseClass::complexIm_; // [RO E1 N S]
    using BaseClass::recon_time_stamp_; // [1 1 1 N S]
    using BaseClass::recon_physio_time_stamp_; // [1 1 1 N S]

    using BaseClass::fullkspace_second_; // [RO E1 dstCHA N S]
    using BaseClass::complexIm_second_; // [RO E1 N S]
    using BaseClass::recon_time_stamp_second_; // [1 1 1 N S]
    using BaseClass::recon_physio_time_stamp_second_; // [1 1 1 N S]

    using BaseClass::gfactor_; // [RO E1 1 or N S]
    using BaseClass::wrap_around_map_; // [RO E1 2 1 or N S]

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

    using BaseClass::CloudComputing_;
    using BaseClass::CloudSize_;
    using BaseClass::gt_cloud_;

    // for 2DT
    using BaseClass::kernel_; // [RO E1 srcCHA dstCHA dstE1 1 or N S]
    using BaseClass::kernelIm_; // [RO E1 srcCHA dstCHA 1 or N S]
    using BaseClass::unmixingCoeffIm_; // [RO E1 srcCHA 1 or N S]
    using BaseClass::coilCompressionCoef_; // [dstCHA srcCHA] matrixes

    // parameters to change the default behavior

    // if true, the actual full kspace is computed, not only the coil combined complex images
    bool recon_kspace_needed_;

    // if true, no coil compression will be performed
    bool coil_compression_;
    // if true, the same coil compression coefficient is computed for all S
    bool same_coil_compression_coeff_allS_;

    // no acceleration
    // if true, the average of all M ref will be used
    // the coil sensitivity will be only estimed once for all N
    bool no_acceleration_averageall_ref_;
    // number of modes kept for ref data
    int no_acceleration_ref_numOfModes_;
    // if true, the same coil combination coefficients will be used for all S
    bool no_acceleration_same_combinationcoeff_allS_;
    // if no_acceleration_same_combinationcoeff_allS_==true, select the S for coil combination coefficient estimation
    size_t no_acceleration_whichS_combinationcoeff_;

    // embedded mode
    // if true, the average of all M ref will be used
    // the kernel/sensitivity will be only estimed once for all N
    bool embedded_averageall_ref_;
    // number of modes kept for ref data
    int embedded_ref_numOfModes_;
    // if true, the coil map will be estimated from the fullkspace_
    bool embedded_fullres_coilmap_;
    // if embedded_averageall_ref_==true && embedded_fullres_coilmap_==true, whether to select the highest signal frame to compute full res coil map
    // if false, the averageall image will be used to compute full res coil map
    bool embedded_fullres_coilmap_useHighestSignal_;
    // if true, the same coil combination coefficients will be used for all S
    bool embedded_same_combinationcoeff_allS_;
    // if embedded_same_combinationcoeff_allS_==true, select the S for coil combination coefficient estimation
    size_t embedded_whichS_combinationcoeff_;
    // if true, the ref lines will be filled back to fullkspace
    bool embedded_ref_fillback_;

    // separate mode
    // if true, the average of all M ref will be used
    // the kernel/sensitivity will be only estimed once for every S
    bool separate_averageall_ref_;
    // number of modes kept for ref data
    int separate_ref_numOfModes_;
    // if true, the coil map will be estimated from the fullkspace_
    bool separate_fullres_coilmap_;
    // if true, the same coil combination coefficients will be used for all S
    bool separate_same_combinationcoeff_allS_;
    // if separate_same_combinationcoeff_allS_==true, select the S for coil combination coefficient estimation
    size_t separate_whichS_combinationcoeff_;

    // interleaved mode
    // if true, the same coil combination coefficients will be used for all S
    bool interleaved_same_combinationcoeff_allS_;
    // if separate_same_combinationcoeff_allS_==true, select the S for coil combination coefficient estimation
    size_t interleaved_whichS_combinationcoeff_;
    // number of modes kept for ref data
    int interleaved_ref_numOfModes_;
};

template <typename T> 
gtPlusReconWorkOrder2DT<T>::gtPlusReconWorkOrder2DT() : BaseClass()
{
    coil_compression_ = true;
    same_coil_compression_coeff_allS_ = false;

    no_acceleration_averageall_ref_ = true;
    no_acceleration_ref_numOfModes_ = 3;
    no_acceleration_same_combinationcoeff_allS_ = false;
    no_acceleration_whichS_combinationcoeff_ = 0;

    embedded_averageall_ref_ = false;
    embedded_ref_numOfModes_ = 3;
    embedded_fullres_coilmap_ = true;
    embedded_fullres_coilmap_useHighestSignal_ = false;
    embedded_same_combinationcoeff_allS_ = false;
    embedded_whichS_combinationcoeff_ = false;
    embedded_ref_fillback_ = true;

    separate_averageall_ref_ = false;
    separate_ref_numOfModes_ = 3;
    separate_fullres_coilmap_ = true;
    separate_same_combinationcoeff_allS_ = false;
    separate_whichS_combinationcoeff_ = false;

    interleaved_same_combinationcoeff_allS_ = false;
    interleaved_whichS_combinationcoeff_ = false;
    interleaved_ref_numOfModes_ = 0;
}

template <typename T> 
gtPlusReconWorkOrder2DT<T>::~gtPlusReconWorkOrder2DT()
{
}

template <typename T> 
bool gtPlusReconWorkOrder2DT<T>::reset()
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
        GADGET_ERROR_MSG("Errors in gtPlusReconWorkOrder2DT<T>::reset() ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorkOrder2DT<T>::enforceConsistency(ISMRMRDDIM& lastDim)
{
    if ( lastDim == DIM_Slice )
    {
        same_coil_compression_coeff_allS_ = false;
        no_acceleration_same_combinationcoeff_allS_ = false;
        embedded_same_combinationcoeff_allS_ = false;
        separate_same_combinationcoeff_allS_ = false;
        interleaved_same_combinationcoeff_allS_ = false;
    }

    return true;
}

template <typename T> 
void gtPlusReconWorkOrder2DT<T>::duplicate(gtPlusReconWorkOrder2DT<T>& worder)
{
    BaseClass::duplicate(worder);

    worder.recon_kspace_needed_ = recon_kspace_needed_;

    worder.coil_compression_ = coil_compression_;
    worder.same_coil_compression_coeff_allS_ = same_coil_compression_coeff_allS_;

    worder.no_acceleration_averageall_ref_ = no_acceleration_averageall_ref_;
    worder.no_acceleration_ref_numOfModes_ = no_acceleration_ref_numOfModes_;
    worder.no_acceleration_same_combinationcoeff_allS_ = no_acceleration_same_combinationcoeff_allS_;
    worder.no_acceleration_whichS_combinationcoeff_ = no_acceleration_whichS_combinationcoeff_;

    worder.embedded_averageall_ref_ = embedded_averageall_ref_;
    worder.embedded_ref_numOfModes_ = embedded_ref_numOfModes_;
    worder.embedded_fullres_coilmap_ = embedded_fullres_coilmap_;
    worder.embedded_fullres_coilmap_useHighestSignal_ = embedded_fullres_coilmap_useHighestSignal_;
    worder.embedded_same_combinationcoeff_allS_ = embedded_same_combinationcoeff_allS_;
    worder.embedded_whichS_combinationcoeff_ = embedded_whichS_combinationcoeff_;
    worder.embedded_ref_fillback_ = embedded_ref_fillback_;

    worder.separate_averageall_ref_ = separate_averageall_ref_;
    worder.separate_ref_numOfModes_ = separate_ref_numOfModes_;
    worder.separate_fullres_coilmap_ = separate_fullres_coilmap_;
    worder.separate_same_combinationcoeff_allS_ = separate_same_combinationcoeff_allS_;
    worder.separate_whichS_combinationcoeff_ = separate_whichS_combinationcoeff_;

    worder.interleaved_same_combinationcoeff_allS_ = interleaved_same_combinationcoeff_allS_;
    worder.interleaved_whichS_combinationcoeff_ = interleaved_whichS_combinationcoeff_;
    worder.interleaved_ref_numOfModes_ = interleaved_ref_numOfModes_;
}

template <typename T> 
void gtPlusReconWorkOrder2DT<T>::printInfo(std::ostream& os) const
{
    using namespace std;
    BaseClass::printInfo(os);

    GADGET_OSTREAM_PRINT(os, recon_kspace_needed_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, coil_compression_);
    GADGET_OSTREAM_PRINT(os, same_coil_compression_coeff_allS_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, no_acceleration_averageall_ref_);
    GADGET_OSTREAM_PRINT(os, no_acceleration_ref_numOfModes_);
    GADGET_OSTREAM_PRINT(os, no_acceleration_same_combinationcoeff_allS_);
    GADGET_OSTREAM_PRINT(os, no_acceleration_whichS_combinationcoeff_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, embedded_averageall_ref_);
    GADGET_OSTREAM_PRINT(os, embedded_ref_numOfModes_);
    GADGET_OSTREAM_PRINT(os, embedded_fullres_coilmap_);
    GADGET_OSTREAM_PRINT(os, embedded_fullres_coilmap_useHighestSignal_);
    GADGET_OSTREAM_PRINT(os, embedded_same_combinationcoeff_allS_);
    GADGET_OSTREAM_PRINT(os, embedded_whichS_combinationcoeff_);
    GADGET_OSTREAM_PRINT(os, embedded_ref_fillback_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, separate_averageall_ref_);
    GADGET_OSTREAM_PRINT(os, separate_ref_numOfModes_);
    GADGET_OSTREAM_PRINT(os, separate_fullres_coilmap_);
    GADGET_OSTREAM_PRINT(os, separate_same_combinationcoeff_allS_);
    GADGET_OSTREAM_PRINT(os, separate_whichS_combinationcoeff_);
    os << std::endl;
    GADGET_OSTREAM_PRINT(os, interleaved_same_combinationcoeff_allS_);
    GADGET_OSTREAM_PRINT(os, interleaved_whichS_combinationcoeff_);
    GADGET_OSTREAM_PRINT(os, interleaved_ref_numOfModes_);
}

template <typename T> 
void gtPlusReconWorkOrder2DT<T>::print(std::ostream& os) const
{
    using namespace std;
    os << "-------------- gtPlusReconWorkOrder2DT ---------------" << endl;
    printInfo(os);
    os << "------------------------------------------------------" << endl;
}

}}
