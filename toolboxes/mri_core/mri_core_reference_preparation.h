
/** \file   mri_core_reference_preparation.h
    \brief  Implementation different calibration reference data preparation strategy for 2D and 3D MRI
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "GadgetronTimer.h"
#include "hoNDArray.h"
#include "mri_core_recon_para.h"

// #include "gtPlusIOAnalyze.h"

namespace Gadgetron
{
    /// ------------------------------------------------------------------------
    /// define the reference preparation base class
    /// ------------------------------------------------------------------------
    template <typename T>
    class EXPORTMRICORE referencePreparer
    {
    public:

        referencePreparer();
        virtual ~referencePreparer();

        /// ref: [RO E1 E2 CHA N S SLC]
        /// traj: [TRAJ E1 E2 CHA N S SLC], can be empty
        /// after the ref preparation, ref_calib_ are filled
        /// ref_calib_ is used for calibration and kernel estimation
        virtual void prep_ref(const hoNDArray<T>& ref, const hoNDArray<float>& traj) = 0;

        /// dump the class status
        virtual void dump(std::ostream& os) const;

        /// reference data for calibration
        /// [RO E1 E2 CHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_calib_;

        /// reference data ready for coil map estimation
        /// [RO E1 E2 CHA Nor1 Sor1 SLC]
        hoNDArray<T> ref_coil_map_;

        /// calibration mode
        ismrmrdCALIBMODE calib_mode_;

        /// sampled range along RO, E1, E2 (for asymmetric echo and partial fourier)
        /// min, max and center
        SamplingLimit sampling_limits_[3];

        // -----------------------------------------------
        /// for debug 
        // -----------------------------------------------

        /// clock for timing
        Gadgetron::GadgetronTimer gt_timer1_;
        Gadgetron::GadgetronTimer gt_timer2_;
        Gadgetron::GadgetronTimer gt_timer3_;

        bool perform_timing_;

        /// exporter
        // Gadgetron::gtPlus::gtPlusIOAnalyze gt_exporter_;

        /// debug folder
        std::string debug_folder_;

        /// whether to print out more information
        bool verbose_;
    };

    /// ------------------------------------------------------------------------
    /// reference preparation for cartesian sampling
    /// ------------------------------------------------------------------------
    template <typename T>
    class EXPORTMRICORE referenceCartesianPreparer : public referencePreparer<T>
    {
    public:

        typedef referencePreparer<T> BaseClass;
        typedef typename realType<T>::Type value_type;

        referenceCartesianPreparer();
        virtual ~referenceCartesianPreparer();

        /// besides filling the ref_calib_, the ref_coil_map_ is also filled
        /// ref_coil_map_ is used for coil map estimation and is padded to the reconed image size
        virtual void prep_ref(const hoNDArray<T>& ref, const hoNDArray<float>& traj);

        /// dump the class status
        virtual void dump(std::ostream& os) const;

        using BaseClass::ref_calib_;
        using BaseClass::ref_coil_map_;
        using BaseClass::calib_mode_;
        using BaseClass::sampling_limits_;
        using BaseClass::verbose_;

        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::perform_timing_;
        // using BaseClass::gt_exporter_;
        using BaseClass::debug_folder_;

        /// ref filter for coil map estimation
        /// these filters will be created when they are first used
        hoNDArray<T> filter_RO_ref_coi_map_;
        hoNDArray<T> filter_E1_ref_coi_map_;
        hoNDArray<T> filter_E2_ref_coi_map_;

        /// reconed image size
        size_t recon_RO_;
        size_t recon_E1_;
        size_t recon_E2_;

        /// whether to average all N for ref generation
        bool average_all_ref_N_;
        /// whether to average all S for ref generation
        bool average_all_ref_S_;
        /// whether to filter ref for coil map estimation
        bool filter_ref_coil_map_;
        /// for interleaved mode, whether the interleaved dimension is along N
        /// if false, then the interleaved dimension is along S
        bool interleave_dim_along_N_;

    protected:
        /// prepare functions for every calibration mode

        /// first, averaging over N or S
        /// second, filter the ref_calib_ to perpare ref_coil_map_
        void prep_ref_no_acceleration(const hoNDArray<T>& ref);

        /// first, averaging over N or S, depending the interleaved dimension
        /// second, filter the ref_calib_ to perpare ref_coil_map_
        /// them, crop the unsampled region from ref_calib_, in case of partial fourier or asymmetric echo is used
        /// after cropping, ref_calib_ only includes the fully sampled data and ready for calibration
        void prep_ref_interleaved(const hoNDArray<T>& ref);

        /// first, averaging over N or S
        /// second, detect sampled region for reference along E1 and E2
        /// set up and filter the fully sampled ref data
        /// padding the filtered ref data to have the recon image size; this is the ref_coil_map_
        /// crop the ref data to make ref_calib_ ready
        void prep_ref_embedded(const hoNDArray<T>& ref);

        /// first, averaging over N or S
        /// second, detect sampled region for reference
        /// set up and filter the fully sampled ref data
        /// padding the filtered ref data to have the recon image size; this is the ref_coil_map_
        /// crop the ref data to make ref_calib_ ready
        void prep_ref_separate(const hoNDArray<T>& ref);

        /// average across N and/or S
        /// results are stored in ref_calib_
        void average_ref(const hoNDArray<T>& ref);
    };

    /// ------------------------------------------------------------------------
    /// reference preparation for non-cartesian sampling
    /// ------------------------------------------------------------------------
    template <typename T>
    class EXPORTMRICORE referenceNonCartesianPreparer : public referencePreparer<T>
    {
    public:

        typedef referencePreparer<T> BaseClass;
        typedef typename realType<T>::Type value_type;

        referenceNonCartesianPreparer();
        virtual ~referenceNonCartesianPreparer();

        /// ref: [RO E1 E2 CHA N S SLC]
        virtual void prep_ref(const hoNDArray<T>& ref, const hoNDArray<float>& traj);

        /// dump the class status
        virtual void dump(std::ostream& os) const;

        using BaseClass::ref_calib_;
        using BaseClass::ref_coil_map_;
        using BaseClass::calib_mode_;
        using BaseClass::sampling_limits_;
        using BaseClass::verbose_;

        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::perform_timing_;
        // using BaseClass::gt_exporter_;
        using BaseClass::debug_folder_;

        /// trajectory for calibration data
        /// [TRAJ E1 E2 CHA Nor1 Sor1 SLC]
        hoNDArray<float> traj_calib_;

        /// the default implementation in this calls maks ref_coil_map_ a shared-memory duplicate of ref_calib_, so is the traj_coil_map_
        hoNDArray<float> traj_coil_map_;

        /// whether to combine (not average) all N for ref generation
        bool combine_all_ref_N_;
        /// whether to combine (not average) all S for ref generation
        bool combine_all_ref_S_;
        /// for interleaved mode, whether the interleaved dimension is along N
        /// if false, then the interleaved dimension is along S
        bool interleave_dim_along_N_;

    protected:
        /// combine ref and traj array along N and S, controlled by combineN and combineS
        void combine_along_N_S(const hoNDArray<T>& ref, const hoNDArray<float>& traj, bool combineN, bool combinedS);
    };
}
