
/** \file   mri_core_reference_preparation.cpp
    \brief  Implementation different calibration reference data preparation strategy for 2D and 3D MRI
    \author Hui Xue
*/

#include "mri_core_reference_preparation.h"
#include "mri_core_utility.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron
{
    template<typename T> 
    referencePreparer<T>::referencePreparer()
    {
        calib_mode_ = ISMRMRD_noacceleration;

        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);
        perform_timing_ = false;

        verbose_ = false;
    }

    template<typename T>
    referencePreparer<T>::~referencePreparer()
    {
    }

    template<typename T>
    void referencePreparer<T>::dump(std::ostream& os) const
    {
        using namespace std;
        os << "-------------------------------------------------------------------------------" << endl;
        os << "referencePreparer will prepare the ref_calib_ array used for calibration and kernel estimation ... " << endl;
        os << "pref_ref should be implemented ... " << endl;
        os << "-------------------------------------------------------------------------------" << endl;
    }

    /// ------------------------------------------------------------------------

    template<typename T>
    referenceCartesianPreparer<T>::referenceCartesianPreparer() : BaseClass()
    {
        recon_RO_ = 0;
        recon_E1_ = 0;
        recon_E2_ = 0;

        average_all_ref_N_ = true;
        average_all_ref_S_ = false;
        filter_ref_coil_map_ = true;
        interleave_dim_along_N_ = true;
    }

    template<typename T>
    referenceCartesianPreparer<T>::~referenceCartesianPreparer()
    {
    }

    template<typename T>
    void referenceCartesianPreparer<T>::average_ref(const hoNDArray<T>& ref)
    {
        try
        {
            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N = ref.get_size(4);
            size_t S = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            // if interleaved mode is used, we have to count the sampled times for every kspace location, in case the kspace has cartesian irregular sampling
            if (calib_mode_ == ISMRMRD_interleaved)
            {
                // average over N
                if (interleave_dim_along_N_)
                {
                    if (N > 1)
                    {
                        // averaging by counting the sampled times
                        Gadgetron::average_kspace_across_N(ref, ref_calib_);
                    }
                    else
                    {
                        ref_calib_ = ref;
                    }

                    if (!debug_folder_.empty())
                    {
                        // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_interleaved_ave_N");
                    }

                    if (average_all_ref_S_) // since S is not the interleaved dimension, the sum can be used
                    {
                        if (S > 1)
                        {
                            hoNDArray<T> ref_recon_buf;
                            Gadgetron::sum_over_dimension(ref_calib_, ref_recon_buf, 5);
                            Gadgetron::scal((value_type)(1.0 / S), ref_recon_buf);
                            ref_calib_ = ref_recon_buf;
                        }

                        if (!debug_folder_.empty())
                        {
                            // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_interleaved_ave_S");
                        }
                    }
                }
                else
                    // average over S
                {
                    if (S > 1)
                    {
                        // averaging by counting the sampled times
                        Gadgetron::average_kspace_across_S(ref, ref_calib_);
                    }
                    else
                    {
                        ref_calib_ = ref;
                    }

                    if (!debug_folder_.empty())
                    {
                        // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_interleaved_ave_S");
                    }

                    if (average_all_ref_N_)
                    {
                        if (N > 1)
                        {
                            hoNDArray<T> ref_recon_buf;
                            Gadgetron::sum_over_dimension(ref_calib_, ref_recon_buf, 4);
                            Gadgetron::scal((value_type)(1.0 / N), ref_recon_buf);
                            ref_calib_ = ref_recon_buf;
                        }

                        if (!debug_folder_.empty())
                        {
                            // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_interleaved_ave_N");
                        }
                    }
                }
            }
            else
            {
                if (average_all_ref_N_)
                {
                    if (N > 1)
                    {
                        Gadgetron::sum_over_dimension(ref, ref_calib_, (size_t)4);
                        Gadgetron::scal((value_type)(1.0 / N), ref_calib_);

                        if (!debug_folder_.empty())
                        {
                            // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_after_averaging_N");
                        }
                    }
                    else
                    {
                        ref_calib_ = ref;
                    }
                }
                else
                {
                    ref_calib_ = ref;
                }

                if (average_all_ref_S_)
                {
                    if (S > 1)
                    {
                        hoNDArray<T> ref_recon_buf;
                        Gadgetron::sum_over_dimension(ref_calib_, ref_recon_buf, 5);
                        Gadgetron::scal((value_type)(1.0 / S), ref_recon_buf);
                        ref_calib_ = ref_recon_buf;

                        if (!debug_folder_.empty())
                        {
                            // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_after_averaging_S");
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in referenceCartesianPreparer<T>::average_ref(...) ... ");
        }
    }

    template<typename T>
    void referenceCartesianPreparer<T>::prep_ref_no_acceleration(const hoNDArray<T>& ref)
    {
        try
        {
            size_t RO  = ref.get_size(0);
            size_t E1  = ref.get_size(1);
            size_t E2  = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N   = ref.get_size(4);
            size_t S   = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            this->average_ref(ref);

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_noacceleration");
            }

            if (filter_ref_coil_map_)
            {
                if (E2 > 1)
                {
                    Gadgetron::apply_kspace_filter_ROE1E2(ref_calib_, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_, ref_coil_map_);
                }
                else
                {
                    Gadgetron::apply_kspace_filter_ROE1(ref_calib_, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, ref_coil_map_);
                }
            }

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_coil_map_, debug_folder_ + "ref_coil_map_noacceleration");
            }

            if (recon_RO_ > RO || recon_E1_ > E1 || recon_E2_ > E2)
            {
                hoNDArray<T> ref_recon_buf;
                Gadgetron::pad(recon_RO_, recon_E1_, recon_E2_, &ref_coil_map_, &ref_recon_buf);
                ref_coil_map_ = ref_recon_buf;

                if (!debug_folder_.empty())
                {
                    // gt_exporter_.exportArrayComplex(ref_coil_map_, debug_folder_ + "ref_coil_map_padded_noacceleration");
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in referenceCartesianPreparer<T>::prep_ref_no_acceleration(...) ... ");
        }
    }

    template<typename T>
    void referenceCartesianPreparer<T>::prep_ref_interleaved(const hoNDArray<T>& ref)
    {
        try
        {
            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N = ref.get_size(4);
            size_t S = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            this->average_ref(ref);

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_interleaved");
            }

            // perform Ref filtering
            if (filter_ref_coil_map_)
            {
                if (E2 > 1)
                {
                    Gadgetron::apply_kspace_filter_ROE1E2(ref_calib_, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_, ref_coil_map_);
                }
                else
                {
                    Gadgetron::apply_kspace_filter_ROE1(ref_calib_, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, ref_coil_map_);
                }
            }

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_coil_map_, debug_folder_ + "ref_coil_map_interleaved");
            }

            // pad the ref_coil_map into the data array
            if (recon_RO_ > RO || recon_E1_ > E1 || recon_E1_ > E2)
            {
                hoNDArray<T> ref_recon_buf;
                Gadgetron::pad(recon_RO_, recon_E1_, recon_E2_, &ref_coil_map_, &ref_recon_buf);
                ref_coil_map_ = ref_recon_buf;

                if (!debug_folder_.empty())
                {
                    // gt_exporter_.exportArrayComplex(ref_coil_map_, debug_folder_ + "ref_coil_map_padded_interleaved");
                }
            }

            // crop the ref_calib_ to make it ready for calibration, because only sampled kspace region is needed in calibration
            size_t len_RO = sampling_limits_[0].max_ - sampling_limits_[0].min_ + 1;
            size_t len_E1 = sampling_limits_[1].max_ - sampling_limits_[1].min_ + 1;
            size_t len_E2 = (E2 > 1) ? (sampling_limits_[2].max_ - sampling_limits_[2].min_ + 1) : 1;

            if (len_RO < RO || len_E1 < E1 || len_E2 < E2)
            {
                vector_td<size_t, 3> crop_offset;
                crop_offset[0] = sampling_limits_[0].min_;
                crop_offset[1] = sampling_limits_[1].min_;
                crop_offset[2] = (E2>1) ? sampling_limits_[2].min_ : 0;

                vector_td<size_t, 3> crop_size;
                crop_size[0] = len_RO;
                crop_size[1] = len_E1;
                crop_size[2] = len_E2;

                hoNDArray<T> ref_recon_buf;
                Gadgetron::crop(crop_offset, crop_size, &ref_calib_, &ref_recon_buf);
                ref_calib_ = ref_recon_buf;

                if (!debug_folder_.empty())
                {
                    // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_after_crop_interleaved");
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in referenceCartesianPreparer<T>::prep_ref_interleaved(...) ... ");
        }
    }

    template<typename T>
    void referenceCartesianPreparer<T>::prep_ref_embedded(const hoNDArray<T>& ref)
    {
        try
        {
            // ref lines should be in the center of kspace
            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N = ref.get_size(4);
            size_t S = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            this->average_ref(ref);

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_after_averaging_embedded");
            }

            // detect sampled region in ref
            size_t start_E1(0), end_E1(0);
            Gadgetron::detect_sampled_region_E1(ref_calib_, start_E1, end_E1);

            size_t start_E2(0), end_E2(0);
            if (E2 > 1)
            {
                Gadgetron::detect_sampled_region_E2(ref_calib_, start_E2, end_E2);
            }

            hoNDArray<T> ref_recon_buf;

            // crop the ref, along E1 and E2
            vector_td<size_t, 3> crop_offset;
            crop_offset[0] = 0;
            crop_offset[1] = start_E1;
            crop_offset[2] = start_E2;

            vector_td<size_t, 3> crop_size;
            crop_size[0] = RO;
            crop_size[1] = end_E1 - start_E1 + 1;
            crop_size[2] = end_E2 - start_E2 + 1;

            Gadgetron::crop(crop_offset, crop_size, &ref_calib_, &ref_recon_buf);

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_recon_buf, debug_folder_ + "ref_recon_buf_after_crop_embedded");
            }

            ref_calib_ = ref_recon_buf;

            // generate ref filter
            if (filter_ref_coil_map_)
            {
                SamplingLimit sampleLimE1, sampleLimE2;
                sampleLimE1.min_ = 0;
                sampleLimE1.max_ = ref_recon_buf.get_size(1)-1;

                sampleLimE2.min_ = 0;
                sampleLimE2.max_ = ref_recon_buf.get_size(2) - 1;

                Gadgetron::generate_ref_filter_for_coil_map(ref_recon_buf, sampling_limits_[0], sampleLimE1, sampleLimE2, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_);

                if (!debug_folder_.empty())
                {
                    // gt_exporter_.exportArrayComplex(filter_RO_ref_coi_map_, debug_folder_ + "filter_RO_ref_coi_map_embedded");
                    // gt_exporter_.exportArrayComplex(filter_E1_ref_coi_map_, debug_folder_ + "filter_E1_ref_coi_map_embedded");
                    if (E2 > 1) 
                    {
                        // gt_exporter_.exportArrayComplex(filter_E2_ref_coi_map_, debug_folder_ + "filter_E2_ref_coi_map_embedded"); 
                    }
                }

                if (E2 > 1)
                {
                    Gadgetron::apply_kspace_filter_ROE1E2(ref_calib_, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_, ref_recon_buf);
                }
                else
                {
                    Gadgetron::apply_kspace_filter_ROE1(ref_calib_, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, ref_recon_buf);
                }

                if (!debug_folder_.empty())
                {
                    // gt_exporter_.exportArrayComplex(ref_recon_buf, debug_folder_ + "ref_recon_buf_ref_filtered");
                }
            }

            // pad ref_calib_ into the data size
            Gadgetron::pad(recon_RO_, recon_E1_, recon_E2_, &ref_recon_buf, &ref_coil_map_);

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_coil_map_, debug_folder_ + "ref_coil_map_embedded");
            }

            // further crop along RO
            if (sampling_limits_[0].min_ > 0 || sampling_limits_[0].max_ < RO - 1)
            {
                crop_offset[0] = sampling_limits_[0].min_;
                crop_offset[1] = 0;
                crop_offset[2] = 0;

                crop_size[0] = sampling_limits_[0].max_ - sampling_limits_[0].min_ + 1;
                crop_size[1] = ref_calib_.get_size(1) - 1;
                crop_size[2] = ref_calib_.get_size(2) - 1;

                Gadgetron::crop(crop_offset, crop_size, &ref_calib_, &ref_recon_buf);

                if (!debug_folder_.empty())
                {
                    // gt_exporter_.exportArrayComplex(ref_recon_buf, debug_folder_ + "ref_recon_buf_after_crop_RO_embedded");
                }

                ref_calib_ = ref_recon_buf;
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in referenceCartesianPreparer<T>::prep_ref_embedded(...) ... ");
        }
    }

    template<typename T>
    void referenceCartesianPreparer<T>::prep_ref_separate(const hoNDArray<T>& ref)
    {
        try
        {
            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N = ref.get_size(4);
            size_t S = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            this->average_ref(ref);

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_after_averaging_separate");
            }

            hoNDArray<T> ref_recon_buf;

            // detect sampled region in ref
            size_t start_E1(0), end_E1(0);
            Gadgetron::detect_sampled_region_E1(ref, start_E1, end_E1);

            size_t start_E2(0), end_E2(0);
            if (E2 > 1)
            {
                Gadgetron::detect_sampled_region_E2(ref, start_E2, end_E2);
            }

            // crop the ref_calib_, along E1 and E2
            vector_td<size_t, 3> crop_offset;
            crop_offset[0] = sampling_limits_[0].min_;
            crop_offset[1] = start_E1;
            crop_offset[2] = start_E2;

            vector_td<size_t, 3> crop_size;
            crop_size[0] = sampling_limits_[0].max_ - sampling_limits_[0].min_ + 1;
            crop_size[1] = end_E1 - start_E1 + 1;
            crop_size[2] = end_E2 - start_E2 + 1;

            Gadgetron::crop(crop_offset, crop_size, &ref_calib_, &ref_recon_buf);
            ref_calib_ = ref_recon_buf;

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_after_crop_separate");
            }

            // create filter if needed
            if (filter_ref_coil_map_)
            {
                if ( filter_RO_ref_coi_map_.get_size(0) != RO 
                    || filter_E1_ref_coi_map_.get_size(0) != ref_calib_.get_size(1) 
                    || ((E2 > 1) && (filter_E2_ref_coi_map_.get_size(0) != ref_calib_.get_size(2))) )
                {
                    SamplingLimit sE1;
                    sE1.min_ = 0;
                    sE1.max_ = ref_calib_.get_size(1) - 1;

                    SamplingLimit sE2;
                    sE2.min_ = 0;
                    sE2.max_ = ref_calib_.get_size(2) - 1;

                    Gadgetron::generate_ref_filter_for_coil_map(ref_calib_, sampling_limits_[0], sE1, sE2, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_);
                }

                if (!debug_folder_.empty())
                {
                    // gt_exporter_.exportArrayComplex(filter_RO_ref_coi_map_, debug_folder_ + "filter_RO_ref_coi_map_separate");
                    // gt_exporter_.exportArrayComplex(filter_E1_ref_coi_map_, debug_folder_ + "filter_E1_ref_coi_map_separate");
                    if (E2 > 1) 
                    {
                        // gt_exporter_.exportArrayComplex(filter_E2_ref_coi_map_, debug_folder_ + "filter_E2_ref_coi_map_separate"); 
                    }
                }
            }

            // filter the ref_coil_map_
            ref_coil_map_ = ref_calib_;

            if (filter_ref_coil_map_)
            {
                if (E2 > 1)
                {
                    Gadgetron::apply_kspace_filter_ROE1E2(ref_coil_map_, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_, ref_recon_buf);
                }
                else
                {
                    Gadgetron::apply_kspace_filter_ROE1(ref_coil_map_, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, ref_recon_buf);
                }

                ref_coil_map_ = ref_recon_buf;
            }

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_coil_map_, debug_folder_ + "ref_coil_map_separate");
            }

            // pad the ref_coil_map into the data array
            Gadgetron::pad(recon_RO_, recon_E1_, recon_E2_, &ref_coil_map_, &ref_recon_buf);
            ref_coil_map_ = ref_recon_buf;

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_coil_map_, debug_folder_ + "ref_coil_map_padded_separate");
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in referenceCartesianPreparer<T>::prep_ref_separate(...) ... ");
        }
    }

    template<typename T>
    void referenceCartesianPreparer<T>::prep_ref(const hoNDArray<T>& ref, const hoNDArray<float>& traj)
    {
        try
        {
            size_t RO  = ref.get_size(0);
            size_t E1  = ref.get_size(1);
            size_t E2  = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N   = ref.get_size(4);
            size_t S   = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            // if needed, generate the ref filter
            if (calib_mode_ == ISMRMRD_noacceleration
                || calib_mode_ == ISMRMRD_interleaved)
            {
                if (filter_ref_coil_map_)
                {
                    if (filter_RO_ref_coi_map_.get_size(0) != RO || filter_E1_ref_coi_map_.get_size(0) != E1 || ((E2 > 1) && (filter_E2_ref_coi_map_.get_size(0) != E2)))
                    {
                        Gadgetron::generate_ref_filter_for_coil_map(ref, sampling_limits_[0], sampling_limits_[1], sampling_limits_[2],
                            filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_);
                    }

                    if (!debug_folder_.empty())
                    {
                        // gt_exporter_.exportArrayComplex(filter_RO_ref_coi_map_, debug_folder_ + "filter_RO_ref_coi_map");
                        // gt_exporter_.exportArrayComplex(filter_E1_ref_coi_map_, debug_folder_ + "filter_E1_ref_coi_map");
                        if (E2 > 1) 
                        { 
                            // gt_exporter_.exportArrayComplex(filter_E2_ref_coi_map_, debug_folder_ + "filter_E2_ref_coi_map"); 
                        }
                    }
                }
            }

            // -------------------------------------------------------------------------------

            if (calib_mode_ == ISMRMRD_noacceleration)
            {
                this->prep_ref_no_acceleration(ref);
            }

            // -------------------------------------------------------------------------------

            else if (calib_mode_ == ISMRMRD_interleaved)
            {
                this->prep_ref_interleaved(ref);
            }

            // -------------------------------------------------------------------------------

            else if (calib_mode_ == ISMRMRD_embedded)
            {
                this->prep_ref_embedded(ref);
            }

            // -------------------------------------------------------------------------------

            else if (calib_mode_ == ISMRMRD_separate || calib_mode_ == ISMRMRD_external || calib_mode_ == ISMRMRD_other)
            {
                this->prep_ref_separate(ref);
            }
            else
            {
                GADGET_THROW("Unrecognized calibration mode ... ");
            }

            // -------------------------------------------------------------------------------

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_after_all_prep");
            }

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_coil_map_, debug_folder_ + "ref_coil_map_after_all_prep");
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in referenceCartesianPreparer<T>::prep_ref() ... ");
        }
    }

    template<typename T>
    void referenceCartesianPreparer<T>::dump(std::ostream& os) const
    {
        using namespace std;
        os << "-------------------------------------------------------------------------------" << endl;
        os << "referenceCartesianPreparer ... " << endl;
        os << "calib_mode_ is " << Gadgetron::get_ismrmrd_calib_mode_name(calib_mode_) << endl;
        os << "recon image size is [" << recon_RO_ << " " << recon_E1_ << " " << recon_E2_ << "]" << endl;
        os << "average_all_ref_N_ is " << average_all_ref_N_ << endl;
        os << "average_all_ref_S_ is " << average_all_ref_S_ << endl;
        os << "filter_ref_coil_map_ is " << filter_ref_coil_map_ << endl;
        os << "interleave_dim_along_N_ is " << interleave_dim_along_N_ << endl;
        os << "verbose_ is " << verbose_ << endl;
        os << "-------------------------------------------------------------------------------" << endl;
    }

    template class EXPORTMRICORE referenceCartesianPreparer < std::complex<float> >;
    template class EXPORTMRICORE referenceCartesianPreparer < std::complex<double> >;

    /// ------------------------------------------------------------------------

    template<typename T>
    referenceNonCartesianPreparer<T>::referenceNonCartesianPreparer() : BaseClass()
    {
        combine_all_ref_N_ = true;
        combine_all_ref_S_ = false;
        interleave_dim_along_N_ = true;
    }

    template<typename T>
    referenceNonCartesianPreparer<T>::~referenceNonCartesianPreparer()
    {
    }

    template<typename T>
    void referenceNonCartesianPreparer<T>::combine_along_N_S(const hoNDArray<T>& ref, const hoNDArray<float>& traj, bool combineN, bool combineS)
    {
        try
        {
            if (!combineN && !combineS)
            {
                ref_calib_ = ref;
                traj_calib_ = traj;
                return;
            }

            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N = ref.get_size(4);
            size_t S = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            size_t TRAJ = traj.get_size(0);

            if (E2 > 1)
            {
                if (combineN && combineS)
                {
                    ref_calib_.create(RO, E1, E2*N*S, CHA, 1, 1, SLC);
                    traj_calib_.create(TRAJ, E1, E2*N*S, CHA, 1, 1, SLC);
                }
                else if (combineN && !combineS)
                {
                    ref_calib_.create(RO, E1, E2*N, CHA, 1, S, SLC);
                    traj_calib_.create(TRAJ, E1, E2*N, CHA, 1, S, SLC);
                }
                else if (!combineN && combineS)
                {
                    ref_calib_.create(RO, E1, E2*S, CHA, N, 1, SLC);
                    traj_calib_.create(TRAJ, E1, E2*S, CHA, N, 1, SLC);
                }
            }
            else
            {
                if (combineN && combineS)
                {
                    ref_calib_.create(RO, E1*N*S, E2, CHA, 1, 1, SLC);
                    traj_calib_.create(TRAJ, E1*N*S, E2, CHA, 1, 1, SLC);
                }
                else if (combineN && !combineS)
                {
                    ref_calib_.create(RO, E1*N, E2, CHA, 1, S, SLC);
                    traj_calib_.create(TRAJ, E1*N, E2, CHA, 1, S, SLC);
                }
                else if (!combineN && combineS)
                {
                    ref_calib_.create(RO, E1*S, E2, CHA, N, 1, SLC);
                    traj_calib_.create(TRAJ, E1*S, E2, CHA, N, 1, SLC);
                }
            }

            const T* pRef = ref.begin();
            const float* pTraj= traj.begin();

            T* pRefCalib = ref_calib_.begin();
            float* pTrajCalib = traj_calib_.begin();

            long long num;
            long long NUM = SLC*S*N*CHA;

#pragma omp parallel for private(num) shared(RO, E1, E2, CHA, N, S, TRAJ, NUM, ref, traj)
            for (num = 0; num < NUM; num++)
            {
                size_t slc = num / (CHA*N*S);
                size_t s = (num - slc*CHA*N*S) / (CHA*N);
                size_t n = (num - slc*CHA*N*S - s*CHA*N) / (CHA);
                size_t cha = num - slc*CHA*N*S - s*CHA*N - n*CHA;

                size_t e1, e2;

                for (e2 = 0; e2 < E2; e2++)
                {
                    for (e1 = 0; e1 < E1; e1++)
                    {
                        const T* pR = &(ref(0, e1, e2, cha, n, s, slc));
                        const float* pT = &(traj(0, e1, e2, cha, n, s, slc));

                        T* pRC;
                        float* pTC;

                        if (combineN && combineS)
                        {
                            if (E2 > 1)
                            {
                                pRC = &(ref_calib_(0, e1, e2 + n*E2 + s*N*E2, cha, 0, 0, slc));
                                pTC = &(traj_calib_(0, e1, e2 + n*E2 + s*N*E2, cha, 0, 0, slc));
                            }
                            else
                            {
                                pRC = &(ref_calib_(0, e1 + n*E1 + s*N*E1, e2, cha, 0, 0, slc));
                                pTC = &(traj_calib_(0, e1 + n*E1 + s*N*E1, e2, cha, 0, 0, slc));
                            }
                        }
                        else if (combineN && !combineS)
                        {
                            if (E2 > 1)
                            {
                                pRC = &(ref_calib_(0, e1, e2 + n*E2, cha, 0, s, slc));
                                pTC = &(traj_calib_(0, e1, e2 + n*E2, cha, 0, s, slc));
                            }
                            else
                            {
                                pRC = &(ref_calib_(0, e1 + n*E1, e2, cha, 0, s, slc));
                                pTC = &(traj_calib_(0, e1 + n*E1, e2, cha, 0, s, slc));
                            }
                        }
                        else // !combineN && combineS
                        {
                            if (E2 > 1)
                            {
                                pRC = &(ref_calib_(0, e1, e2 + s*E2, cha, n, 0, slc));
                                pTC = &(traj_calib_(0, e1, e2 + s*E2, cha, n, 0, slc));
                            }
                            else
                            {
                                pRC = &(ref_calib_(0, e1 + s*E1, e2, cha, n, 0, slc));
                                pTC = &(traj_calib_(0, e1 + s*E1, e2, cha, n, 0, slc));
                            }
                        }

                        memcpy(pRC, pR, sizeof(T)*RO);
                        memcpy(pTC, pT, sizeof(float)*TRAJ);
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in referenceNonCartesianPreparer<T>::combine_along_N_S() ... ");
        }
    }

    template<typename T>
    void referenceNonCartesianPreparer<T>::prep_ref(const hoNDArray<T>& ref, const hoNDArray<float>& traj)
    {
        try
        {
            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N = ref.get_size(4);
            size_t S = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            GADGET_CHECK_THROW(traj.get_size(1) == E1);
            GADGET_CHECK_THROW(traj.get_size(2) == E2);
            GADGET_CHECK_THROW(traj.get_size(3) == CHA);
            GADGET_CHECK_THROW(traj.get_size(4) == N);
            GADGET_CHECK_THROW(traj.get_size(5) == S);
            GADGET_CHECK_THROW(traj.get_size(6) == SLC);

            if (calib_mode_ == ISMRMRD_interleaved)
            {
                if (interleave_dim_along_N_) // must combine along N
                {
                    if (combine_all_ref_S_)
                    {
                        this->combine_along_N_S(ref, traj, true, true);
                    }
                    else
                    {
                        this->combine_along_N_S(ref, traj, true, false);
                    }
                }
                else
                {
                    if (combine_all_ref_N_)
                    {
                        this->combine_along_N_S(ref, traj, true, true);
                    }
                    else
                    {
                        this->combine_along_N_S(ref, traj, false, true);
                    }
                }
            }
            else
            {
                this->combine_along_N_S(ref, traj, combine_all_ref_N_, combine_all_ref_S_);
            }

            if (!debug_folder_.empty())
            {
                // gt_exporter_.exportArrayComplex(ref_calib_, debug_folder_ + "ref_calib_after_combine_N_S");
            }

            // make the ref_coil_map_ and traj_coil_map_

            std::vector<size_t> dim;
            ref_calib_.get_dimensions(dim);
            ref_coil_map_.create(dim, ref_calib_.begin());

            traj_calib_.get_dimensions(dim);
            traj_coil_map_.create(dim, traj_calib_.begin());
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in referenceNonCartesianPreparer<T>::prep_ref() ... ");
        }
    }

    template<typename T>
    void referenceNonCartesianPreparer<T>::dump(std::ostream& os) const
    {
        using namespace std;
        os << "-------------------------------------------------------------------------------" << endl;
        os << "referenceNonCartesianPreparer ... " << endl;
        os << "calib_mode_ is " << Gadgetron::get_ismrmrd_calib_mode_name(calib_mode_) << endl;
        os << "combine_all_ref_N_ is " << combine_all_ref_N_ << endl;
        os << "combine_all_ref_S_ is " << combine_all_ref_S_ << endl;
        os << "interleave_dim_along_N_ is " << interleave_dim_along_N_ << endl;
        os << "verbose_ is " << verbose_ << endl;
        os << "-------------------------------------------------------------------------------" << endl;
    }

    template class EXPORTMRICORE referenceNonCartesianPreparer < std::complex<float> >;
    template class EXPORTMRICORE referenceNonCartesianPreparer < std::complex<double> >;
}
