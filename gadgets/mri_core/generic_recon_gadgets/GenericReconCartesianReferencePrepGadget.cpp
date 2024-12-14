
#include "GenericReconCartesianReferencePrepGadget.h"
#include <iomanip>

#include "hoNDArray_reductions.h"

namespace Gadgetron {

    GenericReconCartesianReferencePrepGadget::GenericReconCartesianReferencePrepGadget() : BaseClass()
    {
    }

    GenericReconCartesianReferencePrepGadget::~GenericReconCartesianReferencePrepGadget()
    {
    }

    int GenericReconCartesianReferencePrepGadget::process_config(const mrd::Header& header)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(header) == GADGET_OK, GADGET_FAIL);

        auto& h = header;

        if (!h.acquisition_system_information) {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        // -------------------------------------------------

        size_t NE = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        calib_mode_.resize(NE, mrd::CalibrationMode::kNoacceleration);
        ref_prepared_.resize(NE, false);

        for (size_t e = 0; e < h.encoding.size(); e++)
        {
            mrd::EncodingSpaceType e_space = h.encoding[e].encoded_space;
            mrd::EncodingSpaceType r_space = h.encoding[e].recon_space;
            mrd::EncodingLimitsType e_limits = h.encoding[e].encoding_limits;

            GDEBUG_CONDITION_STREAM(verbose.value(), "---> Encoding space : " << e << " <---");
            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding matrix size: " << e_space.matrix_size.x << " " << e_space.matrix_size.y << " " << e_space.matrix_size.z);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding field_of_view : " << e_space.field_of_view_mm.x << " " << e_space.field_of_view_mm.y << " " << e_space.field_of_view_mm.z);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Recon matrix size : " << r_space.matrix_size.x << " " << r_space.matrix_size.y << " " << r_space.matrix_size.z);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Recon field_of_view :  " << r_space.field_of_view_mm.x << " " << r_space.field_of_view_mm.y << " " << r_space.field_of_view_mm.z);

            calib_mode_[e] = mrd::CalibrationMode::kNoacceleration;

            if (!h.encoding[e].parallel_imaging)
            {
                GDEBUG_STREAM("Parallel Imaging section not found in header");
            }
            else
            {
                mrd::ParallelImagingType p_imaging = *h.encoding[0].parallel_imaging;
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE1 is " << p_imaging.acceleration_factor.kspace_encoding_step_1);
                GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE2 is " << p_imaging.acceleration_factor.kspace_encoding_step_2);

                if (p_imaging.acceleration_factor.kspace_encoding_step_1 > 1 || p_imaging.acceleration_factor.kspace_encoding_step_2 > 1)
                {
                    calib_mode_[e] = p_imaging.calibration_mode.value_or(mrd::CalibrationMode::kNoacceleration);
                }
            }
        }

        return GADGET_OK;
    }

    int GenericReconCartesianReferencePrepGadget::process(Gadgetron::GadgetContainerMessage< mrd::ReconData >* m1)
    {
        if (perform_timing.value()) { gt_timer_.start("GenericReconCartesianReferencePrepGadget::process"); }

        process_called_times_++;

        mrd::ReconData* recon_data = m1->getObjectPtr();
        if (recon_data->buffers.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_data->buffers.size() << " instead of " << num_encoding_spaces_);
        }

        GadgetContainerMessage<std::vector<mrd::WaveformUint32>>* wav =
            AsContainerMessage<std::vector<mrd::WaveformUint32>>(m1->cont());
        if (wav)
        {
            if (verbose.value())
            {
                GDEBUG_STREAM("Incoming recon_bit with " << wav->getObjectPtr()->size() << " wave form samples ");
            }
        }

        // a data buffer for N and S selection
        hoNDArray< std::complex<float> > ref_selected_N_S;

        // for every encoding space, prepare the recon_bit_->rbit_[e].ref_
        size_t e;
        for (e = 0; e < recon_data->buffers.size(); e++)
        {
            auto & rbit = recon_data->buffers[e];
            std::stringstream os;
            os << "_encoding_" << e;

            // -----------------------------------------
            // no acceleration mode
            // check the availability of ref data
            // -----------------------------------------
            if (prepare_ref_always.value() || !ref_prepared_[e])
            {
                // if no ref data is set, make copy the ref point from the  data
                if (!rbit.ref)
                {
                    rbit.ref = rbit.data;
                    ref_prepared_[e] = true;
                }
            }
            else
            {
                if (ref_prepared_[e])
                {
                    if (rbit.ref)
                    {
                        // remove the ref
                        rbit.ref = std::nullopt;
                    }
                }

                continue;
            }

            if ((!rbit.ref) || (rbit.ref->data.size()==0)) continue;

            if (!debug_folder_full_path_.empty())
            {
                this->gt_exporter_.export_array_complex(rbit.ref->data, debug_folder_full_path_ + "ref_data" + os.str());
            }

            // useful variables
            hoNDArray< std::complex<float> >& ref = (*rbit.ref).data;

            auto sampling_limits = (*rbit.ref).sampling.sampling_limits;

            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N = ref.get_size(4);
            size_t S = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            // -----------------------------------------
            // 1) average or pick the ref according to the input parameters;
            //    if interleaved mode, sampling times for every E1/E2 location is detected and line by line averaging is performed
            //    this is required when irregular cartesian sampling is used or number of frames cannot be divided in full by acceleration factor
            // 2) detect the sampled region and crop the ref data if needed
            // 3) update the sampling_limits
            // -----------------------------------------

            // if embedded mode, fill back ref if required
            if((calib_mode_[e] == mrd::CalibrationMode::kEmbedded) && ref_fill_into_data_embedded.value())
            {
                hoNDArray< std::complex<float> >& data = rbit.data.data;

                GADGET_CHECK_THROW(data.get_size(0) == RO);
                GADGET_CHECK_THROW(data.get_size(1) == E1);
                GADGET_CHECK_THROW(data.get_size(2) == E2);
                GADGET_CHECK_THROW(data.get_size(3) == CHA);
                GADGET_CHECK_THROW(data.get_size(6) == SLC);

                size_t slc, n, s, cha, e2, e1, ro;
                for (slc = 0; slc < SLC; slc++)
                {
                    for (s = 0; s < S; s++)
                    {
                        for (n = 0; n < N; n++)
                        {
                            for (e2 = 0; e2 < E2; e2++)
                            {
                                for (e1 = 0; e1 < E1; e1++)
                                {
                                    if (std::abs(ref(RO/2, e1, e2, 0, n, s, slc))>0 && std::abs(ref(RO/2, e1, e2, CHA-1, n, s, slc))>0)
                                    {
                                        for (cha = 0; cha < CHA; cha++)
                                        {
                                            std::complex<float>* pRef = &(ref(0, e1, e2, cha, n, s, slc));
                                            std::complex<float>* pData = &(data(0, e1, e2, cha, n, s, slc));
                                            memcpy(pData, pRef, sizeof(std::complex<float>)*RO);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                //if (!debug_folder_full_path_.empty()) { this->gt_exporter_.exportArrayComplex(data, debug_folder_full_path_ + "data_after_ref_filled_back" + os.str()); }
            }

            // stored the ref data ready for calibration
            hoNDArray< std::complex<float> > ref_calib;

            // step 1
            bool count_sampling_freq = (calib_mode_[e] == mrd::CalibrationMode::kInterleaved);

            bool valid_N_for_ref = (N_for_ref.value()<N && N_for_ref.value() >= 0);
            bool valid_S_for_ref = (S_for_ref.value()<S && S_for_ref.value() >= 0);

            if (!valid_N_for_ref && !valid_S_for_ref)
            {
                // use average N S
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::compute_averaged_data_N_S(ref, average_all_ref_N.value(), average_all_ref_S.value(), count_sampling_freq, ref_calib), GADGET_FAIL);
            }
            else if(valid_N_for_ref && !valid_S_for_ref)
            {
                // pick N, average S if needed
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::select_data_N_S(ref, true, N_for_ref.value(), false, 0, ref_selected_N_S), GADGET_FAIL);
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::compute_averaged_data_N_S(ref_selected_N_S, false, average_all_ref_S.value(), count_sampling_freq, ref_calib), GADGET_FAIL);
            }
            else if(!valid_N_for_ref && valid_S_for_ref)
            {
                // pick S, average N if needed
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::select_data_N_S(ref, false, 0, true, S_for_ref.value(), ref_selected_N_S), GADGET_FAIL);
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::compute_averaged_data_N_S(ref_selected_N_S, false, average_all_ref_S.value(), count_sampling_freq, ref_calib), GADGET_FAIL);
            }
            else if (valid_N_for_ref && valid_S_for_ref)
            {
                // pick N and S
                GADGET_CHECK_EXCEPTION_RETURN(Gadgetron::select_data_N_S(ref, true, N_for_ref.value(), true, S_for_ref.value(), ref_calib), GADGET_FAIL);
            }

            if (!debug_folder_full_path_.empty())
            {
                this->gt_exporter_.export_array_complex(ref_calib, debug_folder_full_path_ + "ref_calib" + os.str());
            }

            // step 2, detect sampled region in ref, along E1 and E2
            size_t start_E1(0), end_E1(0);
            auto t = Gadgetron::detect_sampled_region_E1(ref);
            start_E1 = std::get<0>(t);
            end_E1 = std::get<1>(t);

            size_t start_E2(0), end_E2(0);
            if (E2 > 1)
            {
                auto t = Gadgetron::detect_sampled_region_E2(ref);
                start_E2 = std::get<0>(t);
                end_E2 = std::get<1>(t);
            }

            // crop the ref_calib, along RO, E1 and E2 for separate or embedded mode
            vector_td<size_t, 3> crop_offset;
            crop_offset[0] = sampling_limits.kspace_encoding_step_0.minimum;
            crop_offset[1] = start_E1;
            crop_offset[2] = start_E2;

            vector_td<size_t, 3> crop_size;
            crop_size[0] = sampling_limits.kspace_encoding_step_0.maximum - sampling_limits.kspace_encoding_step_0.minimum + 1;
            crop_size[1] = end_E1 - start_E1 + 1;
            crop_size[2] = end_E2 - start_E2 + 1;

            if (crop_size[0]> (ref_calib.get_size(0) - crop_offset[0]))
            {
                crop_size[0] = ref_calib.get_size(0) - crop_offset[0];
            }

            if (crop_size[1]> (ref_calib.get_size(1) - crop_offset[1]))
            {
                crop_size[1] = ref_calib.get_size(1) - crop_offset[1];
            }

            if (crop_size[2]> (ref_calib.get_size(2) - crop_offset[2]))
            {
                crop_size[2] = ref_calib.get_size(2) - crop_offset[2];
            }

            hoNDArray< std::complex<float> > ref_recon_buf;
            Gadgetron::crop(crop_offset, crop_size, ref_calib, ref_recon_buf);
            ref_calib = ref_recon_buf;

            if (!debug_folder_full_path_.empty())
            {
                this->gt_exporter_.export_array_complex(ref_calib, debug_folder_full_path_ + "ref_calib_after_crop" + os.str());
            }

            // step 3, update the sampling limits
            sampling_limits.kspace_encoding_step_0.center = (uint16_t)(RO/2);

            sampling_limits.kspace_encoding_step_1.minimum = 0;
            sampling_limits.kspace_encoding_step_1.maximum = (uint16_t)(end_E1 - start_E1);

            sampling_limits.kspace_encoding_step_2.minimum = 0;
            sampling_limits.kspace_encoding_step_2.maximum = (uint16_t)(end_E2 - start_E2);

            if ( (calib_mode_[e] == mrd::CalibrationMode::kInterleaved) || (calib_mode_[e] == mrd::CalibrationMode::kNoacceleration) )
            {
                // need to keep the ref kspace center information
                sampling_limits.kspace_encoding_step_1.center = (uint16_t)(sampling_limits.kspace_encoding_step_1.center - start_E1);
                sampling_limits.kspace_encoding_step_2.center = (uint16_t)(sampling_limits.kspace_encoding_step_2.center - start_E2);
            }
            else
            {
                // separate, embedded mode, the ref center is the kspace center
                sampling_limits.kspace_encoding_step_1.center = (sampling_limits.kspace_encoding_step_1.maximum + 1) / 2;
                sampling_limits.kspace_encoding_step_2.center = (sampling_limits.kspace_encoding_step_2.maximum + 1) / 2;
            }

            if(sampling_limits.kspace_encoding_step_0.maximum>=RO)
            {
                sampling_limits.kspace_encoding_step_0.maximum = RO - 1;
            }

            ref = ref_calib;
            ref_prepared_[e] = true;

            (*rbit.ref).sampling.sampling_limits = sampling_limits;

            if (!debug_folder_full_path_.empty())
            {
                this->gt_exporter_.export_array_complex(rbit.ref->data, debug_folder_full_path_ + "ref_calib_final" + os.str());
            }
        }

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("Put ReconData to Q failed ... ");
            return GADGET_FAIL;
        }

        if (perform_timing.value()) { gt_timer_.stop(); }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GenericReconCartesianReferencePrepGadget)
}
