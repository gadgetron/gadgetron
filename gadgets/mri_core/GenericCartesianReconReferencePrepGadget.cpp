
#include "GenericCartesianReconReferencePrepGadget.h"
#include <iomanip>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron {

    GenericCartesianReconReferencePrepGadget::GenericCartesianReconReferencePrepGadget() : num_encoding_spaces_(1), process_called_times_(0)
    {
        gt_timer_.set_timing_in_destruction(false);
        gt_timer_local_.set_timing_in_destruction(false);
    }

    GenericCartesianReconReferencePrepGadget::~GenericCartesianReconReferencePrepGadget()
    {
    }

    int GenericCartesianReconReferencePrepGadget::process_config(ACE_Message_Block* mb)
    {
        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!h.acquisitionSystemInformation)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        // -------------------------------------------------

        size_t NE = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        calib_mode_.resize(NE, ISMRMRD_noacceleration);

        size_t e;
        for (e = 0; e < h.encoding.size(); e++)
        {
            ISMRMRD::EncodingSpace e_space = h.encoding[e].encodedSpace;
            ISMRMRD::EncodingSpace r_space = h.encoding[e].reconSpace;
            ISMRMRD::EncodingLimits e_limits = h.encoding[e].encodingLimits;

            GDEBUG_CONDITION_STREAM(verbose.value(), "---> Encoding space : " << e << " <---");
            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding matrix size: " << e_space.matrixSize.x << " " << e_space.matrixSize.y << " " << e_space.matrixSize.z);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Encoding field_of_view : " << e_space.fieldOfView_mm.x << " " << e_space.fieldOfView_mm.y << " " << e_space.fieldOfView_mm.z);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Recon matrix size : " << r_space.matrixSize.x << " " << r_space.matrixSize.y << " " << r_space.matrixSize.z);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Recon field_of_view :  " << r_space.fieldOfView_mm.x << " " << r_space.fieldOfView_mm.y << " " << r_space.fieldOfView_mm.z);

            if (!h.encoding[e].parallelImaging)
            {
                GDEBUG("Parallel Imaging section not found in header");
                return GADGET_FAIL;
            }

            ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;
            GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE1 is " << p_imaging.accelerationFactor.kspace_encoding_step_1);
            GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE2 is " << p_imaging.accelerationFactor.kspace_encoding_step_2);

            std::string calib = *p_imaging.calibrationMode;

            bool separate = (calib.compare("separate") == 0);
            bool embedded = (calib.compare("embedded") == 0);
            bool external = (calib.compare("external") == 0);
            bool interleaved = (calib.compare("interleaved") == 0);
            bool other = (calib.compare("other") == 0);

            calib_mode_[e] = Gadgetron::ISMRMRD_noacceleration;
            if (p_imaging.accelerationFactor.kspace_encoding_step_1 > 1 || p_imaging.accelerationFactor.kspace_encoding_step_2 > 1)
            {
                if (interleaved)
                    calib_mode_[e] = Gadgetron::ISMRMRD_interleaved;
                else if (embedded)
                    calib_mode_[e] = Gadgetron::ISMRMRD_embedded;
                else if (separate)
                    calib_mode_[e] = Gadgetron::ISMRMRD_separate;
                else if (external)
                    calib_mode_[e] = Gadgetron::ISMRMRD_external;
                else if (other)
                    calib_mode_[e] = Gadgetron::ISMRMRD_other;
            }
        }

        // ---------------------------------------------------------------------------------------------------------

        /*if (!debug_folder.value().empty())
        {
            Gadgetron::get_debug_folder_path(debug_folder.value(), debug_folder_full_path_);
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is " << debug_folder_full_path_);
        }
        else
        {
            GDEBUG_CONDITION_STREAM(verbose.value(), "Debug folder is not set ... ");
        }*/

        return GADGET_OK;
    }

    int GenericCartesianReconReferencePrepGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        process_called_times_++;

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        // for every encoding space, prepare the recon_bit_->rbit_[e].ref_
        size_t e;
        for (e = 0; e < recon_bit_->rbit_.size(); e++)
        {
            std::stringstream os;
            os << "_encoding_" << e;

            // -----------------------------------------
            // no acceleration mode
            // check the availability of ref data
            // -----------------------------------------
            if (calib_mode_[e] == Gadgetron::ISMRMRD_noacceleration)
            {
                // if no ref data is set, make the ref point to the  data
                if (recon_bit_->rbit_[e].ref_.data_.get_number_of_elements() == 0)
                {
                    std::vector<size_t> dim;
                    recon_bit_->rbit_[e].data_.data_.get_dimensions(dim);
                    recon_bit_->rbit_[e].ref_.data_.create(dim, recon_bit_->rbit_[e].data_.data_.begin());

                    if (recon_bit_->rbit_[e].data_.trajectory_.get_number_of_elements() > 0)
                    {
                        recon_bit_->rbit_[e].data_.trajectory_.get_dimensions(dim);
                        recon_bit_->rbit_[e].ref_.trajectory_.create(dim, recon_bit_->rbit_[e].data_.trajectory_.begin());
                    }

                    recon_bit_->rbit_[e].ref_.headers_ = recon_bit_->rbit_[e].data_.headers_;
                    recon_bit_->rbit_[e].ref_.sampling_ = recon_bit_->rbit_[e].data_.sampling_;
                }
            }

            if (recon_bit_->rbit_[e].ref_.data_.get_number_of_elements() == 0)
            {
                continue;
            }

            // useful variables
            hoNDArray< std::complex<float> >& ref = recon_bit_->rbit_[e].ref_.data_;

            SamplingLimit sampling_limits[3];
            sampling_limits[0] = recon_bit_->rbit_[e].ref_.sampling_.sampling_limits_[0];
            sampling_limits[1] = recon_bit_->rbit_[e].ref_.sampling_.sampling_limits_[1];
            sampling_limits[2] = recon_bit_->rbit_[e].ref_.sampling_.sampling_limits_[2];

            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N = ref.get_size(4);
            size_t S = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            // stored the ref data ready for calibration
            hoNDArray< std::complex<float> > ref_calib;

            /*if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.exportArrayComplex(recon_bit_->rbit_[e].ref_.data_, debug_folder_full_path_ + "ref" + os.str());
            }

            if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].ref_.trajectory_.get_number_of_elements() > 0)
            {
                gt_exporter_.exportArray(recon_bit_->rbit_[e].ref_.trajectory_, debug_folder_full_path_ + "ref_traj" + os.str());
            }*/

            // -----------------------------------------
            // 1) average the ref according to the input parameters; 
            //    if interleaved mode, sampling times for every E1/E2 location is detected and line by line averaging is performed 
            //    this is required when irregular cartesian sampling is used or number of frames cannot be divided in full by acceleration factor
            // 2) detect the sampled region and crop the ref data if needed
            // 3) update the sampling_limits
            // -----------------------------------------

            hoNDArray< std::complex<float> > ref_recon_buf;

            // step 1
            if (average_all_ref_N.value())
            {
                if (N > 1)
                {
                    if (calib_mode_[e] == ISMRMRD_interleaved)
                    {
                        hoNDArray<bool> sampled = Gadgetron::detect_readout_sampling_status(ref);
                        Gadgetron::sum_over_dimension(ref, ref_calib, 4);

                        // for every E1/E2 location, count how many times it is sampled for all N
                        size_t ro, e1, e2, cha, n, s, slc;
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (cha = 0; cha < CHA; cha++)
                            {
                                for (e2 = 0; e2 < E2; e2++)
                                {
                                    for (e1 = 0; e1 < E1; e1++)
                                    {
                                        float freq = 0;

                                        for (s = 0; s < S; s++)
                                        {
                                            for (n = 0; n < N; n++)
                                            {
                                                if (sampled(e1, e2, n, s, slc)) freq += 1;
                                            }

                                            if (freq > 1)
                                            {
                                                float freq_reciprocal = (float)(1.0 / freq);
                                                std::complex<float>* pAve = &(ref_calib(0, e1, e2, cha, 0, s, slc));
                                                for (ro = 0; ro < RO; ro++) pAve[ro] *= freq_reciprocal;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        Gadgetron::sum_over_dimension(ref, ref_calib, (size_t)4);
                        Gadgetron::scal((float)(1.0 / N), ref_calib);
                    }
                }
                else
                {
                    ref_calib = ref;
                }
            }
            else
            {
                ref_calib = ref;
            }

            if (average_all_ref_S.value())
            {
                if (S > 1)
                {
                    Gadgetron::sum_over_dimension(ref_calib, ref_recon_buf, 5);
                    Gadgetron::scal((float)(1.0 / S), ref_recon_buf);
                    ref_calib = ref_recon_buf;
                }
            }

            /*if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.exportArrayComplex(ref_calib, debug_folder_full_path_ + "ref_calib_after_averaging" + os.str());
            }*/

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

            // crop the ref_calib, along RO, E1 and E2
            vector_td<size_t, 3> crop_offset;
            crop_offset[0] = sampling_limits[0].min_;
            crop_offset[1] = start_E1;
            crop_offset[2] = start_E2;

            vector_td<size_t, 3> crop_size;
            crop_size[0] = sampling_limits[0].max_ - sampling_limits[0].min_ + 1;
            crop_size[1] = end_E1 - start_E1 + 1;
            crop_size[2] = end_E2 - start_E2 + 1;

            Gadgetron::crop(crop_offset, crop_size, &ref_calib, &ref_recon_buf);
            ref_calib = ref_recon_buf;

            /*if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.exportArrayComplex(ref_calib, debug_folder_full_path_ + "ref_calib_after_cropping" + os.str());
            }*/

            // step 3, update the sampling limits
            sampling_limits[0].center_ = (uint16_t)(RO/2);

            if ( (calib_mode_[e] == Gadgetron::ISMRMRD_interleaved) || (calib_mode_[e] == Gadgetron::ISMRMRD_noacceleration) )
            {
                // need to keep the ref kspace center information
                sampling_limits[1].min_ = (uint16_t)(start_E1);
                sampling_limits[1].max_ = (uint16_t)(end_E1);

                sampling_limits[2].min_ = (uint16_t)(start_E2);
                sampling_limits[2].max_ = (uint16_t)(end_E2);

                sampling_limits[1].center_ = (uint16_t)(E1 / 2);
                sampling_limits[2].center_ = (uint16_t)(E2 / 2);
            }
            else
            {
                // sepearate, embedded mode, the ref center is the kspace center
                sampling_limits[1].min_ = 0;
                sampling_limits[1].max_ = (uint16_t)(end_E1 - start_E1);

                sampling_limits[2].min_ = 0;
                sampling_limits[2].max_ = (uint16_t)(end_E2 - start_E2);

                sampling_limits[1].center_ = (sampling_limits[1].max_ + 1) / 2;
                sampling_limits[2].center_ = (sampling_limits[2].max_ + 1) / 2;
            }

            ref = ref_calib;
            recon_bit_->rbit_[e].ref_.sampling_.sampling_limits_[0] = sampling_limits[0];
            recon_bit_->rbit_[e].ref_.sampling_.sampling_limits_[1] = sampling_limits[1];
            recon_bit_->rbit_[e].ref_.sampling_.sampling_limits_[2] = sampling_limits[2];

            /*if (!debug_folder_full_path_.empty())
            {
                gt_exporter_.exportArrayComplex(recon_bit_->rbit_[e].ref_.data_, debug_folder_full_path_ + "ref_after_prep" + os.str());
            }*/
        }

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("Put IsmrmrdReconData to Q failed ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GenericCartesianReconReferencePrepGadget)
}
