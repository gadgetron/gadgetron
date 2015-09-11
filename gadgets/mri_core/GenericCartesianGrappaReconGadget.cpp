
#include "GenericCartesianGrappaReconGadget.h"
#include <iomanip>

#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "mri_core_coil_map_estimation.h"
#include "mri_core_grappa.h"
#include "mri_core_utility.h"

namespace Gadgetron {

    GenericCartesianGrappaReconGadget::GenericCartesianGrappaReconGadget() : num_encoding_spaces_(1), process_called_times_(0)
    {
    }

    GenericCartesianGrappaReconGadget::~GenericCartesianGrappaReconGadget()
    {
    }

    int GenericCartesianGrappaReconGadget::process_config(ACE_Message_Block* mb)
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

        meas_max_idx_.resize(NE);
        acceFactorE1_.resize(NE, 1);
        acceFactorE2_.resize(NE, 1);
        calib_mode_.resize(NE, ISMRMRD_noacceleration);

        recon_obj_.resize(NE);

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

            meas_max_idx_[e].kspace_encode_step_1 = (uint16_t)e_space.matrixSize.y - 1;
            meas_max_idx_[e].set = (e_limits.set && (e_limits.set->maximum > 0)) ? e_limits.set->maximum : 0;
            meas_max_idx_[e].phase = (e_limits.phase && (e_limits.phase->maximum > 0)) ? e_limits.phase->maximum : 0;

            meas_max_idx_[e].kspace_encode_step_2 = (uint16_t)e_space.matrixSize.z - 1;

            meas_max_idx_[e].contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum : 0;
            meas_max_idx_[e].slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;
            meas_max_idx_[e].repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;
            meas_max_idx_[e].slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;
            meas_max_idx_[e].average = e_limits.average ? e_limits.average->maximum : 0;
            meas_max_idx_[e].segment = 0;

            if (!h.encoding[e].parallelImaging)
            {
                GDEBUG("Parallel Imaging section not found in header");
                return GADGET_FAIL;
            }

            ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;

            acceFactorE1_[e] = p_imaging.accelerationFactor.kspace_encoding_step_1;
            acceFactorE2_[e] = p_imaging.accelerationFactor.kspace_encoding_step_2;
            GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE1 is " << acceFactorE1_[e]);
            GDEBUG_CONDITION_STREAM(verbose.value(), "acceFactorE2 is " << acceFactorE2_[e]);

            std::string calib = *p_imaging.calibrationMode;

            bool separate = (calib.compare("separate") == 0);
            bool embedded = (calib.compare("embedded") == 0);
            bool external = (calib.compare("external") == 0);
            bool interleaved = (calib.compare("interleaved") == 0);
            bool other = (calib.compare("other") == 0);

            calib_mode_[e] = Gadgetron::ISMRMRD_noacceleration;
            if (acceFactorE1_[e] > 1 || acceFactorE2_[e] > 1)
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

        return GADGET_OK;
    }

    int GenericCartesianGrappaReconGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        process_called_times_++;

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        // for every encoding space
        size_t e;
        for (e = 0; e < recon_bit_->rbit_.size(); e++)
        {
            this->ensure_recon_size(recon_bit_->rbit_[e], recon_obj_[e], e);

            if (calib_mode_[e] == Gadgetron::ISMRMRD_noacceleration && recon_bit_->rbit_[e].ref_.data_.get_number_of_elements()==0)
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

            if (recon_bit_->rbit_[e].ref_.data_.get_number_of_elements() > 0)
            {
                // after this step, the recon_obj_[e].ref_calib_ and recon_obj_[e].ref_coil_map_ are set
                this->prepare_ref(recon_bit_->rbit_[e], recon_obj_[e], e);

                // after this step, coil map is computed and stored in recon_obj_[e].coil_map_
                this->perform_coil_map_estimation(recon_bit_->rbit_[e], recon_obj_[e], e);

                // after this step, recon_obj_[e].kernel_, recon_obj_[e].kernelIm_, recon_obj_[e].unmixing_coeff_ are filled
                // gfactor is computed too
                this->perform_calib(recon_bit_->rbit_[e], recon_obj_[e], e);

                recon_bit_->rbit_[e].ref_.data_.clear();
                recon_bit_->rbit_[e].ref_.trajectory_.clear();
            }

            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0)
            {
                this->perform_unwrapping(recon_bit_->rbit_[e], recon_obj_[e], e);

                this->compute_image_header(recon_bit_->rbit_[e], recon_obj_[e], e);
                this->send_out_image_array(recon_bit_->rbit_[e], recon_obj_[e].recon_res_, e, image_series.value() + ((int)e + 1), GADGETRON_IMAGE_REGULAR);

                if (send_out_gfactor.value() && recon_obj_[e].gfactor_.get_number_of_elements()>0)
                {
                    IsmrmrdImageArray res;
                    Gadgetron::real_to_complex(recon_obj_[e].gfactor_, res.data_);
                    res.headers_ = recon_obj_[e].recon_res_.headers_;
                    res.meta_ = recon_obj_[e].recon_res_.meta_;

                    this->send_out_image_array(recon_bit_->rbit_[e], res, e, image_series.value() + 10 * ((int)e + 1), GADGETRON_IMAGE_GFACTOR);
                }
            }

            recon_obj_[e].recon_res_.data_.clear();
            recon_obj_[e].gfactor_.clear();
        }

        m1->release();
        return GADGET_OK;
    }

    size_t GenericCartesianGrappaReconGadget::compute_image_number(ISMRMRD::ImageHeader& imheader, size_t encoding, size_t CHA, size_t cha, size_t E2)
    {
        if (encoding >= meas_max_idx_.size())
        {
            GWARN_STREAM("encoding >= meas_max_idx_.size()");
            encoding = 0;
        }

        size_t SET = meas_max_idx_[encoding].set + 1;
        size_t REP = meas_max_idx_[encoding].repetition + 1;
        size_t PHS = meas_max_idx_[encoding].phase + 1;
        size_t SLC = meas_max_idx_[encoding].slice + 1;
        size_t CON = meas_max_idx_[encoding].contrast + 1;
        if (E2 == 0) E2 = 1;

        size_t imageNum = imheader.average*REP*SET*PHS*CON*SLC*E2*CHA + imheader.repetition*SET*PHS*CON*SLC*E2*CHA 
            + imheader.set*PHS*CON*SLC*E2*CHA + imheader.phase*CON*SLC*E2*CHA + imheader.contrast*SLC*E2*CHA + imheader.slice*E2*CHA + cha + 1;

        return imageNum;
    }

    int GenericCartesianGrappaReconGadget::send_out_image_array(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t encoding, int series_num, const std::string& data_role)
    {
        try
        {
            size_t RO = res.data_.get_size(0);
            size_t E1 = res.data_.get_size(1);
            size_t E2 = res.data_.get_size(2);
            size_t CHA = res.data_.get_size(3);
            size_t N = res.data_.get_size(4);
            size_t S = res.data_.get_size(5);
            size_t SLC = res.data_.get_size(6);

            GDEBUG_CONDITION_STREAM(true, "sending out image array, acquisition boundary [RO E1 E2 CHA N S SLC] = [" << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC << "] ");

            // compute image numbers and fill the image meta
            size_t n, s, slc;
            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        ISMRMRD::ImageHeader header = res.headers_(n, s, slc);

                        if (header.measurement_uid == 0) continue;

                        res.headers_(n, s, slc).image_index = (uint16_t)this->compute_image_number(res.headers_(n, s, slc), encoding, CHA, 0, E2);
                        res.headers_(n, s, slc).image_series_index = series_num;

                        size_t offset = n + s*N + slc*N*S;
                        res.meta_[offset].set(GADGETRON_IMAGENUMBER, (long)res.headers_(n, s, slc).image_index);
                        res.meta_[offset].set(GADGETRON_IMAGEPROCESSINGHISTORY, "GT");

                        if (data_role == GADGETRON_IMAGE_REGULAR)
                        {
                            res.headers_(n, s, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                            res.meta_[offset].append(GADGETRON_IMAGECOMMENT, "GT");

                            res.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, "_GT");
                            res.meta_[offset].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_REGULAR);
                        }
                        else if (data_role == GADGETRON_IMAGE_GFACTOR)
                        {
                            res.headers_(n, s, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                            res.meta_[offset].append(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_GFACTOR);
                            res.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_GFACTOR);
                            res.meta_[offset].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_GFACTOR);
                        }

                        if (verbose.value())
                        {
                            for (size_t cha = 0; cha < CHA; cha++)
                            {
                                GDEBUG_STREAM("sending out " << data_role << " image [CHA SLC CON PHS REP SET AVE] = [" << cha << " "<< res.headers_(n, s, slc).slice << " " << res.headers_(n, s, slc).contrast << " "<< res.headers_(n, s, slc).phase << " " << res.headers_(n, s, slc).repetition << " " << res.headers_(n, s, slc).set << " " << res.headers_(n, s, slc).average << " " << "] "<< " -- Image number -- " << res.headers_(n, s, slc).image_index); 
                            }
                        }
                    }
                }
            }

            // send out the images
            Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>* cm1 = new Gadgetron::GadgetContainerMessage<IsmrmrdImageArray>();
            *(cm1->getObjectPtr()) = res;

            if (this->next()->putq(cm1) < 0)
            {
                GERROR_STREAM("Put image array to Q failed ... ");
                return GADGET_FAIL;
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors in GenericCartesianGrappaReconGadget::send_out_image_array(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    void GenericCartesianGrappaReconGadget::ensure_recon_size(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding)
    {
        try
        {
            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);

            size_t N = recon_bit.data_.data_.get_size(4);
            size_t S = recon_bit.data_.data_.get_size(5);
            size_t SLC = recon_bit.data_.data_.get_size(6);

            // compensate for sampling limits under acceleration
            if (E1-1 - recon_bit.data_.sampling_.sampling_limits_[1].max_ < acceFactorE1_[encoding])
            {
                recon_bit.data_.sampling_.sampling_limits_[1].max_ = (uint16_t)E1 - 1;
            }

            if ( (E2>1) && (E2-1 - recon_bit.data_.sampling_.sampling_limits_[2].max_ < acceFactorE2_[encoding]))
            {
                recon_bit.data_.sampling_.sampling_limits_[2].max_ = (uint16_t)E2 - 1;
            }

            float spacingE1 = recon_bit.data_.sampling_.recon_FOV_[1] / recon_bit.data_.sampling_.recon_matrix_[1];
            size_t encodingE1 = (size_t)std::floor(recon_bit.data_.sampling_.encoded_FOV_[1] / spacingE1 + 0.5);

            if (encodingE1 > E1)
            {
                GDEBUG_STREAM("recon_squared_pixel is true; change encoding E1 to be " << encodingE1);

                // pad the data
                hoNDArray< std::complex<float> > dataPadded;
                Gadgetron::pad(RO, encodingE1, &recon_bit.data_.data_, &dataPadded);
                recon_bit.data_.data_ = dataPadded;

                // update the sampling_limits
                uint16_t offsetE1 = (uint16_t)(encodingE1 / 2 - E1 / 2);

                recon_bit.data_.sampling_.sampling_limits_[1].min_ += offsetE1;
                recon_bit.data_.sampling_.sampling_limits_[1].max_ += offsetE1;

                // update image headers
                size_t headerE1 = recon_bit.data_.headers_.get_size(0);
                size_t headerE2 = recon_bit.data_.headers_.get_size(1);

                size_t n, s, slc, e1, e2;

                for (slc = 0; slc < SLC; slc++)
                {
                    for (s = 0; s < S; s++)
                    {
                        for (n = 0; n < N; n++)
                        {
                            for (e2 = 0; e2 < headerE2; e2++)
                            {
                                for (e1 = 0; e1 < headerE1; e1++)
                                {
                                    if (recon_bit.data_.headers_(e1, e2, n, s, slc).measurement_uid != 0)
                                    {
                                        recon_bit.data_.headers_(e1, e2, n, s, slc).idx.kspace_encode_step_1 += offsetE1;
                                    }
                                }
                            }
                        }
                    }
                }

                // if calib_mode is embedded, pad the ref
                if (calib_mode_[encoding] == Gadgetron::ISMRMRD_embedded)
                {
                    if (recon_bit.ref_.data_.get_size(0) == RO && recon_bit.ref_.data_.get_size(1) == E1)
                    {
                        Gadgetron::pad(RO, encodingE1, &recon_bit.ref_.data_, &dataPadded);
                        recon_bit.ref_.data_ = dataPadded;

                        recon_bit.ref_.sampling_.sampling_limits_[1].min_ += offsetE1;
                        recon_bit.ref_.sampling_.sampling_limits_[1].max_ += offsetE1;
                    }
                }
            }
            else
            {
                GDEBUG_STREAM("it is not required to change encoding E1 ... ");
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericCartesianGrappaReconGadget::ensure_recon_size(...) ... ");
        }
    }

    void GenericCartesianGrappaReconGadget::prepare_ref(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t encoding)
    {
        try
        {
            hoNDArray< std::complex<float> >& ref = recon_bit.ref_.data_;
            hoNDArray< std::complex<float> >& ref_calib = recon_obj.ref_calib_;
            hoNDArray< std::complex<float> >& ref_coil_map = recon_obj.ref_coil_map_;

            // sampling limits
            SamplingLimit sampling_limits[3];
            sampling_limits[0] = recon_bit.ref_.sampling_.sampling_limits_[0];
            sampling_limits[1] = recon_bit.ref_.sampling_.sampling_limits_[1];
            sampling_limits[2] = recon_bit.ref_.sampling_.sampling_limits_[2];

            // recon size
            size_t recon_RO = recon_bit.data_.data_.get_size(0);
            size_t recon_E1 = recon_bit.data_.data_.get_size(1);
            size_t recon_E2 = recon_bit.data_.data_.get_size(2);

            // filter ref coil map
            bool filter_ref_coil_map = true;

            // ref array size
            size_t RO = ref.get_size(0);
            size_t E1 = ref.get_size(1);
            size_t E2 = ref.get_size(2);
            size_t CHA = ref.get_size(3);
            size_t N = ref.get_size(4);
            size_t S = ref.get_size(5);
            size_t SLC = ref.get_size(6);

            if (calib_mode_[encoding] == ISMRMRD_separate || calib_mode_[encoding] == ISMRMRD_external || calib_mode_[encoding] == ISMRMRD_other)
            {
                if (average_all_ref_N.value())
                {
                    if (N > 1)
                    {
                        Gadgetron::sum_over_dimension(ref, ref_calib, (size_t)4);
                        Gadgetron::scal((float)(1.0 / N), ref_calib);
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
                        hoNDArray< std::complex<float> > ref_recon_buf;
                        Gadgetron::sum_over_dimension(ref_calib, ref_recon_buf, 5);
                        Gadgetron::scal((float)(1.0 / S), ref_recon_buf);
                        ref_calib = ref_recon_buf;
                    }
                }

                hoNDArray< std::complex<float> > ref_recon_buf;

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
                crop_offset[0] = sampling_limits[0].min_;
                crop_offset[1] = start_E1;
                crop_offset[2] = start_E2;

                vector_td<size_t, 3> crop_size;
                crop_size[0] = sampling_limits[0].max_ - sampling_limits[0].min_ + 1;
                crop_size[1] = end_E1 - start_E1 + 1;
                crop_size[2] = end_E2 - start_E2 + 1;

                Gadgetron::crop(crop_offset, crop_size, &ref_calib, &ref_recon_buf);
                ref_calib = ref_recon_buf;

                // create filter if needed
                if (filter_ref_coil_map)
                {
                    if (filter_RO_ref_coi_map_.get_size(0) != RO
                        || filter_E1_ref_coi_map_.get_size(0) != ref_calib.get_size(1)
                        || ((E2 > 1) && (filter_E2_ref_coi_map_.get_size(0) != ref_calib.get_size(2))))
                    {
                        SamplingLimit sE1;
                        sE1.min_ = 0;
                        sE1.max_ = (uint16_t)ref_calib.get_size(1) - 1;

                        SamplingLimit sE2;
                        sE2.min_ = 0;
                        sE2.max_ = (uint16_t)ref_calib.get_size(2) - 1;

                        Gadgetron::generate_ref_filter_for_coil_map(ref_calib, sampling_limits[0], sE1, sE2, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_);
                    }
                }

                // filter the ref_coil_map_
                ref_coil_map = ref_calib;

                if (filter_ref_coil_map)
                {
                    if (E2 > 1)
                    {
                        Gadgetron::apply_kspace_filter_ROE1E2(ref_coil_map, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_, ref_recon_buf);
                    }
                    else
                    {
                        Gadgetron::apply_kspace_filter_ROE1(ref_coil_map, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, ref_recon_buf);
                    }

                    ref_coil_map = ref_recon_buf;
                }

                // pad the ref_coil_map into the data array
                Gadgetron::pad(recon_RO, recon_E1, recon_E2, &ref_coil_map, &ref_recon_buf);
                ref_coil_map = ref_recon_buf;
            }
            else
            {
                GADGET_THROW("To be implemented ... ");
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericCartesianGrappaReconGadget::prepare_ref(...) ... ");
        }
    }

    void GenericCartesianGrappaReconGadget::perform_coil_map_estimation(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            recon_obj.coil_map_ = recon_obj.ref_coil_map_;
            Gadgetron::clear(recon_obj.coil_map_);

            size_t E2 = recon_obj.ref_coil_map_.get_size(2);
            if (E2 > 1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_obj.ref_coil_map_, complex_im_recon_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_obj.ref_coil_map_, complex_im_recon_buf_);
            }

            this->compute_coil_map(complex_im_recon_buf_, recon_obj.coil_map_);
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericCartesianGrappaReconGadget::perform_coil_map_estimation(...) ... ");
        }
    }

    void GenericCartesianGrappaReconGadget::perform_calib(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            typedef  std::complex<float> T;

            if (acceFactorE1_[e] <= 1 && acceFactorE2_[e] <= 1) return;

            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);

            hoNDArray<T>& src = recon_obj.ref_calib_;
            hoNDArray<T>& dst = recon_obj.ref_calib_;

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);
            size_t ref_N = src.get_size(4);
            size_t ref_S = src.get_size(5);
            size_t ref_SLC = src.get_size(6);

            size_t dstCHA = dst.get_size(3);

            // allocate buffer for kernels
            size_t kRO = grappa_kSize_RO.value();
            size_t kNE1 = grappa_kSize_E1.value();
            size_t kNE2 = grappa_kSize_E2.value();

            size_t convKRO, convKE1, convKE2;

            if (E2 > 1)
            {
                std::vector<int> kE1, oE1;
                std::vector<int> kE2, oE2;
                bool fitItself = true;

                grappa3d_kerPattern(kE1, oE1, kE2, oE2, convKRO, convKE1, convKE2, (size_t)acceFactorE1_[e], (size_t)acceFactorE2_[e], kRO, kNE1, kNE2, fitItself);

                recon_obj.kernel_.create(convKRO, convKE1, convKE2, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
                recon_obj.unmixing_coeff_.create(RO, E1, E2, srcCHA, ref_N, ref_S, ref_SLC);
                recon_obj.gfactor_.create(RO, E1, E2, 1, ref_N, ref_S, ref_SLC);
            }
            else
            {
                std::vector<int> kE1, oE1;
                bool fitItself = true;

                Gadgetron::grappa2d_kerPattern(kE1, oE1, convKRO, convKE1, (size_t)acceFactorE1_[e], kRO, kNE1, fitItself);

                recon_obj.kernel_.create(convKRO, convKE1, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
                recon_obj.kernelIm_.create(RO, E1, 1, srcCHA, dstCHA, ref_N, ref_S, ref_SLC);
                recon_obj.unmixing_coeff_.create(RO, E1, 1, srcCHA, ref_N, ref_S, ref_SLC);
                recon_obj.gfactor_.create(RO, E1, 1, 1, ref_N, ref_S, ref_SLC);
            }

            Gadgetron::clear(recon_obj.kernel_);
            Gadgetron::clear(recon_obj.kernelIm_);
            Gadgetron::clear(recon_obj.unmixing_coeff_);
            Gadgetron::clear(recon_obj.gfactor_);

            // estimate kernel

            long long num = ref_N*ref_S*ref_SLC;

            long long ii;

#pragma omp parallel for default(none) private(ii) shared(src, dst, recon_obj, e, num, ref_N, ref_S, ref_RO, ref_E1, ref_E2, RO, E1, E2, dstCHA, srcCHA)
            for (ii = 0; ii < num; ii++)
            {
                size_t slc = ii / (ref_N*ref_S);
                size_t s = (ii - slc*ref_N*ref_S) / (ref_N);
                size_t n = ii - slc*ref_N*ref_S - s*ref_N;

                T* pSrc = &(src(0, 0, 0, 0, n, s, slc));
                hoNDArray<T> ref_src(ref_RO, ref_E1, ref_E2, srcCHA, pSrc);

                T* pDst = &(dst(0, 0, 0, 0, n, s, slc));
                hoNDArray<T> ref_dst(ref_RO, ref_E1, ref_E2, dstCHA, pDst);

                this->perform_calib_impl(n, s, slc, e, ref_src, ref_dst, recon_obj);
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericCartesianGrappaReconGadget::perform_calib(...) ... ");
        }
    }

    void GenericCartesianGrappaReconGadget::perform_calib_impl(size_t n, size_t s, size_t slc, size_t e, const hoNDArray< std::complex<float> >& src, const hoNDArray< std::complex<float> >& dst, ReconObjType& recon_obj)
    {
        try
        {
            typedef std::complex<float> T;

            size_t RO = recon_obj.ref_coil_map_.get_size(0);
            size_t E1 = recon_obj.ref_coil_map_.get_size(1);
            size_t E2 = recon_obj.ref_coil_map_.get_size(2);

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);

            size_t dstCHA = dst.get_size(3);

            // allocate buffer for kernels
            size_t kRO = grappa_kSize_RO.value();
            size_t kNE1 = grappa_kSize_E1.value();
            size_t kNE2 = grappa_kSize_E2.value();

            size_t convKRO = recon_obj.kernel_.get_size(0);
            size_t convKNE1 = recon_obj.kernel_.get_size(1);

            if (E2 > 1)
            {
                size_t convKNE2 = recon_obj.kernel_.get_size(2);

                hoNDArray<T> ker(convKRO, convKNE1, convKNE2, srcCHA, dstCHA, &(recon_obj.kernel_(0, 0, 0, 0, 0, n, s, slc)));
                hoNDArray<T> coilMap(RO, E1, E2, dstCHA, &(recon_obj.coil_map_(0, 0, 0, 0, n, s, slc)));

                hoNDArray<T> unmixC(RO, E1, E2, srcCHA, &(recon_obj.unmixing_coeff_(0, 0, 0, 0, n, s, slc)));
                hoNDArray<float> gFactor(RO, E1, E2, 1, &(recon_obj.gfactor_(0, 0, 0, 0, n, s, slc)));

                Gadgetron::grappa3d_calib_convolution_kernel(src, dst, (size_t)acceFactorE1_[e], (size_t)acceFactorE2_[e], grappa_reg_lamda.value(), grappa_calib_over_determine_ratio.value(), kRO, kNE1, kNE2, ker);
                Gadgetron::grappa3d_unmixing_coeff(ker, coilMap, (size_t)acceFactorE1_[e], (size_t)acceFactorE2_[e], unmixC, gFactor);
            }
            else
            {
                hoNDArray<T> acsSrc(ref_RO, ref_E1, srcCHA, const_cast<T*>(src.begin()));
                hoNDArray<T> acsDst(ref_RO, ref_E1, dstCHA, const_cast<T*>(dst.begin()));

                hoNDArray<T> convKer(convKRO, convKNE1, srcCHA, dstCHA, &(recon_obj.kernel_(0, 0, 0, 0, n, s, slc)));
                hoNDArray<T> kIm(RO, E1, srcCHA, dstCHA, &(recon_obj.kernelIm_(0, 0, 0, 0, 0, n, s, slc)));
                hoNDArray<T> coilMap(RO, E1, dstCHA, &(recon_obj.coil_map_(0, 0, 0, 0, n, s, slc)));
                hoNDArray<T> unmixC(RO, E1, srcCHA, &(recon_obj.unmixing_coeff_(0, 0, 0, 0, n, s, slc)));

                Gadgetron::grappa2d_calib_convolution_kernel(acsSrc, acsDst, (size_t)acceFactorE1_[e], grappa_reg_lamda.value(), kRO, kNE1, convKer);
                Gadgetron::grappa2d_image_domain_kernel(convKer, RO, E1, kIm);

                hoNDArray<float> gFactor;
                Gadgetron::grappa2d_unmixing_coeff(kIm, coilMap, (size_t)acceFactorE1_[e], unmixC, gFactor);

                memcpy(&(recon_obj.gfactor_(0, 0, 0, 0, n, s, slc)), gFactor.begin(), gFactor.get_number_of_bytes());
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericCartesianGrappaReconGadget::perform_calib_impl(n, s, slc) ... ");
        }
    }

    void GenericCartesianGrappaReconGadget::perform_unwrapping(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            typedef std::complex<float> T;

            size_t RO = recon_bit.data_.data_.get_size(0);
            size_t E1 = recon_bit.data_.data_.get_size(1);
            size_t E2 = recon_bit.data_.data_.get_size(2);
            size_t dstCHA = recon_bit.data_.data_.get_size(3);
            size_t N = recon_bit.data_.data_.get_size(4);
            size_t S = recon_bit.data_.data_.get_size(5);
            size_t SLC = recon_bit.data_.data_.get_size(6);

            hoNDArray<T>& src = recon_obj.ref_calib_;
            hoNDArray<T>& dst = recon_obj.ref_calib_;

            size_t ref_RO = src.get_size(0);
            size_t ref_E1 = src.get_size(1);
            size_t ref_E2 = src.get_size(2);
            size_t srcCHA = src.get_size(3);
            size_t ref_N = src.get_size(4);
            size_t ref_S = src.get_size(5);
            size_t ref_SLC = src.get_size(6);

            size_t convkRO = recon_obj.kernel_.get_size(0);
            size_t convkE1 = recon_obj.kernel_.get_size(1);
            size_t convkE2 = recon_obj.kernel_.get_size(2);

            recon_obj.recon_res_.data_.create(RO, E1, E2, 1, N, S, SLC);

            // compute aliased images
            data_recon_buf_.create(RO, E1, E2, dstCHA, N, S, SLC);

            if (E2>1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(recon_bit.data_.data_, complex_im_recon_buf_, data_recon_buf_);
            }

            // SNR unit scaling
            float effectiveAcceFactor = acceFactorE1_[e] * acceFactorE2_[e];
            if (effectiveAcceFactor > 1)
            {
                float fftCompensationRatio = (float)(1.0 / std::sqrt(effectiveAcceFactor));
                Gadgetron::scal(fftCompensationRatio, complex_im_recon_buf_);
            }

            // unwrapping

            long long num = N*S*SLC;

            long long ii;

            if (effectiveAcceFactor > 1)
            {
                if (E2 > 1)
                {
#pragma omp parallel default(none) private(ii) shared(num, N, S, RO, E1, E2, srcCHA, convkRO, convkE1, convkE2, ref_N, ref_S, recon_obj, dstCHA, e)
                    {
#pragma omp for 
                        for (ii = 0; ii < num; ii++)
                        {
                            size_t slc = ii / (N*S);
                            size_t s = (ii - slc*N*S) / N;
                            size_t n = ii - slc*N*S - s*N;

                            // combined channels
                            T* pIm = &(complex_im_recon_buf_(0, 0, 0, 0, n, s, slc));
                            hoNDArray<T> aliasedIm(RO, E1, E2, srcCHA, 1, pIm);

                            size_t usedN = n;
                            if (n >= ref_N) usedN = ref_N - 1;

                            size_t usedS = s;
                            if (s >= ref_S) usedS = ref_S - 1;

                            T* pKer = &(recon_obj.kernel_(0, 0, 0, 0, 0, usedN, usedS, slc));
                            hoNDArray<T> ker(convkRO, convkE1, convkE2, srcCHA, dstCHA, pKer);

                            T* pUnmix = &(recon_obj.unmixing_coeff_(0, 0, 0, 0, usedN, usedS, slc));
                            hoNDArray<T> unmixing(RO, E1, E2, srcCHA, pUnmix);

                            T* pRes = &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc));
                            hoNDArray<T> res(RO, E1, E2, 1, pRes);

                            Gadgetron::apply_unmix_coeff_aliased_image_3D(aliasedIm, unmixing, res);
                        }
                    }
                }
                else
                {
#pragma omp parallel default(none) private(ii) shared(num, N, S, RO, E1, E2, srcCHA, ref_N, ref_S, recon_obj, dstCHA, e)
                    {
#pragma omp for 
                        for (ii = 0; ii < num; ii++)
                        {
                            size_t slc = ii / (N*S);
                            size_t s = (ii - slc*N*S) / N;
                            size_t n = ii - slc*N*S - s*N;

                            // combined channels
                            T* pIm = &(complex_im_recon_buf_(0, 0, 0, 0, n, s, slc));
                            hoNDArray<T> aliasedIm(RO, E1, srcCHA, pIm);

                            size_t usedN = n;
                            if (n >= ref_N) usedN = ref_N - 1;

                            size_t usedS = s;
                            if (s >= ref_S) usedS = ref_S - 1;

                            T* pUnmix = &(recon_obj.unmixing_coeff_(0, 0, 0, 0, usedN, usedS, slc));
                            hoNDArray<T> unmixing(RO, E1, srcCHA, pUnmix);

                            T* pRes = &(recon_obj.recon_res_.data_(0, 0, 0, 0, n, s, slc));
                            hoNDArray<T> res(RO, E1, 1, pRes);

                            Gadgetron::apply_unmix_coeff_aliased_image(aliasedIm, unmixing, res);
                        }
                    }
                }
            }
            else
            {
                this->coil_combination(complex_im_recon_buf_, recon_obj.coil_map_, recon_obj.recon_res_.data_);
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericCartesianGrappaReconGadget::perform_unwrapping(...) ... ");
        }
    }

    void GenericCartesianGrappaReconGadget::compute_coil_map(const hoNDArray< std::complex<float> >& complexIm, hoNDArray< std::complex<float> >& coilMap)
    {
        try
        {
            typedef  std::complex<float> T;

            size_t RO = complexIm.get_size(0);
            size_t E1 = complexIm.get_size(1);
            size_t E2 = complexIm.get_size(2);
            size_t CHA = complexIm.get_size(3);

            if (!complexIm.dimensions_equal(&coilMap))
            {
                coilMap = complexIm;
            }

            if (CHA <= 1)
            {
                GWARN_STREAM("coilMapMakerInati<T>::make_coil_map, CHA <= 1");
                return;
            }

            size_t ks = 7;
            size_t power = 3;
            size_t num = complexIm.get_number_of_elements() / (RO*E1*E2*CHA);

            long long n;

            if (E2 > 1)
            {
                for (n = 0; n < (long long)num; n++)
                {

                    hoNDArray<T> im(RO, E1, E2, CHA, const_cast<T*>(complexIm.begin() + n*RO*E1*E2*CHA));
                    hoNDArray<T> cmap(RO, E1, E2, CHA, coilMap.begin() + n*RO*E1*E2*CHA);

                    Gadgetron::coil_map_3d_Inati(im, cmap, ks, power);
                }
            }
            else
            {
#pragma omp parallel for default(none) private(n) shared(num, RO, E1, CHA, complexIm, coilMap, ks, power) if(num>8)
                for (n = 0; n < (long long)num; n++)
                {
                    hoNDArray<T> im(RO, E1, CHA, const_cast<T*>(complexIm.begin()) + n*RO*E1*CHA);
                    hoNDArray<T> cmap(RO, E1, CHA, coilMap.begin() + n*RO*E1*CHA);

                    Gadgetron::coil_map_2d_Inati(im, cmap, ks, power);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericCartesianGrappaReconGadget::compute_coil_map(...) ... ")
        }
    }

    void GenericCartesianGrappaReconGadget::coil_combination(const hoNDArray< std::complex<float> >& complex_im, const hoNDArray< std::complex<float> >& coil_map, hoNDArray< std::complex<float> >& res)
    {
        try
        {
            size_t RO = complex_im.get_size(0);
            size_t E1 = complex_im.get_size(1);
            size_t E2 = complex_im.get_size(2);
            size_t CHA = complex_im.get_size(3);
            size_t N = complex_im.get_size(4);
            size_t S = complex_im.get_size(5);
            size_t SLC = complex_im.get_size(6);

            size_t dstCHA = 1;

            size_t coilCHA = coil_map.get_size(3);

            size_t coil_N = coil_map.get_size(4);
            size_t coil_S = coil_map.get_size(5);

            typedef std::complex<float> T;

            T* pIm = const_cast<T*>(complex_im.begin());
            T* pCoilMap = const_cast<T*>(coil_map.begin());

            res.create(RO, E1, E2, dstCHA, N, S, SLC);

            long long num = SLC*S*N;

            long long ii;

#pragma omp parallel default(none) private(ii) shared(RO, E1, E2, CHA, N, S, SLC, coilCHA, dstCHA, num, coil_N, coil_S, pIm, pCoilMap, res) if(num>8)
            {
                hoNDArray<T> dataBuf;
                dataBuf.create(RO, E1, E2, coilCHA);

                hoNDArray<T> dataBufCombined;
                dataBufCombined.create(RO, E1, E2, 1);

#pragma omp for 
                for (ii = 0; ii < num; ii++)
                {
                    size_t slc = ii / (N*S);
                    size_t s = (ii - slc*N*S) / N;
                    size_t n = ii - slc*N*S - s*N;

                    size_t ns = s;
                    if (ns >= coil_S) ns = coil_S - 1;

                    size_t nn = n;
                    if (nn >= coil_N) nn = coil_N - 1;

                    hoNDArray<T> imCurr(RO, E1, E2, CHA, pIm + slc*RO*E1*E2*CHA*N*S + s*RO*E1*E2*CHA*N + n*RO*E1*E2*CHA);
                    hoNDArray<T> coilMapCurr(RO, E1, E2, CHA, pCoilMap + slc*RO*E1*E2*CHA*coil_N*coil_S + ns*RO*E1*E2*CHA*coil_N + nn*RO*E1*E2*CHA);

                    Gadgetron::multiplyConj(imCurr, coilMapCurr, dataBuf);
                    Gadgetron::sum_over_dimension(dataBuf, dataBufCombined, 3);

                    memcpy(res.begin() + slc*RO*E1*E2*dstCHA*N*S + s*RO*E1*E2*dstCHA*N + n*RO*E1*E2*dstCHA, dataBufCombined.begin(), sizeof(T)*RO*E1*E2);
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericCartesianGrappaReconGadget::coil_combination(complex_im) ... ");
        }
    }

    void GenericCartesianGrappaReconGadget::compute_image_header(IsmrmrdReconBit& recon_bit, ReconObjType& recon_obj, size_t e)
    {
        try
        {
            size_t RO = recon_obj.recon_res_.data_.get_size(0);
            size_t E1 = recon_obj.recon_res_.data_.get_size(1);
            size_t E2 = recon_obj.recon_res_.data_.get_size(2);
            size_t CHA = recon_obj.recon_res_.data_.get_size(3);
            size_t N = recon_obj.recon_res_.data_.get_size(4);
            size_t S = recon_obj.recon_res_.data_.get_size(5);
            size_t SLC = recon_obj.recon_res_.data_.get_size(6);

            GADGET_CHECK_THROW(N == recon_bit.data_.headers_.get_size(2));
            GADGET_CHECK_THROW(S == recon_bit.data_.headers_.get_size(3));

            recon_obj.recon_res_.headers_.create(N, S, SLC);
            recon_obj.recon_res_.meta_.resize(N*S*SLC);

            size_t n, s, slc;

            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        size_t header_E1 = recon_bit.data_.headers_.get_size(0);
                        size_t header_E2 = recon_bit.data_.headers_.get_size(1);

                        // for every kspace, find the recorded header which is closest to the kspace center [E1/2 E2/2]
                        ISMRMRD::AcquisitionHeader acq_header;

                        long long bestE1 = E1 + 1;
                        long long bestE2 = E2 + 1;

                        size_t e1, e2;
                        for (e2 = 0; e2 < header_E2; e2++)
                        {
                            for (e1 = 0; e1 < header_E1; e1++)
                            {
                                ISMRMRD::AcquisitionHeader& curr_header = recon_bit.data_.headers_(e1, e2, n, s, slc);

                                if (curr_header.measurement_uid != 0) // a valid header
                                {
                                    if (E2 > 1)
                                    {
                                        if (std::abs((long long)curr_header.idx.kspace_encode_step_1 - (long long)(E1 / 2)) < bestE1
                                            && std::abs((long long)curr_header.idx.kspace_encode_step_2 - (long long)(E2 / 2)) < bestE2)
                                        {
                                            bestE1 = std::abs((long long)curr_header.idx.kspace_encode_step_1 - (long long)E1 / 2);
                                            bestE2 = std::abs((long long)curr_header.idx.kspace_encode_step_2 - (long long)E2 / 2);

                                            acq_header = curr_header;
                                        }
                                    }
                                    else
                                    {
                                        if (std::abs((long long)curr_header.idx.kspace_encode_step_1 - (long long)(E1 / 2)) < bestE1)
                                        {
                                            bestE1 = std::abs((long long)curr_header.idx.kspace_encode_step_1 - (long long)E1 / 2);

                                            acq_header = curr_header;
                                        }
                                    }
                                }
                            }
                        }

                        if (acq_header.measurement_uid == 0)
                        {
                            std::ostringstream ostr;
                            ostr << "Cannot create valid image header : n = " << n << ", s = " << s << ", slc = " << slc;
                            GADGET_THROW(ostr.str());
                        }
                        else
                        {
                            ISMRMRD::ImageHeader& im_header = recon_obj.recon_res_.headers_(n, s, slc);
                            ISMRMRD::MetaContainer& meta = recon_obj.recon_res_.meta_[n + s*N + slc*N*S];

                            im_header.version = acq_header.version;
                            im_header.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                            im_header.flags = acq_header.flags;
                            im_header.measurement_uid = acq_header.measurement_uid;

                            im_header.matrix_size[0] = (uint16_t)RO;
                            im_header.matrix_size[1] = (uint16_t)E1;
                            im_header.matrix_size[2] = (uint16_t)E2;

                            im_header.field_of_view[0] = recon_bit.data_.sampling_.recon_FOV_[0];
                            im_header.field_of_view[1] = recon_bit.data_.sampling_.recon_FOV_[1];
                            im_header.field_of_view[2] = recon_bit.data_.sampling_.recon_FOV_[2];

                            im_header.channels = (uint16_t)CHA;

                            im_header.position[0] = acq_header.position[0];
                            im_header.position[1] = acq_header.position[1];
                            im_header.position[2] = acq_header.position[2];

                            im_header.read_dir[0] = acq_header.read_dir[0];
                            im_header.read_dir[1] = acq_header.read_dir[1];
                            im_header.read_dir[2] = acq_header.read_dir[2];

                            im_header.phase_dir[0] = acq_header.phase_dir[0];
                            im_header.phase_dir[1] = acq_header.phase_dir[1];
                            im_header.phase_dir[2] = acq_header.phase_dir[2];

                            im_header.slice_dir[0] = acq_header.slice_dir[0];
                            im_header.slice_dir[1] = acq_header.slice_dir[1];
                            im_header.slice_dir[2] = acq_header.slice_dir[2];

                            im_header.patient_table_position[0] = acq_header.patient_table_position[0];
                            im_header.patient_table_position[1] = acq_header.patient_table_position[1];
                            im_header.patient_table_position[2] = acq_header.patient_table_position[2];

                            im_header.average = acq_header.idx.average;
                            im_header.slice = acq_header.idx.slice;
                            im_header.contrast = acq_header.idx.contrast;
                            im_header.phase = acq_header.idx.phase;
                            im_header.repetition = acq_header.idx.repetition;
                            im_header.set = acq_header.idx.set;

                            im_header.acquisition_time_stamp = acq_header.acquisition_time_stamp;

                            im_header.physiology_time_stamp[0] = acq_header.physiology_time_stamp[0];
                            im_header.physiology_time_stamp[1] = acq_header.physiology_time_stamp[1];
                            im_header.physiology_time_stamp[2] = acq_header.physiology_time_stamp[2];

                            im_header.image_type = ISMRMRD::ISMRMRD_IMTYPE_COMPLEX;
                            im_header.image_index = (uint16_t)(n + s*N + slc*N*S);
                            im_header.image_series_index = 0;

                            memcpy(im_header.user_int, acq_header.user_int, sizeof(int32_t)*ISMRMRD::ISMRMRD_USER_INTS);
                            memcpy(im_header.user_float, acq_header.user_float, sizeof(float)*ISMRMRD::ISMRMRD_USER_FLOATS);

                            im_header.attribute_string_len = 0;
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericCartesianGrappaReconGadget::compute_image_header(...) ... ");
        }
    }

    GADGET_FACTORY_DECLARE(GenericCartesianGrappaReconGadget)
}
