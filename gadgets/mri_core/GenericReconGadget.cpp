
#include "GenericReconGadget.h"
#include "mri_core_kspace_filter.h"
#include "hoNDArray_reductions.h"

namespace Gadgetron {

    GenericReconGadget::GenericReconGadget() : BaseClass()
    {
    }

    GenericReconGadget::~GenericReconGadget()
    {
    }

    int GenericReconGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

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
            meas_max_idx_[e].set                  = (e_limits.set) ? e_limits.set->maximum : 0;
            meas_max_idx_[e].phase                = (e_limits.phase) ? e_limits.phase->maximum : 0;

            meas_max_idx_[e].kspace_encode_step_2 = (uint16_t)e_space.matrixSize.z - 1;

            meas_max_idx_[e].contrast             = (e_limits.contrast) ? e_limits.contrast->maximum : 0;
            meas_max_idx_[e].slice                = (e_limits.slice) ? e_limits.slice->maximum : 0;
            meas_max_idx_[e].repetition           = (e_limits.repetition) ? e_limits.repetition->maximum : 0;
            meas_max_idx_[e].slice                = (e_limits.slice) ? e_limits.slice->maximum : 0;
            meas_max_idx_[e].average              = (e_limits.average) ? e_limits.average->maximum : 0;
            meas_max_idx_[e].segment              = 0;

            if (!h.encoding[e].parallelImaging)
            {
                GDEBUG_STREAM("Parallel Imaging section not found in header for encoding space " << e);
                calib_mode_[e] = ISMRMRD_noacceleration;
                acceFactorE1_[e] = 1;
                acceFactorE2_[e] = 1;
            }
            else
            {
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
        }

        return GADGET_OK;
    }

    int GenericReconGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        process_called_times_++;

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        // for every encoding space
        for (size_t e = 0; e < recon_bit_->rbit_.size(); e++)
        {
            std::stringstream os;
            os << "_encoding_" << e;

            GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");

            // ---------------------------------------------------------------
            // export incoming data

/**
            if (!debug_folder_full_path_.empty())
            {
                if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0)
                {
                    gt_exporter_.exportArrayComplex(recon_bit_->rbit_[e].data_.data_, debug_folder_full_path_ + "data" + os.str());
                }
            }

            if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].data_.trajectory_)
            {
                if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0)
                {
                    gt_exporter_.exportArray(*(recon_bit_->rbit_[e].data_.trajectory_), debug_folder_full_path_ + "data_traj" + os.str());
                }
            }
**/

            // ---------------------------------------------------------------
            // export incoming ref

/**
            if (recon_bit_->rbit_[e].ref_)
            {
                if (!debug_folder_full_path_.empty())
                {
                    gt_exporter_.exportArrayComplex(recon_bit_->rbit_[e].ref_->data_, debug_folder_full_path_ + "ref" + os.str());
                }

                if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].ref_->trajectory_)
                {
                    if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0)
                    {
                        gt_exporter_.exportArray(*(recon_bit_->rbit_[e].ref_->trajectory_), debug_folder_full_path_ + "ref_traj" + os.str());
                    }
                }
            }
**/

            // ---------------------------------------------------------------
            // add recon code here ...
        }

        m1->release();
        return GADGET_OK;
    }

    size_t GenericReconGadget::compute_image_number(ISMRMRD::ImageHeader& imheader, size_t encoding, size_t CHA, size_t cha, size_t E2)
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

    int GenericReconGadget::send_out_image_array(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t encoding, int series_num, const std::string& data_role)
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

                            // set the skip processing flag, so gfactor map will not be processed during e.g. partial fourier handling or kspace filter gadgets
                            res.meta_[offset].set(GADGETRON_SKIP_PROCESSING_AFTER_RECON, (long)1);
                        }
                        else if (data_role == GADGETRON_IMAGE_SNR_MAP)
                        {
                            res.headers_(n, s, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

                            res.meta_[offset].append(GADGETRON_IMAGECOMMENT, GADGETRON_IMAGE_SNR_MAP);
                            res.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, GADGETRON_IMAGE_SNR_MAP);
                            res.meta_[offset].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_SNR_MAP);
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
            GERROR_STREAM("Errors in GenericReconGadget::send_out_image_array(...) ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    // ----------------------------------------------------------------------------------------

    void GenericReconGadget::make_ref_coil_map(IsmrmrdDataBuffered& ref_, std::vector<size_t>  recon_dims, hoNDArray< std::complex<float> >& ref_calib, hoNDArray< std::complex<float> >& ref_coil_map, size_t encoding)
    {

        try
        {
            hoNDArray< std::complex<float> >& ref_data = ref_.data_;

            // sampling limits
            size_t sRO = ref_.sampling_.sampling_limits_[0].min_;
            size_t eRO = ref_.sampling_.sampling_limits_[0].max_;
            size_t cRO = ref_.sampling_.sampling_limits_[0].center_;

            size_t sE1 = ref_.sampling_.sampling_limits_[1].min_;
            size_t eE1 = ref_.sampling_.sampling_limits_[1].max_;
            size_t cE1 = ref_.sampling_.sampling_limits_[1].center_;

            size_t sE2 = ref_.sampling_.sampling_limits_[2].min_;
            size_t eE2 = ref_.sampling_.sampling_limits_[2].max_;
            size_t cE2 = ref_.sampling_.sampling_limits_[2].center_;

            // recon size
            size_t recon_RO = recon_dims[0];
            size_t recon_E1 = recon_dims[1];
            size_t recon_E2 = recon_dims[2];

            // ref array size
            size_t ref_data_RO = ref_data.get_size(0);
            size_t ref_data_E1 = ref_data.get_size(1);
            size_t ref_data_E2 = ref_data.get_size(2);
            size_t CHA = ref_data.get_size(3);
            size_t N = ref_data.get_size(4);
            size_t S = ref_data.get_size(5);
            size_t SLC = ref_data.get_size(6);

            // determine the ref_coil_map size
            size_t RO = 2 * cRO;
            if (sRO>0 || eRO<RO - 1)
            {
                RO = 2 * std::max(cRO - sRO, eRO - cRO+1);
                if (RO>recon_RO) RO = recon_RO;
            }

            size_t E1 = eE1 - sE1 + 1;
            size_t E2 = eE2 - sE2 + 1;

            // cut the center region for ref coil map
            if ((calib_mode_[encoding] == Gadgetron::ISMRMRD_interleaved) || (calib_mode_[encoding] == Gadgetron::ISMRMRD_noacceleration))
            {
                E1 = 2 * std::min(cE1 - sE1, eE1 - cE1+1) - 1;
                if (E1>recon_E1) E1 = recon_E1;

                if (E2 > 1)
                {
                    E2 = 2 * std::min(cE2 - sE2, eE2 - cE2 + 1) - 1;
                    if (E2 > recon_E2) E2 = recon_E2;
                }
            }

            ref_coil_map.create(RO, E1, E2, CHA, N, S, SLC);
            Gadgetron::clear(ref_coil_map);

            // make sure center aligned
            size_t ref_sE1 = cE1 - E1 / 2;
            size_t ref_eE1 = ref_sE1 + E1 - 1;

            size_t ref_sE2 = cE2 - E2 / 2;
            size_t ref_eE2 = ref_sE2 + E2 - 1;

            size_t slc, s, n, cha, e2, e1;
            for (slc = 0; slc < SLC; slc++)
            {
                for (s = 0; s < S; s++)
                {
                    for (n = 0; n < N; n++)
                    {
                        for (cha = 0; cha < CHA; cha++)
                        {
                            for (e2 = ref_sE2; e2 <= ref_eE2; e2++)
                            {
                                for (e1 = ref_sE1; e1 <= ref_eE1; e1++)
                                {
                                    std::complex<float>* pSrc = &(ref_data(0, e1, e2, cha, n, s, slc));
                                    std::complex<float>* pDst = &(ref_coil_map(0, e1- ref_sE1, e2- ref_sE2, cha, n, s, slc));

                                    memcpy(pDst + sRO, pSrc, sizeof(std::complex<float>)*(eRO - sRO + 1));
                                }
                            }
                        }
                    }
                }
            }

/**
            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << encoding;

                gt_exporter_.exportArrayComplex(ref_coil_map, debug_folder_full_path_ + "ref_coil_map_before_filtering_" + os.str());
            }
**/

            // filter the ref_coil_map
            if (filter_RO_ref_coi_map_.get_size(0) != RO)
            {
                Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(0), ref_.sampling_.sampling_limits_[0].min_, ref_.sampling_.sampling_limits_[0].max_, filter_RO_ref_coi_map_);

/**
                if (!debug_folder_full_path_.empty())
                {
                    std::stringstream os;
                    os << "encoding_" << encoding;

                    gt_exporter_.exportArrayComplex(filter_RO_ref_coi_map_, debug_folder_full_path_ + "filter_RO_ref_coi_map_" + os.str());
                }
**/
            }

            if (filter_E1_ref_coi_map_.get_size(0) != E1)
            {
                Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(1), 0, E1-1, filter_E1_ref_coi_map_);

/**
                if (!debug_folder_full_path_.empty())
                {
                    std::stringstream os;
                    os << "encoding_" << encoding;

                    gt_exporter_.exportArrayComplex(filter_E1_ref_coi_map_, debug_folder_full_path_ + "filter_E1_ref_coi_map_" + os.str());
                }
**/
            }

            if ( (E2 > 1) && (filter_E2_ref_coi_map_.get_size(0) != E2) )
            {
                Gadgetron::generate_symmetric_filter_ref(ref_coil_map.get_size(2), 0, E2-1, filter_E2_ref_coi_map_);

/**
                if (!debug_folder_full_path_.empty())
                {
                    std::stringstream os;
                    os << "encoding_" << encoding;

                    gt_exporter_.exportArrayComplex(filter_E2_ref_coi_map_, debug_folder_full_path_ + "filter_E2_ref_coi_map_" + os.str());
                }
**/
            }

            hoNDArray< std::complex<float> > ref_recon_buf;

            if (E2 > 1)
            {
                Gadgetron::apply_kspace_filter_ROE1E2(ref_coil_map, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, filter_E2_ref_coi_map_, ref_recon_buf);
            }
            else
            {
                Gadgetron::apply_kspace_filter_ROE1(ref_coil_map, filter_RO_ref_coi_map_, filter_E1_ref_coi_map_, ref_recon_buf);
            }

/**
            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << encoding;

                gt_exporter_.exportArrayComplex(ref_recon_buf, debug_folder_full_path_ + "ref_coil_map_after_filtering_" + os.str());
            }
**/

            // pad the ref_coil_map into the data array
            Gadgetron::pad(recon_RO, recon_E1, recon_E2, &ref_recon_buf, &ref_coil_map);

            std::vector<size_t> dim = *ref_data.get_dimensions();
            ref_calib.create(dim, ref_data.begin());

/**
            if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << encoding;

                gt_exporter_.exportArrayComplex(ref_coil_map, debug_folder_full_path_ + "ref_coil_map_" + os.str());
                gt_exporter_.exportArrayComplex(ref_calib, debug_folder_full_path_ + "ref_calib_" + os.str());
            }
**/
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconGadget::make_ref_coil_map(...) ... ");
        }
    }

    void GenericReconGadget::perform_coil_map_estimation(const hoNDArray< std::complex<float> >& ref_coil_map, hoNDArray< std::complex<float> >& coil_map, size_t e)
    {
        try
        {
            coil_map = ref_coil_map;
            Gadgetron::clear(coil_map);

            size_t E2 = ref_coil_map.get_size(2);
            if (E2 > 1)
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft3c(ref_coil_map, complex_im_recon_buf_);
            }
            else
            {
                Gadgetron::hoNDFFT<float>::instance()->ifft2c(ref_coil_map, complex_im_recon_buf_);
            }

            /*if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << e;

                gt_exporter_.exportArrayComplex(complex_im_recon_buf_, debug_folder_full_path_ + "complex_im_for_coil_map_" + os.str());
            }*/

            if (coil_map_algorithm.value() == "Inati")
            {
                size_t ks = 7;
                size_t kz = 5;
                size_t power = 3;

                Gadgetron::coil_map_Inati(complex_im_recon_buf_, coil_map, ks, kz, power);
            }
            else
            {
                size_t ks = 7;
                size_t kz = 5;
                size_t iterNum = 5;
                float thres = 0.001;

                Gadgetron::coil_map_Inati_Iter(complex_im_recon_buf_, coil_map, ks, kz, iterNum, thres);
            }

            /*if (!debug_folder_full_path_.empty())
            {
                std::stringstream os;
                os << "encoding_" << e;

                gt_exporter_.exportArrayComplex(coil_map, debug_folder_full_path_ + "coil_map_" + os.str());
            }*/
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconGadget::perform_coil_map_estimation(...) ... ");
        }
    }

    void GenericReconGadget::compute_image_header(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t e)
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

            GADGET_CHECK_THROW(N == recon_bit.data_.headers_.get_size(2));
            GADGET_CHECK_THROW(S == recon_bit.data_.headers_.get_size(3));

            res.headers_.create(N, S, SLC);
            res.meta_.resize(N*S*SLC);

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

                        ISMRMRD::ImageHeader& im_header = res.headers_(n, s, slc);
                        ISMRMRD::MetaContainer& meta = res.meta_[n + s*N + slc*N*S];

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

                        im_header.image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;
                        im_header.image_index = (uint16_t)(n + s*N + slc*N*S);
                        im_header.image_series_index = 0;

                        memcpy(im_header.user_int, acq_header.user_int, sizeof(int32_t)*ISMRMRD::ISMRMRD_USER_INTS);
                        memcpy(im_header.user_float, acq_header.user_float, sizeof(float)*ISMRMRD::ISMRMRD_USER_FLOATS);

                        im_header.attribute_string_len = 0;

                        meta.set("encoding", (long)e);

                        meta.set("encoding_FOV"         , recon_bit.data_.sampling_.encoded_FOV_[0]);
                        meta.append("encoding_FOV"      , recon_bit.data_.sampling_.encoded_FOV_[1]);
                        meta.append("encoding_FOV"      , recon_bit.data_.sampling_.encoded_FOV_[2]);

                        meta.set("recon_FOV"            , recon_bit.data_.sampling_.recon_FOV_[0]);
                        meta.append("recon_FOV"         , recon_bit.data_.sampling_.recon_FOV_[1]);
                        meta.append("recon_FOV"         , recon_bit.data_.sampling_.recon_FOV_[2]);

                        meta.set("encoded_matrix"       , (long)recon_bit.data_.sampling_.encoded_matrix_[0]);
                        meta.append("encoded_matrix"    , (long)recon_bit.data_.sampling_.encoded_matrix_[1]);
                        meta.append("encoded_matrix"    , (long)recon_bit.data_.sampling_.encoded_matrix_[2]);

                        meta.set("recon_matrix"         , (long)recon_bit.data_.sampling_.recon_matrix_[0]);
                        meta.append("recon_matrix"      , (long)recon_bit.data_.sampling_.recon_matrix_[1]);
                        meta.append("recon_matrix"      , (long)recon_bit.data_.sampling_.recon_matrix_[2]);

                        meta.set("sampling_limits_RO"   , (long)recon_bit.data_.sampling_.sampling_limits_[0].min_);
                        meta.append("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].center_);
                        meta.append("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].max_);

                        meta.set("sampling_limits_E1"   , (long)recon_bit.data_.sampling_.sampling_limits_[1].min_);
                        meta.append("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].center_);
                        meta.append("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].max_);

                        meta.set("sampling_limits_E2"   , (long)recon_bit.data_.sampling_.sampling_limits_[2].min_);
                        meta.append("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].center_);
                        meta.append("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].max_);
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in GenericReconGadget::compute_image_header(...) ... ");
        }
    }

    GADGET_FACTORY_DECLARE(GenericReconGadget)
}
