//
// Created by dchansen on 10/23/18.
//




template<class Derived> size_t Gadgetron::ImageArraySendMixin<Derived>::compute_image_number(const ISMRMRD::ImageHeader& imheader, size_t encoding, size_t CHA, size_t cha, size_t E2) const
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

template<class Derived> void Gadgetron::ImageArraySendMixin<Derived>::prep_image_header_send_out(IsmrmrdImageArray& res, size_t n, size_t s, size_t slc, size_t encoding, int series_num, const std::string& data_role) const
{

        size_t RO = res.data_.get_size(0);
        size_t E1 = res.data_.get_size(1);
        size_t E2 = res.data_.get_size(2);
        size_t CHA = res.data_.get_size(3);
        size_t N = res.data_.get_size(4);
        size_t S = res.data_.get_size(5);
        size_t SLC = res.data_.get_size(6);

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
        else if (data_role == GADGETRON_IMAGE_BINNED)
        {
            res.headers_(n, s, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

            res.meta_[offset].append(GADGETRON_IMAGECOMMENT, "Binned_Image");

            res.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, "Binned_Image");
            res.meta_[offset].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_BINNED);
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
        else if (data_role == GADGETRON_IMAGE_RETRO)
        {
            res.headers_(n, s, slc).image_type = ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE;

            res.meta_[offset].append(GADGETRON_IMAGECOMMENT, "RETRO");
            res.meta_[offset].append(GADGETRON_SEQUENCEDESCRIPTION, "RETRO");
            res.meta_[offset].set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_RETRO);
        }

}

template<class Derived> void Gadgetron::ImageArraySendMixin<Derived>::prepare_image_array(IsmrmrdImageArray& res, size_t encoding, int series_num, const std::string& data_role) const
{
        size_t RO = res.data_.get_size(0);
        size_t E1 = res.data_.get_size(1);
        size_t E2 = res.data_.get_size(2);
        size_t CHA = res.data_.get_size(3);
        size_t N = res.data_.get_size(4);
        size_t S = res.data_.get_size(5);
        size_t SLC = res.data_.get_size(6);

        GDEBUG_CONDITION_STREAM(true, "sending out image array, acquisition boundary [RO E1 E2 CHA N S SLC] = [" << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC << "] ");
        const Derived* derived = static_cast<const Derived*>(this);
        // compute image numbers and fill the image meta
        size_t n, s, slc;
        for (slc = 0; slc < SLC; slc++)
        {
            for (s = 0; s < S; s++)
            {
                for (n = 0; n < N; n++)
                {
                    this->prep_image_header_send_out(res, n, s, slc, encoding, series_num, data_role);

                    if (derived->verbose)
                    {
                        for (size_t cha = 0; cha < CHA; cha++)
                        {
                            GDEBUG_STREAM("sending out " << data_role << " image [CHA SLC CON PHS REP SET AVE] = [" << cha << " "<< res.headers_(n, s, slc).slice << " " << res.headers_(n, s, slc).contrast << " "<< res.headers_(n, s, slc).phase << " " << res.headers_(n, s, slc).repetition << " " << res.headers_(n, s, slc).set << " " << res.headers_(n, s, slc).average << " " << "] " << " -- Image series -- " << res.headers_(n, s, slc).image_series_index << " -- Image number -- " << res.headers_(n, s, slc).image_index);
                        }
                    }
                }
            }
        }


}

template<class Derived>
void Gadgetron::ImageArraySendMixin<Derived>::initialize_encoding_space_limits(const ISMRMRD::IsmrmrdHeader &h) {

    meas_max_idx_ = std::vector<ISMRMRD::EncodingCounters>(h.encoding.size());
    for (size_t e = 0; e < h.encoding.size(); e++) {
        ISMRMRD::EncodingSpace e_space = h.encoding[e].encodedSpace;
        ISMRMRD::EncodingSpace r_space = h.encoding[e].reconSpace;
        ISMRMRD::EncodingLimits e_limits = h.encoding[e].encodingLimits;

        Derived *derived = static_cast<Derived *>(this);
        GDEBUG_CONDITION_STREAM(derived->verbose, "---> Encoding space : " << e << " <---");
        GDEBUG_CONDITION_STREAM(derived->verbose,
                                "Encoding matrix size: " << e_space.matrixSize.x << " " << e_space.matrixSize.y << " "
                                                         << e_space.matrixSize.z);
        GDEBUG_CONDITION_STREAM(derived->verbose, "Encoding field_of_view : " << e_space.fieldOfView_mm.x << " "
                                                                              << e_space.fieldOfView_mm.y << " "
                                                                              << e_space.fieldOfView_mm.z);
        GDEBUG_CONDITION_STREAM(derived->verbose,
                                "Recon matrix size : " << r_space.matrixSize.x << " " << r_space.matrixSize.y << " "
                                                       << r_space.matrixSize.z);
        GDEBUG_CONDITION_STREAM(derived->verbose,
                                "Recon field_of_view :  " << r_space.fieldOfView_mm.x << " " << r_space.fieldOfView_mm.y
                                                          << " " << r_space.fieldOfView_mm.z);

        meas_max_idx_[e].kspace_encode_step_1 = (uint16_t) e_space.matrixSize.y - 1;
        meas_max_idx_[e].set = (e_limits.set && (e_limits.set->maximum > 0)) ? e_limits.set->maximum : 0;
        meas_max_idx_[e].phase = (e_limits.phase && (e_limits.phase->maximum > 0)) ? e_limits.phase->maximum : 0;

        meas_max_idx_[e].kspace_encode_step_2 = (uint16_t) e_space.matrixSize.z - 1;

        meas_max_idx_[e].contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum
                                                                                            : 0;
        meas_max_idx_[e].slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;
        meas_max_idx_[e].repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;
        meas_max_idx_[e].slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;
        meas_max_idx_[e].average = e_limits.average ? e_limits.average->maximum : 0;
        meas_max_idx_[e].segment = 0;
    }

}
