//
// Created by dchansen on 10/23/18.
//

template<class Derived> size_t Gadgetron::ImageArraySendMixin<Derived>::compute_image_number(const mrd::ImageHeader& imheader, size_t encoding, size_t CHA, size_t cha, size_t E2) const
{
    if (encoding >= meas_max_idx_.size())
    {
        GWARN_STREAM("encoding >= meas_max_idx_.size()");
        encoding = 0;
    }

    size_t SET = meas_max_idx_[encoding].set.value_or(0) + 1;
    size_t REP = meas_max_idx_[encoding].repetition.value_or(0) + 1;
    size_t PHS = meas_max_idx_[encoding].phase.value_or(0) + 1;
    size_t SLC = meas_max_idx_[encoding].slice.value_or(0) + 1;
    size_t CON = meas_max_idx_[encoding].contrast.value_or(0) + 1;
    if (E2 == 0) E2 = 1;

    size_t imageNum = (
        imheader.average.value_or(0)*REP*SET*PHS*CON*SLC*E2*CHA +
        imheader.repetition.value_or(0)*SET*PHS*CON*SLC*E2*CHA +
        imheader.set.value_or(0)*PHS*CON*SLC*E2*CHA +
        imheader.phase.value_or(0)*CON*SLC*E2*CHA +
        imheader.contrast.value_or(0)*SLC*E2*CHA +
        imheader.slice.value_or(0)*E2*CHA + cha + 1
    );

    return imageNum;
}

template<class Derived> void Gadgetron::ImageArraySendMixin<Derived>::prep_image_header_send_out(mrd::ImageArray& res, size_t n, size_t s, size_t slc, size_t encoding, int series_num, const std::string& data_role) const
{
        size_t RO = res.data.get_size(0);
        size_t E1 = res.data.get_size(1);
        size_t E2 = res.data.get_size(2);
        size_t CHA = res.data.get_size(3);
        size_t N = res.data.get_size(4);
        size_t S = res.data.get_size(5);
        size_t SLC = res.data.get_size(6);

        auto& header = res.headers(n, s, slc);
        auto& meta = res.meta(n, s, slc);

        auto image_index = this->compute_image_number(header, encoding, CHA, 0, E2);
        header.image_index = image_index;
        header.image_series_index = series_num;

        meta[GADGETRON_IMAGENUMBER] = {(long)image_index};
        meta[GADGETRON_IMAGEPROCESSINGHISTORY] = {"GT"};

        if (data_role == GADGETRON_IMAGE_REGULAR)
        {
            header.image_type = mrd::ImageType::kMagnitude;

            meta[GADGETRON_IMAGECOMMENT].push_back("GT");

            meta[GADGETRON_SEQUENCEDESCRIPTION].push_back("_GT");
            meta[GADGETRON_DATA_ROLE] = { GADGETRON_IMAGE_REGULAR };
        }
        else if (data_role == GADGETRON_IMAGE_GFACTOR)
        {
            header.image_type = mrd::ImageType::kMagnitude;

            meta[GADGETRON_IMAGECOMMENT].push_back(GADGETRON_IMAGE_GFACTOR);
            meta[GADGETRON_SEQUENCEDESCRIPTION].push_back(GADGETRON_IMAGE_GFACTOR);
            meta[GADGETRON_DATA_ROLE] = { GADGETRON_IMAGE_GFACTOR };

            // set the skip processing flag, so gfactor map will not be processed during e.g. partial fourier handling or kspace filter gadgets
            meta[GADGETRON_SKIP_PROCESSING_AFTER_RECON] = {(long)1};
        }
        else if (data_role == GADGETRON_IMAGE_SNR_MAP)
        {
            header.image_type = mrd::ImageType::kMagnitude;

            meta[GADGETRON_IMAGECOMMENT].push_back(GADGETRON_IMAGE_SNR_MAP);
            meta[GADGETRON_SEQUENCEDESCRIPTION].push_back(GADGETRON_IMAGE_SNR_MAP);
            meta[GADGETRON_DATA_ROLE] = { GADGETRON_IMAGE_SNR_MAP };
        }
        else if (data_role == GADGETRON_IMAGE_RETRO)
        {
            header.image_type = mrd::ImageType::kMagnitude;

            meta[GADGETRON_IMAGECOMMENT].push_back("RETRO");
            meta[GADGETRON_SEQUENCEDESCRIPTION].push_back("RETRO");
            meta[GADGETRON_DATA_ROLE] = { GADGETRON_IMAGE_RETRO };
        }
}

template<class Derived> void Gadgetron::ImageArraySendMixin<Derived>::prepare_image_array(mrd::ImageArray& res, size_t encoding, int series_num, const std::string& data_role) const
{
        size_t RO = res.data.get_size(0);
        size_t E1 = res.data.get_size(1);
        size_t E2 = res.data.get_size(2);
        size_t CHA = res.data.get_size(3);
        size_t N = res.data.get_size(4);
        size_t S = res.data.get_size(5);
        size_t SLC = res.data.get_size(6);

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
                        auto& header = res.headers(n, s, slc);
                        for (size_t cha = 0; cha < CHA; cha++)
                        {
                            GDEBUG_STREAM("sending out " << data_role << " image [CHA SLC CON PHS REP SET AVE] = [" << cha << " "
                                << header.slice.value_or(0) << " " << header.contrast.value_or(0) << " "<< header.phase.value_or(0) << " "
                                << header.repetition.value_or(0) << " " << header.set.value_or(0) << " " << header.average.value_or(0) << " "
                                << "] " << " -- Image series -- " << header.image_series_index.value_or(0) << " -- Image number -- " << header.image_index.value_or(0));
                        }
                    }
                }
            }
        }
}

template<class Derived>
void Gadgetron::ImageArraySendMixin<Derived>::initialize_encoding_space_limits(const mrd::Header &h)
{
    meas_max_idx_ = std::vector<mrd::EncodingCounters>(h.encoding.size());
    for (size_t e = 0; e < h.encoding.size(); e++) {
        mrd::EncodingSpaceType e_space = h.encoding[e].encoded_space;
        mrd::EncodingSpaceType r_space = h.encoding[e].recon_space;
        mrd::EncodingLimitsType e_limits = h.encoding[e].encoding_limits;

        Derived *derived = static_cast<Derived *>(this);
        GDEBUG_CONDITION_STREAM(derived->verbose, "---> Encoding space : " << e << " <---");
        GDEBUG_CONDITION_STREAM(derived->verbose,
                                "Encoding matrix size: " << e_space.matrix_size.x << " " << e_space.matrix_size.y << " "
                                                         << e_space.matrix_size.z);
        GDEBUG_CONDITION_STREAM(derived->verbose, "Encoding field_of_view : " << e_space.field_of_view_mm.x << " "
                                                                              << e_space.field_of_view_mm.y << " "
                                                                              << e_space.field_of_view_mm.z);
        GDEBUG_CONDITION_STREAM(derived->verbose,
                                "Recon matrix size : " << r_space.matrix_size.x << " " << r_space.matrix_size.y << " "
                                                       << r_space.matrix_size.z);
        GDEBUG_CONDITION_STREAM(derived->verbose,
                                "Recon field_of_view :  " << r_space.field_of_view_mm.x << " " << r_space.field_of_view_mm.y
                                                          << " " << r_space.field_of_view_mm.z);

        meas_max_idx_[e].kspace_encode_step_1 = (uint16_t) e_space.matrix_size.y - 1;
        meas_max_idx_[e].set = (e_limits.set && (e_limits.set->maximum > 0)) ? e_limits.set->maximum : 0;
        meas_max_idx_[e].phase = (e_limits.phase && (e_limits.phase->maximum > 0)) ? e_limits.phase->maximum : 0;

        meas_max_idx_[e].kspace_encode_step_2 = (uint16_t) e_space.matrix_size.z - 1;

        meas_max_idx_[e].contrast = (e_limits.contrast && (e_limits.contrast->maximum > 0)) ? e_limits.contrast->maximum
                                                                                            : 0;
        meas_max_idx_[e].slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;
        meas_max_idx_[e].repetition = e_limits.repetition ? e_limits.repetition->maximum : 0;
        meas_max_idx_[e].slice = (e_limits.slice && (e_limits.slice->maximum > 0)) ? e_limits.slice->maximum : 0;
        meas_max_idx_[e].average = e_limits.average ? e_limits.average->maximum : 0;
        meas_max_idx_[e].segment = 0;
    }

}
