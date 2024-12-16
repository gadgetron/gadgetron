//
// Created by dchansen on 10/23/18.
//
#pragma once

#include <vector>


namespace Gadgetron {

    template<class Derived> class ImageArraySendMixin
    {
    protected:
        void initialize_encoding_space_limits(const mrd::Header& header);
        // encoding space limits

        // --------------------------------------------------
        // utility functions
        // --------------------------------------------------
        // compute image number
        size_t compute_image_number(const mrd::ImageHeader& imheader, size_t encoding = 0, size_t CHA = 1, size_t cha = 0, size_t E2 = 1) const ;

        // prepare header for sending out
        void prep_image_header_send_out(mrd::ImageArray& res, size_t n, size_t s, size_t slc, size_t encoding, int series_num, const std::string& data_role) const ;

        // send out the recon results
        void prepare_image_array(mrd::ImageArray& res, size_t encoding, int series_num, const std::string& data_role) const ;

    private:
        std::vector<mrd::EncodingCounters> meas_max_idx_;
    };
}

#include "ImageArraySendMixin.hpp"
