/** \file       cmr_ismrmrd_util.h
    \brief      Util functions for ismrmrd
    \author     Hui Xue
*/

#pragma once

#define NOMINMAX

#include <cstdio>
#include <complex>

#include "cmr_export.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/meta.h"

#include "hoNDArray.h"

namespace Gadgetron
{
    EXPORTCMR void set_attrib_from_ismrmrd_header(const ISMRMRD::ISMRMRD_ImageHeader& header, ISMRMRD::MetaContainer& attrib);

    // get a vector of values from ismrmrd meta
    EXPORTCMR bool getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<long>& v);
    EXPORTCMR bool getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<double>& v);
    EXPORTCMR bool getISMRMRMetaValues(const ISMRMRD::MetaContainer& attrib, const std::string& name, std::vector<std::string>& v);

    template <typename T> bool setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v);
    EXPORTCMR bool setISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v);

    template <typename T> bool appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<T>& v);
    EXPORTCMR bool appendISMRMRMetaValues(ISMRMRD::MetaContainer& attrib, const std::string& name, const std::vector<std::string>& v);

    EXPORTCMR void create_image_name_from_header(const ISMRMRD::ImageHeader& head, std::string& name);
    EXPORTCMR void create_image_name_from_header(const ISMRMRD::ImageHeader& head, ISMRMRD::MetaContainer& attrib, std::string& name);

    /// initialize image header from acquisition header
    EXPORTCMR void initialize_image_header_from_acq_header(const ISMRMRD::ISMRMRD_AcquisitionHeader& acq_h, ISMRMRD::ISMRMRD_ImageHeader& im_h);
}
