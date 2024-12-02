/** \file   dicom_ismrmrd_utility.h
    \brief  Implement some utility functions to convert ismrmrd image into dicom image
    \author Hui Xue
*/

#pragma once

#include <string>
#include <map>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "Gadget.h"
#include "GadgetMRIHeaders.h"

#include "dcmtk/config/osconfig.h"
#include "dcmtk/ofstd/ofstdinc.h"
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcostrmb.h"

#include "gadgetron_dicom_export.h"
#include "hoNDArray.h"
#include "hoNDImage.h"
#include "ismrmrd/meta.h"
#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/xml.h"
#include "mri_core_def.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------------
    /// fill dicom image from ismrmrd header
    // --------------------------------------------------------------------------
    EXPORTGADGETSDICOM void fill_dicom_image_from_ismrmrd_header(const ISMRMRD::IsmrmrdHeader& h, DcmFileFormat& dcmFile);

    // --------------------------------------------------------------------------
    /// write a key and its value into dicom image
    // --------------------------------------------------------------------------
    EXPORTGADGETSDICOM void write_dcm_string(DcmDataset *dataset, DcmTagKey& key, const char* s);

    // --------------------------------------------------------------------------
    /// write ismrmrd image into a dcm image
    // --------------------------------------------------------------------------
    template<typename T> EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<T>& m2, std::string& seriesIUID, DcmFileFormat& dcmFile);
    // with image attribute
    template<typename T> EXPORTGADGETSDICOM void write_ismrmd_image_into_dicom(const ISMRMRD::ImageHeader& m1, const hoNDArray<T>& m2, ISMRMRD::IsmrmrdHeader& h, ISMRMRD::MetaContainer& attrib, std::string& seriesIUID, DcmFileFormat& dcmFile);
}
