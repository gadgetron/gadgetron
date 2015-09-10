/** \file   mri_core_recon_para.h
    \brief  Define the structures storing the reconstruction parameters for different algorithms
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"

#include "mri_core_export.h"
#include "mri_core_data.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------------
    /// naming conversion functions for ISMRMRD dimensions
    // --------------------------------------------------------------------------
    EXPORTMRICORE std::string get_ismrmrd_dimension_name(IsmrmrdCONDITION v);
    EXPORTMRICORE IsmrmrdCONDITION get_ismrmrd_dimension(const std::string& name);

    // --------------------------------------------------------------------------
    /// define the calibration mode of ISMRMRD
    // --------------------------------------------------------------------------
    enum ismrmrdCALIBMODE
    {
        ISMRMRD_embedded,
        ISMRMRD_interleaved,
        ISMRMRD_separate,
        ISMRMRD_external,
        ISMRMRD_other,
        ISMRMRD_noacceleration
    };

    EXPORTMRICORE std::string get_ismrmrd_calib_mode_name(ismrmrdCALIBMODE v);
    EXPORTMRICORE ismrmrdCALIBMODE get_ismrmrd_calib_mode(const std::string& name);
}
