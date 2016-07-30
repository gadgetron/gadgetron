
/** \file   mri_core_dependencies.h
    \brief  Implementation useful utility functionalities for MRI dependency scan handling
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"
#include "mri_core_data.h"
#include "ismrmrd/xml.h"

namespace Gadgetron
{
    /// data: [RO E1 E2 SLC PHS CON REP SET SEG AVE]

    /// ismrmd_header : ismrmrd protocol
    /// scc_array and body_array: surface coil and body coil data, can be empty
    /// scc_header and body_header: acquisition header for surface coil and body coil data
    /// filename: file full name to save the data
    template <typename T> EXPORTMRICORE void save_dependency_data(const std::string& ismrmd_header, 
        const hoNDArray<T>& scc_array, const ISMRMRD::AcquisitionHeader& scc_header, 
        const hoNDArray<T>& body_array, const ISMRMRD::AcquisitionHeader& body_header,
        const std::string& filename);

    template <typename T> EXPORTMRICORE void load_dependency_data(const std::string& filename, 
        std::string& ismrmd_header, hoNDArray<T>& scc_array, ISMRMRD::AcquisitionHeader& scc_header, 
        hoNDArray<T>& body_array, ISMRMRD::AcquisitionHeader& body_header);
}
