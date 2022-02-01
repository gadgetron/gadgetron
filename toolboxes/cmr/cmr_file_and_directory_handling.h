/** \file       cmr_slice_geometry.h
    \brief      Place to store util files
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
#include "hoNDArray_elemwise.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_reductions.h"
#include "hoNDImage.h"
#include "hoMRImage.h"

namespace Gadgetron
{
    // --------------------------------------------------------
    // create folder
    // -------------------------------------------------------------------
    EXPORTCMR void create_folder_with_all_permissions(const std::string& workingdirectory);

    // --------------------------------------------------------
    // clean items older than specified hours
    // all files created older than "hours" (in hour) are deleted
    // -------------------------------------------------------------------
    EXPORTCMR void clean_items_older_than(const std::string& workingdirectory, double hours);

    // --------------------------------------------------------
    // list items from a session
    // -------------------------------------------------------------------
    EXPORTCMR void list_items_from_session(const std::string& workingdirectory, const std::string& session_id, std::vector<std::string>& items_list);

    // --------------------------------------------------------
    // load items
    // -------------------------------------------------------------------
    template<typename T>
    EXPORTCMR void load_items(const std::string& workingdirectory, const std::vector<std::string>& items_list, std::vector<hoNDArray<T> >& items);

    // --------------------------------------------------------
    // save items
    // -------------------------------------------------------------------
    template<typename T>
    EXPORTCMR void save_item(const std::string& workingdirectory, const std::string& session_id, const hoNDArray<T>& item);
}
