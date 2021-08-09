#pragma once

#include <boost/filesystem/path.hpp>

namespace Gadgetron::Server {
    const boost::filesystem::path default_working_folder();
    const boost::filesystem::path default_gadgetron_home();
    const boost::filesystem::path default_database_folder();
    const boost::filesystem::path default_storage_folder();
}


