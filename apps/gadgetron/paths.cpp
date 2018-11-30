#include <boost/filesystem/operations.hpp>

#include "gadgetron_home.h"
#include "log.h"
#include "paths.h"

#ifdef _WIN32
const boost::filesystem::path default_working_folder() {
    return "c:/temp/gadgetron/";
}
#else
const boost::filesystem::path default_working_folder() {
    return "/tmp/gadgetron/";
}
#endif

const boost::filesystem::path default_gadgetron_home() {
    return Gadgetron::get_gadgetron_home();
}

