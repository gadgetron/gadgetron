#ifndef GADGETRON_GADGETRON_HOME_H
#define GADGETRON_GADGETRON_HOME_H

#include "gadgetbase_export.h"
#include <boost/filesystem.hpp>

namespace Gadgetron {
    boost::filesystem::path EXPORTGADGETBASE get_gadgetron_home();
}

#endif //GADGETRON_GADGETRON_HOME_H
