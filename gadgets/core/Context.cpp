
#include "Context.h"

Gadgetron::Core::Context::Paths::Paths(const boost::filesystem::path &gadgetron_home,
                                       const boost::filesystem::path &working_folder)
                                       : gadgetron_home(gadgetron_home)
                                       , working_folder(working_folder) {}


