#include "Storage.h"

namespace Gadgetron::Core {
    StorageSpace::StorageSpace(std::shared_ptr<StreamProvider> provider) : provider(std::move(provider)) {}
}
