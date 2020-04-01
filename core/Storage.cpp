#include "Storage.h"


namespace Gadgetron::Core {
    StorageSpace::StorageSpace(std::unique_ptr<StreamProvider> provider) : provider(std::move(provider)) {}
}