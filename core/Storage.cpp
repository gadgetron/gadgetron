#include "Storage.h"
#include <cpprest/http_client.h>
#include <cpprest/interopstream.h>

namespace Gadgetron::Core {
    StorageSpace::StorageSpace(std::shared_ptr<StreamProvider> provider) : provider(std::move(provider)) {}
}

