#include "SaveImageToStorageServer.h"

#include "hoNDArray_math.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {
    template<class T>
    void save_image(const Image<T> &image, const Gadgetron::StorageSpaces &storageSpace, std::string uri, int duration) {
        auto data = std::get<hoNDArray<T>>(image);
        storageSpace.session->store(uri, data, std::chrono::seconds(duration));
    }
}

namespace Gadgetron::Examples {

    SaveImageToStorageServer::SaveImageToStorageServer(const Core::Context& context, const Core::GadgetProperties& props) : PureGadget(context,props) {
        auto current_ismrmrd_header = (context.header);
        storageSpace = context.storage;
    }

    AnyImage SaveImageToStorageServer::process_function(AnyImage image) const {

        // Lambda, performs image saving
        auto save = [&](auto image){ 
            save_image(image, storageSpace, storage_uri, storage_duration);
        };

        visit(save, image);
        return image;
    }

    GADGETRON_GADGET_EXPORT(SaveImageToStorageServer);
}

