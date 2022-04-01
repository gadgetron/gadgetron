#pragma once

#include "PureGadget.h"
#include "Types.h"

namespace Gadgetron::Examples {
    class SaveImageToStorageServer : public Core::PureGadget<Core::AnyImage, Core::AnyImage> {
    public:
      using Core::PureGadget<Core::AnyImage,Core::AnyImage>::PureGadget;
        SaveImageToStorageServer(const Core::Context& context, const Core::GadgetProperties& props);
        ~SaveImageToStorageServer() override = default;
        Core::AnyImage process_function(Core::AnyImage image) const override;
    protected:
      Gadgetron::StorageSpaces storageSpace;
      NODE_PROPERTY(storage_uri, std::string, "URI at which to save the image", "SaveImageSample");
      NODE_PROPERTY(storage_duration, int, "Duration to store the image in seconds", 600);
    };
}
