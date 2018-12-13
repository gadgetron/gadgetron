#include "ImageFinishGadget.h"

namespace Gadgetron {

    void ImageFinishGadget::process(std::shared_ptr<Core::InputChannel<Core::Message>> in,
                                    std::shared_ptr<Core::OutputChannel> out) {

        for (auto message : *in) {
            out->push_message(std::move(message));
        }
        out->close();
    }

    GADGETRON_GADGET_EXPORT(ImageFinishGadget);
}
