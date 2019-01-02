#include "ImageFinishGadget.h"

namespace Gadgetron {

    void ImageFinishGadget::process(Core::InputChannel& in,
                                    Core::OutputChannel& out) {

        for (auto message : in) {
            out.push_message(std::move(message));
        }
    }

    GADGETRON_GADGET_EXPORT(ImageFinishGadget);
}
