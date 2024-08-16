#include "ImageFinishGadget.h"

namespace Gadgetron {

    void ImageFinishGadget::process(Core::GenericInputChannel& in,
                                    Core::OutputChannel& out) {

        for (auto message : in) {
            // since we don't really do anything with text messages, received texts are printed out for the record
            if (Gadgetron::Core::convertible_to<std::string>(message))
            {
                std::string str = Gadgetron::Core::force_unpack<std::string>(std::move(message));
                GDEBUG_STREAM("Receive text message : " << str);
            }
            else
                out.push_message(std::move(message));
        }
    }

    GADGETRON_GADGET_EXPORT(ImageFinishGadget);
}
