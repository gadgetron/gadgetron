#include "ImageFinishGadget.h"

namespace Gadgetron{

    void ImageFinishGadget::process(std::shared_ptr<Core::InputChannel<Core::Message>> in, std::shared_ptr<Core::OutputChannel> out) {

        try {
//        for(auto message : *in ){
            while (true) {
                auto message = in->pop();
                GDEBUG("ImageFinish got here");
                out->push_message(std::move(message));
            }
        } catch(Core::ChannelClosedError err){
            out->close();

        }

    }


    GADGETRON_GADGET_EXPORT(ImageFinishGadget);
}
