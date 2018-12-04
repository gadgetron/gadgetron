#include "Message.h"
#include "GadgetContainerMessage.h"

Gadgetron::GadgetContainerMessageBase *Gadgetron::Core::MessageTuple::to_container_message() {

    auto front = messages_.front()->to_container_message();

    std::accumulate(messages_.begin()++, messages_.end(), front,
                    [](auto &container_message, auto &old_message) {

                        auto container_message2 = old_message->to_container_message();
                        container_message->cont(container_message2);
                        return container_message2;
                    });

    messages_.clear();
    return front;
}

