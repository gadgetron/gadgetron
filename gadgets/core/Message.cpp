#include "Message.h"
#include "GadgetContainerMessage.h"

Gadgetron::GadgetContainerMessageBase *Gadgetron::Core::MessageTuple::to_container_message() {

    auto back = messages_.back()->to_container_message();

    auto result = std::accumulate(++messages_.rbegin(), messages_.rend(), back,
                    [](auto &container_message, auto &old_message) {
                        auto container_message2 = old_message->to_container_message();
                        container_message2->cont(container_message);
                        return container_message2;
                    });

    messages_.clear();
    return result;
}

