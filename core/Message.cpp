#include "Message.h"
#include "GadgetContainerMessage.h"

Gadgetron::GadgetContainerMessageBase *Gadgetron::Core::Message::to_container_message() {

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



Gadgetron::Core::Message::Message(std::vector<std::unique_ptr<Gadgetron::Core::MessageChunk>> message_vector)
        : messages_(std::move(message_vector)) {

}

const std::vector<std::unique_ptr<Gadgetron::Core::MessageChunk>> &Gadgetron::Core::Message::messages() const {
    return messages_;
}

std::vector<std::unique_ptr<Gadgetron::Core::MessageChunk>> Gadgetron::Core::Message::take_messages() {
    return std::move(messages_);
}

