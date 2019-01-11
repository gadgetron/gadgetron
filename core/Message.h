#pragma once

#include <memory>
#include <vector>
#include <typeindex>
#include <numeric>

namespace Gadgetron {

    class GadgetContainerMessageBase;
    template<class T> class GadgetContainerMessage;
    class LegacyGadgetNode;

    namespace Core {
        class MessageTuple;

        class Message {
        public:
            virtual ~Message() = default;


        protected:

            virtual GadgetContainerMessageBase* to_container_message() = 0;
            friend LegacyGadgetNode;
            friend MessageTuple;
        };


        template<class T>
        class TypedMessage : public Message {
        public:


            explicit TypedMessage(const T& input) : data(input) {

            }

            explicit TypedMessage(T&& input) : data(std::move(input)) {

            }

            TypedMessage(TypedMessage &&other) noexcept  = default;
            TypedMessage(const TypedMessage& other) = default;
            TypedMessage& operator=(const TypedMessage& other) = default;
            TypedMessage& operator=(TypedMessage&& other) = default;

            GadgetContainerMessageBase* to_container_message() override;

            ~TypedMessage() override =  default;

            T data;

        };

        class MessageTuple : public Message {
        public:
            template<class ...ARGS>
            MessageTuple(ARGS &&...  args) : Message(){
                static_assert(sizeof...(ARGS) > 1);
                add_messages(std::forward<ARGS>(args)...);

            }

            explicit MessageTuple(std::vector<std::unique_ptr<Message>>&& message_vector) : messages_(std::move(message_vector)){

            }

            const std::vector<std::unique_ptr<Message>> &messages() const {
                return messages_;
            }

            std::vector<std::unique_ptr<Message>>&& take_messages() {
                return std::move(messages_);
            }

            GadgetContainerMessageBase* to_container_message() override;

        private:

            template<class T> static std::unique_ptr<Message> make_message(T&& input){
                return std::make_unique<TypedMessage<T>>(std::forward<T>(input));
            }

            std::vector<std::unique_ptr<Message>> messages_;

            template<class T, class ...REST> void add_messages(T&& arg, REST&&... args){
                messages_.emplace_back(make_message(std::forward<T>(arg)));
                add_messages(std::forward<REST>(args)...);
            }

            template<class T> void add_messages(T&& arg){
                messages_.emplace_back(make_message(std::forward<T>(arg)));
            }







        };


//         template<class ...REST>
//         bool convertible_to(const Message &message);
//
//         template<class ...ARGS>
//         std::tuple<std::unique_ptr<ARGS>...> unpack(std::unique_ptr<Message> &message);
    }
}

#include "Message.hpp"