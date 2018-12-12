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
            virtual ~Message(){};


        protected:

            virtual GadgetContainerMessageBase* to_container_message() = 0;
            friend LegacyGadgetNode;
            friend MessageTuple;
        };

         template<class T, typename = std::enable_if_t<!std::is_convertible_v<T*,Message*>>>
        class TypedMessage;

        template<class T>
        class TypedMessage<T> : public Message {
        public:


            TypedMessage(const T& input) : data(std::make_unique<T>(input)) {

            }

            TypedMessage(T&& input) : data(std::make_unique<T>(std::move(input))) {

            }

            TypedMessage(std::unique_ptr<T> &&input_ptr) : data(std::move(input_ptr)) {

            }

            TypedMessage(TypedMessage &&other) : data(other.get_data()) {

            }

            TypedMessage(TypedMessage& other) = delete;

            std::unique_ptr<T> &&take_data() {
                return std::move(data);
            }

            virtual GadgetContainerMessageBase* to_container_message() override;

            virtual ~TypedMessage() {};

        protected:
            std::unique_ptr<T> data;

        };

        class MessageTuple : public Message {
        public:


            template<class ...ARGS>
            MessageTuple(ARGS &&...  args) : Message(){
                static_assert(sizeof...(ARGS) > 1);
                add_messages(std::move(args)...);

            }

            explicit MessageTuple(std::vector<std::unique_ptr<Message>>&& message_vector) : messages_(std::move(message_vector)){

            }

            const std::vector<std::unique_ptr<Message>> &messages() const {
                return messages_;
            }

            std::vector<std::unique_ptr<Message>> &&take_messages() {
                return std::move(messages_);
            }



            virtual GadgetContainerMessageBase* to_container_message() override;

        private:

            template<class T> static std::unique_ptr<Message> make_message(std::unique_ptr<T>&& input){
                return std::make_unique<TypedMessage<T>>(std::move(input));
            }

            std::vector<std::unique_ptr<Message>> messages_;

            template<class T, class ...REST> void add_messages(T&& arg, REST&&... args){
                messages_.emplace_back(make_message(std::move(arg)));
                add_messages(std::move(args)...);
            }

            template<class T> void add_messages(T&& arg){
                messages_.emplace_back(make_message(std::move(arg)));
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