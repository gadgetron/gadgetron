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

            virtual std::type_index type() = 0;

        protected:

            virtual GadgetContainerMessageBase* to_container_message() = 0;
            friend LegacyGadgetNode;
            friend MessageTuple;
        };


        template<class T>
        class TypedMessage : public Message {
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

            virtual std::type_index type() override final {
                return std::type_index(typeid(T));
            }


        protected:
            std::unique_ptr<T> data;

        };

        class MessageTuple : public Message {
        public:
            template<class ...ARGS>
            explicit MessageTuple(ARGS &&...  args) : Message(){
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

            virtual std::type_index type() override {
                return std::type_index(typeid(messages_));
            }

            virtual GadgetContainerMessageBase* to_container_message() override;

        private:

            template<class T> static std::unique_ptr<Message> make_message(T&& input){
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
    }
}

#include "Message.hpp"