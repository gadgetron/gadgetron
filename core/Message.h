#pragma once

#include <memory>
#include <vector>
#include <typeindex>
#include <numeric>
#include "Types.h"

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

            template<class... ARGS>
            explicit TypedMessage(ARGS&& ... xs) : data(std::forward<ARGS>(xs)...) {}

            TypedMessage(TypedMessage &&other) = default;
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
                return std::make_unique<TypedMessage<std::remove_reference_t<T>>>(std::forward<T>(input));
            }

            std::vector<std::unique_ptr<Message>> messages_;

          template<class... VARGS, class ...REST> void add_message(variant<VARGS...> var, REST&&... args){
              boost::apply_visitor([&](auto val){add_messages(val,std::forward<REST>(args)...);},var);
          }

          template<class... TARGS, class ...REST> void add_messages(tuple<TARGS...> opt, REST&&... args){
                Core::apply([&](auto... targs){add_messages(std::move(targs)...,std::forward<REST...>(args)...); },opt);
            }

            template<class T, class ...REST> void add_messages(optional<T> opt, REST&&... args){
                if (opt) messages_.emplace_back(make_message(*opt));
                add_messages(std::forward<REST>(args)...);
            }

            template<class T, class ...REST> void add_messages(T&& arg, REST&&... args){
                messages_.emplace_back(make_message(std::forward<T>(arg)));
                add_messages(std::forward<REST>(args)...);
            }

            void add_messages(){
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