#pragma once

#include <memory>
#include <vector>
#include <typeindex>
#include <numeric>
#include "Types.h"

namespace Gadgetron {

    class GadgetContainerMessageBase;

    template<class T>
    class GadgetContainerMessage;

    class LegacyGadgetNode;

    namespace Core {
        class Message;

        class MessageChunk {
        public:
            virtual ~MessageChunk() = default;
            virtual std::unique_ptr<MessageChunk> clone() const = 0;
        protected:
            virtual GadgetContainerMessageBase *to_container_message() = 0;

            friend Message;
        };


        class Message {
        public:
            template<class ...ARGS>
            explicit Message(ARGS &&...  args);

            explicit Message(std::vector<std::unique_ptr<MessageChunk>> message_vector);

            Message(const Message&) = delete;
            Message(Message&&) = default;

            Message& operator=(Message&&) = default;

            const std::vector<std::unique_ptr<MessageChunk>> &messages() const;

            std::vector<std::unique_ptr<MessageChunk>> take_messages();

            GadgetContainerMessageBase *to_container_message();

            Message clone();

        private:
            std::vector<std::unique_ptr<MessageChunk>> messages_;
        };

        template<class... ARGS>
        bool convertible_to(const Message &);

        template<class ...ARGS>
        std::enable_if_t<(sizeof...(ARGS) > 1), std::tuple<ARGS...>>
        force_unpack(Message message);


        template<class T>
        T force_unpack(Message message);


        template<class ...ARGS>
        std::enable_if_t<(sizeof...(ARGS) > 1), optional<std::tuple<ARGS...>>>
        unpack(Message &&message);

        template<class T>
        optional<T> unpack(Message &&message);


        template<class T>
        class TypedMessageChunk : public MessageChunk {
        public:

            template<class... ARGS>
            explicit TypedMessageChunk(ARGS &&... xs) : data(std::forward<ARGS>(xs)...) {}

            TypedMessageChunk(TypedMessageChunk &&other) = default;

            TypedMessageChunk(const TypedMessageChunk &other) = default;

            TypedMessageChunk &operator=(const TypedMessageChunk &other) = default;

            TypedMessageChunk &operator=(TypedMessageChunk &&other) = default;

            GadgetContainerMessageBase *to_container_message() override;

            std::unique_ptr<MessageChunk> clone() const override;

            ~TypedMessageChunk() override = default;

            T data;

        };
    }
}

#include "Message.hpp"