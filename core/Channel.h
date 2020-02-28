#pragma once

#include <condition_variable>
#include <list>
#include <memory>
#include <mutex>

#include "MPMCChannel.h"
#include "Message.h"
#include "Types.h"

#include "ChannelIterator.h"

namespace Gadgetron { namespace Core {
    class GenericInputChannel;

    class OutputChannel;
    struct ChannelPair;

    class Channel {
    public:
        virtual ~Channel() = default;

        friend GenericInputChannel;
        friend OutputChannel;

    protected:
        virtual Message pop() = 0;

        virtual optional<Message> try_pop() = 0;

        virtual void push_message(Message) = 0;

        virtual void close() = 0;

        class Closer;
    };

    class MessageChannel;
    template <class ChannelType=MessageChannel, class... ARGS> ChannelPair make_channel(ARGS&&... args);
    /**
     * The end of a channel which provides output. Only constructible through make_channel(args)
     */
    class OutputChannel {
    public:
        OutputChannel(OutputChannel&& other) noexcept = default;
        OutputChannel& operator=(OutputChannel&& other) noexcept = default;

        /// Pushes a message of type ARGS to the channel
        template <class... ARGS> void push(ARGS&&... ptrs);

        /// Pushes a message to the channel
        void push_message(Message);

        ChannelIterator<OutputChannel> begin();

    private:
        OutputChannel(const OutputChannel&) = default;

        template <class ChannelType, class... ARGS> friend ChannelPair make_channel(ARGS&&... args);

        friend OutputChannel split(const OutputChannel& channel);

        explicit OutputChannel(std::shared_ptr<Channel>);

        std::shared_ptr<Channel> channel;
        std::shared_ptr<Channel::Closer> closer;
    };
} }

#include "Channel.hpp"

namespace Gadgetron { namespace Core {

    GenericInputChannel split(const GenericInputChannel& channel);
    OutputChannel split(const OutputChannel& channel);

    /**
     * The end of a channel which provides input
     */
    class GenericInputChannel : public ChannelRange<GenericInputChannel> {
    public:
        GenericInputChannel(GenericInputChannel&& other) noexcept = default;
        GenericInputChannel& operator=(GenericInputChannel&& other) noexcept = default;

        /// Blocks until it can take a message from the channel
        Message pop();

        /// Nonblocking method returning a message if one is available, or None otherwise
        optional<Message> try_pop();

    private:
        GenericInputChannel(const GenericInputChannel&) = default;

        template <class ChannelType, class... ARGS> friend ChannelPair make_channel(ARGS&&... args);

        friend GenericInputChannel split(const GenericInputChannel& channel);

        explicit GenericInputChannel(std::shared_ptr<Channel>);

        std::shared_ptr<Channel::Closer> closer;
        std::shared_ptr<Channel> channel;
    };

    template <class CHANNEL> class ChannelIterator;

    struct ChannelPair {
        GenericInputChannel input;
        OutputChannel output;
    };
    class MessageChannel : public Channel {

    protected:
        Message pop() override;

        optional<Message> try_pop() override;

        void close() override;

        void push_message(Message) override;

        MPMCChannel<Message> channel;
    };

    /***
     * Creates a ChannelPair
     * @tparam ChannelType Type of Channel, typically MessageChannel
     * @tparam ARGS
     * @param args
     * @return Input and output of a channel.
     */
    template <class ChannelType, class... ARGS> ChannelPair make_channel(ARGS&&... args) {
        auto channel = std::make_shared<ChannelType>(std::forward<ARGS>(args)...);
        return { GenericInputChannel(channel), OutputChannel(channel) };
    }

    ChannelIterator<OutputChannel> begin(OutputChannel&);

    /***
     * A wrapper around a GenericInputChannel. Filters the content of an Inputchannel based on the specified typelist
     * @tparam ARGS
     */
    template <class... TYPELIST> class InputChannel : public ChannelRange<InputChannel<TYPELIST...>> {
    public:
        InputChannel(GenericInputChannel& input, OutputChannel& bypass) : in(input), bypass(bypass){};
        InputChannel(InputChannel&& other) noexcept = default;
        InputChannel& operator=(InputChannel&& other) noexcept = default;

        decltype(auto) pop() {
            Message message = in.pop();
            while (!convertible_to<TYPELIST...>(message)) {
                bypass.push_message(std::move(message));
                message = in.pop();
            }
            return force_unpack<TYPELIST...>(std::move(message));
        }

        optional<decltype(force_unpack<TYPELIST...>(Message{}))> try_pop() {

            optional<Message> message = in.try_pop();

            while (message && !convertible_to<TYPELIST...>(*message)) {
                bypass.push_message(std::move(*message));
                message = in.try_pop();
            }

            if (!message)
                return none;

            return force_unpack<TYPELIST...>(std::move(*message));
        }

    private:
        GenericInputChannel& in;
        OutputChannel& bypass;
    };

}}
