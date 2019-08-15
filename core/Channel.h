#pragma once

#include <condition_variable>
#include <list>
#include <memory>
#include <mutex>

#include "Message.h"
#include "Types.h"
#include "MPMCChannel.h"

namespace Gadgetron::Core {
class InputChannel;

class OutputChannel;
struct ChannelPair;

class Channel {
public:
    virtual ~Channel() = default;

    friend InputChannel;
    friend OutputChannel;

protected:
    virtual Message pop() = 0;

    virtual optional<Message> try_pop() = 0;

    virtual void push_message(Message) = 0;

    virtual void close() = 0;

    class Closer;
};

InputChannel split(const InputChannel& channel);
OutputChannel split(const OutputChannel& channel);

/**
 * The end of a channel which provides input
 */
class InputChannel {
public:
    InputChannel(InputChannel&& other) noexcept = default;
    InputChannel& operator=(InputChannel&& other) noexcept = default;

    /// Blocks until it can take a message from the channel
    Message pop();

    /// Nonblocking method returning a message if one is available, or None otherwise
    optional<Message> try_pop();

private:
    InputChannel(const InputChannel&) = default;

    template <class ChannelType, class... ARGS>
    friend ChannelPair make_channel(ARGS&&... args);

    friend InputChannel split(const InputChannel& channel);

    explicit InputChannel(std::shared_ptr<Channel>);

    std::shared_ptr<Channel::Closer> closer;
    std::shared_ptr<Channel> channel;
};

template <class CHANNEL>
class ChannelIterator;

ChannelIterator<InputChannel> begin(InputChannel&);

ChannelIterator<InputChannel> end(InputChannel&);

/**
 * The end of a channel which provides output. Only constructible through make_channel(args)
 */
class OutputChannel {
public:
    OutputChannel(OutputChannel&& other) noexcept = default;
    OutputChannel& operator=(OutputChannel&& other) noexcept = default;

    /// Pushes a message of type ARGS to the channel
    template <class... ARGS>
    void push(ARGS&&... ptrs);

    /// Pushes a message to the channel
    void push_message(Message);

private:
    OutputChannel(const OutputChannel&) = default;

    template <class ChannelType, class... ARGS>
    friend ChannelPair make_channel(ARGS&&... args);

    friend OutputChannel split(const OutputChannel& channel);

    explicit OutputChannel(std::shared_ptr<Channel>);

    std::shared_ptr<Channel> channel;
    std::shared_ptr<Channel::Closer> closer;
};

struct ChannelPair {
    InputChannel input;
    OutputChannel output;
};

/***
 * Creates a ChannelPair
 * @tparam ChannelType Type of Channel, typically MessageChannel
 * @tparam ARGS
 * @param args
 * @return Input and output of a channel.
 */
template <class ChannelType, class... ARGS>
ChannelPair make_channel(ARGS&&... args)
{
    auto channel = std::make_shared<ChannelType>(std::forward<ARGS>(args)...);
    return { InputChannel(channel), OutputChannel(channel) };
}

ChannelIterator<OutputChannel> begin(OutputChannel&);


class MessageChannel : public Channel {

protected:
    Message pop() override;

    optional<Message> try_pop() override;

    void close() override;

    void push_message(Message) override;

    MPMCChannel<Message> channel;

};

/***
 * A wrapper around an InputChannel. Filters the content of an Inputchannel based on the specified typelist
 * @tparam ARGS
 */
template <class... TYPELIST>
class TypedInputChannel {
public:
    TypedInputChannel(InputChannel& input, OutputChannel& bypass)
        : in(input)
        , bypass(bypass) {};
    TypedInputChannel(TypedInputChannel&& other) noexcept = default;
    TypedInputChannel& operator=(TypedInputChannel&& other) noexcept = default;

    decltype(auto) pop()
    {
        Message message = in.pop();
        while (!convertible_to<TYPELIST...>(message)) {
            bypass.push_message(std::move(message));
            message = in.pop();
        }
        return force_unpack<TYPELIST...>(std::move(message));
    }

    optional<decltype(force_unpack<TYPELIST...>(Message {}))> try_pop()
    {

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
    InputChannel& in;
    OutputChannel& bypass;
};

template <class... ARGS>
ChannelIterator<TypedInputChannel<ARGS...>> begin(TypedInputChannel<ARGS...>&);

template <class... ARGS>
ChannelIterator<TypedInputChannel<ARGS...>> end(TypedInputChannel<ARGS...>&);


}

#include "Channel.hpp"
