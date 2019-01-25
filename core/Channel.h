#pragma once

#include <list>
#include <memory>
#include <mutex>
#include <condition_variable>

#include "Message.h"
#include "Types.h"


namespace Gadgetron::Core {

    class InputChannel {
    public:
        virtual Message pop() = 0;
        virtual optional<Message> try_pop() = 0;

        virtual ~InputChannel() = default;
    };


    template<class CHANNEL>
    class ChannelIterator;

    ChannelIterator<InputChannel> begin(InputChannel &);
    ChannelIterator<InputChannel> end(InputChannel &);


    class OutputChannel {
    public:

        template<class ...ARGS>
        void push(ARGS&&... ptrs);
        virtual void push_message(Message) = 0;

        virtual ~OutputChannel() = default;

    };

    ChannelIterator<OutputChannel> begin(OutputChannel &);


    class Channel : public OutputChannel, public InputChannel {
    public:
        virtual void close() = 0;

        ~Channel() override  = default;
    };


    class MessageChannel : public Channel {
    public:
        Message pop() override;
        optional<Message> try_pop() override;

        void close() override;

        void push_message(Message) override;

    protected:
        Message pop_impl(std::unique_lock<std::mutex> lock);

        std::list<Message> queue;
        std::mutex m;
        std::condition_variable cv;
        bool closed = false;


    };


    template<class ...ARGS>
    class TypedInputChannel {
    public:
        TypedInputChannel(InputChannel &input, OutputChannel &bypass) : in(input), bypass(bypass) {};

        decltype(auto) pop() {
            Message message = in.pop();
            while (!convertible_to<ARGS...>(message)) {
                bypass.push_message(std::move(message));
                message = in.pop();
            }
            return force_unpack<ARGS...>(message);
        }

        optional<decltype(force_unpack<ARGS...>(Message{}))> try_pop() {

            optional<Message> message = in.try_pop();

            while(message && !convertible_to<ARGS...>(*message)) {
                bypass.push_message(std::move(*message));
                message = in.try_pop();
            }

            if (!message) return none;

            return force_unpack<ARGS...>(*message);
        }

    private:
        InputChannel &in;
        OutputChannel &bypass;
    };

    template<class ...ARGS>
    ChannelIterator<TypedInputChannel<ARGS...>> begin(TypedInputChannel<ARGS...> &);

    template<class ...ARGS>
    ChannelIterator<TypedInputChannel<ARGS...>> end(TypedInputChannel<ARGS...> &);

    class ChannelClosed : public std::runtime_error {
    public:
        ChannelClosed() : std::runtime_error("Channel was closed") {};
    };

}


#include "Channel.hpp"
