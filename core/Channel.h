#pragma once

#include "Message.h"
#include <thread>
#include <list>
#include <memory>
#include <mutex>
#include <condition_variable>

namespace Gadgetron::Core {


    class InputChannel {
    public:
        virtual Message pop() = 0;

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

        void close() override;

        void push_message(Message) override;

    protected:

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
