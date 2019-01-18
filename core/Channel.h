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
        virtual std::unique_ptr<Message> pop() = 0;

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


        template<class ... TARGS>
        void push(tuple<TARGS...> tuple);


        virtual void push_message(std::unique_ptr<Message> &&) = 0;

        virtual ~OutputChannel() = default;

    };

    ChannelIterator<OutputChannel> begin(OutputChannel &);


    class Channel : public OutputChannel, public InputChannel {
    public:
        virtual void close() = 0;

        virtual ~Channel() = default;
    };


    class MessageChannel : public Channel {
    public:
        virtual std::unique_ptr<Message> pop() override;

        virtual void close() override;

        void push_message(std::unique_ptr<Message> &&) override;

    protected:

        std::list<std::unique_ptr<Message>> queue;
        std::mutex m;
        std::condition_variable cv;
        bool closed = false;


    };


    template<class ...ARGS>
    class TypedInputChannel {
    public:
        TypedInputChannel(InputChannel &input, OutputChannel &bypass) : in(input), bypass(bypass) {};

        decltype(auto) pop() {
            std::unique_ptr<Message> message = in.pop();
            while (!convertible_to<ARGS...>(*message)) {
                bypass.push_message(std::move(message));
                message = in.pop();
            }
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
