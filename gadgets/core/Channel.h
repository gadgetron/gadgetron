#pragma once

#include "Message.h"
#include <thread>
#include <list>
#include <memory>
#include <mutex>
#include <condition_variable>

namespace Gadgetron::Core {


    template<class T>
    class InputChannel {
    public:
        virtual std::unique_ptr<T> pop() = 0;


    public:
        class Iterator;
    };

    template<class T>
    typename InputChannel<T>::Iterator begin(InputChannel<T> &);

    template<class T>
    typename InputChannel<T>::Iterator end(InputChannel<T> &);


    class OutputChannel {
    public:
        template<class T, class U = std::enable_if<!std::is_same_v<T,Message>>>
        void push(std::unique_ptr<T> &&);

        void push(std::unique_ptr<Message>&&);

        virtual void close() = 0;

    protected:
        virtual void push_message(std::unique_ptr<Message> &&) = 0;

    public:
        class Iterator;

        friend Iterator;
    };

    OutputChannel::Iterator begin(OutputChannel &);


    class MessageChannel : public OutputChannel, public InputChannel<Message> {
    public:
        virtual std::unique_ptr<Message> pop() override;

        virtual void close() override;

    protected:
        virtual void push_message(std::unique_ptr<Message> &&) override;

        std::list<std::unique_ptr<Message>> queue;
        std::mutex m;
        std::condition_variable cv;
        bool closed = false;


    };

    template<class T>
    class TypedInputChannel : public InputChannel<T> {
    public:
        TypedInputChannel(std::shared_ptr<InputChannel<Message>> input, std::shared_ptr<OutputChannel> bypass);

        virtual std::unique_ptr<T> pop();

    private:
        std::shared_ptr<InputChannel<Message>> in;
        std::shared_ptr<OutputChannel> bypass;
    };


    class ChannelClosedError : public std::runtime_error {
    public:
        ChannelClosedError() : std::runtime_error("Channel was closed") {};
    };


}


#include "Channel.hpp"
