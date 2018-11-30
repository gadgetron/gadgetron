#pragma once

#include "Message.h"
#include <thread>
#include <list>
#include <memory>
#include <mutex>
#include <condition_variable>

namespace Gadgetron::Core {

    template<class ...ARGS> class InputChannel {
    public:
        virtual std::tuple<std::unique_ptr<ARGS>...> pop() = 0;

    public:
        class Iterator;
    };

    template<class T>
    class InputChannel<T> {
    public:
        virtual std::unique_ptr<T> pop() = 0;


    public:
        class Iterator;
    };

    template<class ...ARGS>
    typename InputChannel<ARGS...>::Iterator begin(InputChannel<ARGS...> &);

    template<class ...ARGS>
    typename InputChannel<ARGS...>::Iterator end(InputChannel<ARGS...> &);


    class OutputChannel {
    public:
        template<class T, class U = std::enable_if<!std::is_same_v<T,Message>>>
        void push(std::unique_ptr<T> &&);


        template<class T>
        void push(std::unique_ptr<TypedMessage<T>>&&);

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



    template<class ...ARGS>
    class TypedInputChannel : public InputChannel<ARGS...> {
    public:
        TypedInputChannel(std::shared_ptr<InputChannel<Message>> input, std::shared_ptr<OutputChannel> bypass);

        virtual std::tuple<std::unique_ptr<ARGS>...> pop() override;

    private:
        std::shared_ptr<InputChannel<Message>> in;
        std::shared_ptr<OutputChannel> bypass;
    };

    template<class T>
    class TypedInputChannel<T> : public InputChannel<T> {
    public:
        TypedInputChannel(std::shared_ptr<InputChannel<Message>> input, std::shared_ptr<OutputChannel> bypass);

        virtual std::unique_ptr<T> pop() override;

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
