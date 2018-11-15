#pragma once

#include "Message.h"
#include <thread>
#include <list>
#include <memory>
#include <mutex>
#include <condition_variable>

namespace Gadgetron {
    namespace Core {


        template<class T>
        class InputChannel {
        public:
            virtual std::unique_ptr <T> &&pop() = 0;


        protected:
            virtual bool is_open() = 0;

        public:
            class Iterator;

            friend Iterator;
        };

        template<class T>
        typename InputChannel<T>::Iterator begin(InputChannel<T> &);

        template<class T>
        typename InputChannel<T>::Iterator end(InputChannel<T> &);

        class OutputChannel {
        public:
            template<class T>
            void push(std::unique_ptr <T> &&);

        protected:
            virtual void push_message(std::unique_ptr <Message> &&) = 0;

        public:
            class Iterator;

            friend Iterator;
        };

        OutputChannel::Iterator begin(OutputChannel &);


        class MessageChannel : public OutputChannel, public InputChannel<Message> {
        public:
            virtual std::unique_ptr<Message>&& pop() override;
            
        protected:
            virtual bool is_open() override;
            virtual void push_message(std::unique_ptr<Message>&& ) override;

            std::list<std::unique_ptr<Message>> queue;
            std::mutex m;
            std::condition_variable cv;
            bool closed;


        };


    }
}
#include "Channel.hpp"
