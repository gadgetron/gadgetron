#pragma once

#include "Message.h"

namespace Gadgetron { namespace Core {



    template<class T> class InputChannel {
    public:
        virtual std::unique_ptr<T>&& pop() = 0;


    protected:
        virtual bool is_open() = 0;

    public:
        class Iterator;
        friend Iterator;
    };

    template<class T> InputChannel<T>::Iterator begin(InputChannel<T>&);
    template<class T> InputChannel<T>::Iterator end(InputChannel<T>&);

    class OutputChannel {
    public:
        template<class T> void push(std::unique_ptr<T>&& );

    protected:
        virtual void push_message(std::unique_ptr<Message>&&) = 0;

    public:
        class Iterator;
        friend Iterator;
    };

    OutputChannel::Iterator begin(OutputChannel&);


    template<class T> class InputChannel<T>::Iterator {
    public:
        Iterator(InputChannel<T>* c) : channel(*c){

        }
        Iterator() : channel(nullptr){}

        Iterator& operator++(){
            if (channel->is_open()){
                element = channel->pop();
            } else {
                channel = nullptr;
            }
            return *this;
        };


        bool operator==(const Iterator& other) const{
            return this->channel == other.channel;
        }

        bool operator!=(const Iterator& other) const{
            return this->channel != other.channel;
        }
        std::unique_ptr<T> operator*(){
            return std::move(element);
        }

    private:
        InputChannel* channel;
        std::unique_ptr<T> element;
    };


    template<class T> InputChannel<T>::Iterator begin(InputChannel<T>& channel){
        return InputChannel<T>::Iterator(&c);
    }
    template<class T> InputChannel<T>::Iterator end(InputChannel<T>&) {
        return InputChannel<T>::Iterator();
    }

    class OutputChannel::Iterator {
    public:

    private:
        OutputChannel* channel;

    };


}}
