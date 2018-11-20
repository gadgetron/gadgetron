#pragma once

namespace Gadgetron::Core {
    template<class T>
    class InputChannel<T>::Iterator {
    public:
        Iterator(InputChannel<T> *c) : channel(*c) {
            this->operator*();
        }

        Iterator() : channel(nullptr) {}

        Iterator &operator++() {
            try {
                element = channel->pop();
            } catch (ChannelClosedError err) {
                channel = nullptr;
            }
            return *this;
        };


        bool operator==(const Iterator &other) const {
            return this->channel == other.channel;
        }

        bool operator!=(const Iterator &other) const {
            return this->channel != other.channel;
        }

        std::unique_ptr <T> operator*() {
            return std::move(element);
        }

    private:
        InputChannel *channel;
        std::unique_ptr <T> element;
    };


    template<class T>
    typename InputChannel<T>::Iterator begin(InputChannel<T> &channel) {
        return InputChannel<T>::Iterator(&channel);
    }

    template<class T>
    typename InputChannel<T>::Iterator end(InputChannel<T> &) {
        return InputChannel<T>::Iterator();
    }

    class OutputChannel::Iterator {
    public:
        Iterator(OutputChannel *c) : channel(c) {

        }

        template<class T>
        void operator=(T &data) {
            channel->push(std::make_unique<T>(data));
        }

        template<class T>
        void operator=(std::unique_ptr <T> &&ptr) {
            channel->push(ptr);
        }


        Iterator &operator++() {
            return *this;
        }

        Iterator &operator++(int) {
            return *this;
        }


        Iterator &operator*() {
            return *this;
        }

    private:
        OutputChannel *channel;


    private:
        Iterator *it;

    };

    OutputChannel::Iterator begin(OutputChannel &channel) {
        return OutputChannel::Iterator(&channel);
    }

    template<class T>
    void OutputChannel::push(std::unique_ptr <T> &&ptr) {
        this->push_message(std::unique_ptr<Message>(new TypedMessage<T>(ptr)));

    }


    template<class T>
    InputMessageChannel<T>::InputMessageChannel(std::shared_ptr <InputChannel<Message>> input,
                                                std::shared_ptr <Gadgetron::Core::OutputChannel> output):
            in(input), out(output) {}

    template<class T>
    std::unique_ptr <T> InputMessageChannel<T>::pop() {

        std::unique_ptr <Message> message = in->pop();

        while (typeid(*message) != typeid(T)) {
            out->push(std::move(message));
            message = in->pop();
        }

        return message;
    }

}
