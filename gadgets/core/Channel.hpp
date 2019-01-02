#pragma once

#include<typeinfo>

namespace Gadgetron::Core {

    class InputChannel::Iterator {
    public:
        Iterator(InputChannel *c) : channel(c) {
            this->operator++();
        }

        Iterator() : channel(nullptr) {}

        Iterator &operator++() {
            try {
                element = channel->pop();
            } catch (ChannelClosed err) {
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

        auto operator*() {
            return std::move(element);
        }

    private:
        InputChannel *channel;
        decltype(channel->pop()) element;
    };





    inline typename InputChannel::Iterator begin(InputChannel &channel) {
        return typename InputChannel::Iterator(&channel);
    }

    inline typename InputChannel::Iterator end(InputChannel &) {
        return typename InputChannel::Iterator();
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
        void operator=(std::unique_ptr<T> &&ptr) {
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

    inline OutputChannel::Iterator begin(OutputChannel &channel) {
        return OutputChannel::Iterator(&channel);
    }


    namespace {
        namespace gadgetron_detail {

        template<class T>
        std::unique_ptr<TypedMessage<T>>
        make_message(std::unique_ptr<T> && ptr) {
                return std::make_unique<TypedMessage < T>>(
                        std::move(ptr)
                        );
        }

        template<class ...ARGS>
        std::enable_if_t<(sizeof...(ARGS) > 1), std::unique_ptr<MessageTuple>>
        make_message(std::unique_ptr<ARGS> &&... ptrs) {
            return std::make_unique<MessageTuple>(ptrs...);

        }

    }
}

template<class ...ARGS>
inline void OutputChannel::push(std::unique_ptr<ARGS> &&... ptr) {
    this->push_message(gadgetron_detail::make_message<ARGS...>(std::move(ptr)...));

}
/*

template<class ...ARGS>
TypedInputChannel<ARGS...>::TypedInputChannel(
        std::shared_ptr<InputChannel < Message>>

input,
std::shared_ptr<Gadgetron::Core::OutputChannel> output
):

in (input), bypass(output) {

}

template<class T>
std::unique_ptr<T> TypedInputChannel<T>::pop() {

    std::unique_ptr<Message> message = in->pop();
    while (!convertible_to<T>(message)) {
        bypass->push(std::move(message));
        message = in->pop();
    }

    return force_unpack<T>(message);
}


template<class ...ARGS>
std::tuple<std::unique_ptr<ARGS>...> TypedInputChannel<ARGS...>::pop() {
    std::unique_ptr<Message> message = in->pop();
    while (!convertible_to<ARGS...>(message)) {
        bypass->push(std::move(message));
        message = in->pop();
    }

    return force_unpack<ARGS...>(message);

}
 */
}
