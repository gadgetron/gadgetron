#pragma once

namespace Gadgetron::Core {
    template<class ...ARGS>
    class InputChannel<ARGS...>::Iterator {
    public:
        Iterator(InputChannel<ARGS...> *c) : channel(*c) {
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

        auto operator*() {
            return std::move(element);
        }

    private:
        InputChannel *channel;
        decltype(channel->pop()) element;
    };


    template<class ...ARGS>
    typename InputChannel<ARGS...>::Iterator begin(InputChannel<ARGS...> &channel) {
        return InputChannel<ARGS...>::Iterator(&channel);
    }

    template<class ...ARGS>
    typename InputChannel<ARGS...>::Iterator end(InputChannel<ARGS...> &) {
        return InputChannel<ARGS...>::Iterator();
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

    template<class T, class U>
    inline void OutputChannel::push(std::unique_ptr<T> &&ptr) {
        this->push_message(std::make_unique<TypedMessage<T>>(std::move(ptr));
    }

    template<class T>
    inline void OutputChannel::push(std::unique_ptr<TypedMessage < T>>

    && ptr) {
    this->

    push_message (std::move(ptr));
}


template<class ...ARGS>
TypedInputChannel<ARGS...>::TypedInputChannel(
        std::shared_ptr<InputChannel < Message>>  input,
        std::shared_ptr<Gadgetron::Core::OutputChannel> output ):
in (input), bypass(output) {

    }

template<class T>
std::unique_ptr<T> TypedInputChannel<T>::pop() {

    std::unique_ptr<Message> message = in->pop();

    while (typeid(*message) != typeid(T)) {
        bypass->push(std::move(message));
        message = in->pop();
    }

    return std::unique_ptr<T>(static_cast<T *>(message.release()));
}

namespace {

    //Used for indices trick;
    template<int ...>
    struct seq {
    };
    template<int N, int ...S>
    struct gens : gens<N - 1, N - 1, S...> {
    };
    template<int ...S>
    struct gens<0, S...> {
        typedef seq<S...> type;
    };

    template<unsigned int I, class T, class ...REST>
    bool convertible_to_tuple(MessageTuple *messagetuple) {
        if (typeid(*messagetuple->messages()[I]) == typeid(TypedMessage < T > )) {
            return convertible_to_tuple<I + 1, REST...>(messagetuple);
        } else {
            return false;
        }
    }

    template<unsigned int I, class T>
    bool convertible_to_tuple(MessageTuple *messagetuple) {
        if (typeid(*messagetuple->messages()[I]) == typeid(TypedMessage < T > )) {
            return true;
        } else {
            return false;
        }
    }

    template<class T>
    std::unique_ptr<T> reinterpret_message(std::unique_ptr<Message> &&message) {
        return reinterpret_cast<TypedMessage <T> *>(message.get())->take_data();
    }

    template<class ...ARGS>
    struct tuple_maker {
        template<int ...S>
        auto from_messages(std::vector<std::unique_ptr<Message>> &messages, seq<S...>) {
            return std::make_tuple(reinterpret_message<ARGS>(messages[S])...);
        }
    };

    template<class ...ARGS>
    std::tuple<std::unique_ptr<ARGS>...> messageTuple_to_tuple(MessageTuple *messageTuple) {
        auto messages = messageTuple->take_messages();
        return tuple_maker<ARGS...>::make_tuple_from_messages(messages, typename gens<sizeof...(ARGS)>::type());
    }

    template<class ...ARGS>
    std::tuple<std::unique_ptr<ARGS>...> unpack(std::unique_ptr<Message> &message) {

        if (typeid(*message) == typeid(MessageTuple)) {
            auto *messagetuple = reinterpret_cast<MessageTuple *>(message.get());
            if (convertible_to_tuple<0, ARGS...>(messagetuple)) {
                message.release();
                return messageTuple_to_tuple(messagetuple);
            }
        }
        return std::tuple<std::unique_ptr<ARGS>...>();
    }
}


template<class ...ARGS>
std::tuple<std::unique_ptr<ARGS>...> TypedInputChannel<ARGS...>::pop() {

    while (true) {
        std::unique_ptr<Message> message = in->pop();
        auto result = unpack<ARGS...>(message);
        if (std::get<0>(result)) {
            return std::move(result);
        }
        bypass->push(std::move(message));
    }
}
}
