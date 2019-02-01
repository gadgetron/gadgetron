#pragma once

#include<typeinfo>


namespace Gadgetron::Core {

    template<class INPUTCHANNEL>
    class ChannelIterator {
    public:
        explicit ChannelIterator(INPUTCHANNEL *c) : channel(c) {
            this->operator++();
        }

        ChannelIterator() : channel(nullptr) {}

        ChannelIterator &operator++() {
            try {
                element = channel->pop();
            } catch (const ChannelClosed& err) {
                channel = nullptr;
            }
            return *this;
        };


        bool operator==(const ChannelIterator &other) const {
            return this->channel == other.channel;
        }

        bool operator!=(const ChannelIterator &other) const {
            return this->channel != other.channel;
        }

        auto operator*() {
            return std::move(element);
        }

    private:
        INPUTCHANNEL *channel;
        decltype(channel->pop()) element;
    };


    inline ChannelIterator<InputChannel> begin(InputChannel &channel) {
        return ChannelIterator<InputChannel>(&channel);
    }

    inline ChannelIterator<InputChannel> end(InputChannel &) {
        return ChannelIterator<InputChannel>();
    }

    template<>
    class ChannelIterator<OutputChannel> {
    public:
        explicit ChannelIterator(OutputChannel *c) : channel(c) {

        }

        template<class T>
        void operator=(T &&data) {
            channel->push(std::forward<T>(data));
        }

        void operator=(Message&& message) {
            channel->push_message(std::move(message));
        }


        ChannelIterator &operator++() {
            return *this;
        }

        ChannelIterator &operator++(int) {
            return *this;
        }


        ChannelIterator &operator*() {
            return *this;
        }

    private:
        OutputChannel *channel;



    };

    inline ChannelIterator<OutputChannel> begin(OutputChannel &channel) {
        return ChannelIterator<OutputChannel>(&channel);
    }

template<class ...ARGS>
inline void OutputChannel::push(ARGS&& ... ptr) {
    this->push_message(Message(std::forward<ARGS>(ptr)...));
}




    template<class... ARGS>
ChannelIterator<TypedInputChannel<ARGS...>> begin(
        TypedInputChannel<ARGS...> &channel) {
    return ChannelIterator < TypedInputChannel < ARGS...>>(&channel);
}

template<class... ARGS>
ChannelIterator <TypedInputChannel<ARGS...>> end(
        TypedInputChannel<ARGS...> &channel) {
    return ChannelIterator < TypedInputChannel < ARGS...>>();
}
}
