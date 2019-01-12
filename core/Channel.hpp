#pragma once

#include<typeinfo>

namespace Gadgetron::Core {

    template<class INPUTCHANNEL>
    class ChannelIterator {
    public:
        ChannelIterator(INPUTCHANNEL *c) : channel(c) {
            this->operator++();
        }

        ChannelIterator() : channel(nullptr) {}

        ChannelIterator &operator++() {
            try {
                element = channel->pop();
            } catch (ChannelClosed err) {
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
        ChannelIterator(OutputChannel *c) : channel(c) {

        }

        template<class T>
        void operator=(T &data) {
            channel->push(std::make_unique<T>(data));
        }

        template<class T>
        void operator=(std::unique_ptr<T> &&ptr) {
            channel->push(ptr);
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


    private:
        ChannelIterator *it;

    };

    inline ChannelIterator<OutputChannel> begin(OutputChannel &channel) {
        return ChannelIterator<OutputChannel>(&channel);
    }

    namespace {
        namespace gadgetron_detail {

            template<class T>
            std::unique_ptr<TypedMessage < std::remove_reference_t<T>>>
            make_message(T && data) {
            return std::make_unique<TypedMessage < std::remove_reference_t<T>>>(
            std::forward<T>(data)
            );
        }

        template<class... ARGS>
        std::enable_if_t<(sizeof...(ARGS) > 1), std::unique_ptr<MessageTuple>>
        make_message(ARGS &&... data) {
            return std::make_unique<MessageTuple>(std::forward<ARGS>(data)...);
        }

    }  // namespace gadgetron_detail
}  // namespace

template<class ...ARGS>
inline void OutputChannel::push(ARGS&& ... ptr) {
    this->push_message(gadgetron_detail::make_message<ARGS...>(std::forward<ARGS>(ptr)...));

}


template<class ...ARGS>
inline void OutputChannel::push(std::unique_ptr<ARGS>&& ... ptr) {
        penguin(ptr...);

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
