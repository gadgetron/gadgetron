#pragma once

#include <typeinfo>

#include "ChannelIterator.h"

namespace Gadgetron { namespace Core {

    template <class CHANNEL> class ChannelRange {
    public:
        ChannelIterator<CHANNEL> begin() {
            return ChannelIterator<CHANNEL>(reinterpret_cast<CHANNEL*>(this));
        }

        ChannelIterator<CHANNEL> end() {
            return ChannelIterator<CHANNEL>();
        }
    };

    template <> class ChannelIterator<OutputChannel> {
    public:
        explicit ChannelIterator(OutputChannel* c) : channel(c) {}

        template <class T> void operator=(T&& data) {
            channel->push(std::forward<T>(data));
        }

        void operator=(Message&& message) {
            channel->push_message(std::move(message));
        }

        ChannelIterator& operator++() {
            return *this;
        }

        ChannelIterator& operator++(int) {
            return *this;
        }

        ChannelIterator& operator*() {
            return *this;
        }

    private:
        OutputChannel* channel;
    };

    inline ChannelIterator<OutputChannel> begin(OutputChannel& channel) {
        return ChannelIterator<OutputChannel>(&channel);
    }

    template <class... ARGS> inline void OutputChannel::push(ARGS&&... ptr) {
        this->push_message(Message(std::forward<ARGS>(ptr)...));
    }

} }
