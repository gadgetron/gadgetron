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
        using difference_type = std::ptrdiff_t;
        using value_type = void;
        using reference = void;
        using iterator_category = std::output_iterator_tag;

        explicit ChannelIterator(OutputChannel* c) : channel(c) {}
        ChannelIterator() = default;

        ChannelIterator(const ChannelIterator& other) = default;
        ChannelIterator(ChannelIterator&&) = default;
        ChannelIterator& operator=(const ChannelIterator&) = default;
        ChannelIterator& operator=(ChannelIterator&&) = default;

        template <class T> ChannelIterator& operator=(T&& data) {
            channel->push(std::forward<T>(data));
            return *this;
        }

        ChannelIterator& operator=(Message&& message) {
            channel->push_message(std::move(message));
            return *this;
        }

        ChannelIterator& operator++() {
            return *this;
        }

        ChannelIterator operator++(int) {
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

template<>
class std::iterator_traits<Gadgetron::Core::ChannelIterator<Gadgetron::Core::OutputChannel>> {
public:
    using difference_type = std::ptrdiff_t;
    using value_type = void;
    using reference = void;
    using iterator_category = std::output_iterator_tag;
};

