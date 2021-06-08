#pragma once

namespace Gadgetron { namespace Core {

    template <class INPUTCHANNEL> class ChannelIterator {
    public:
        explicit ChannelIterator(INPUTCHANNEL* c) : channel(c) {
            this->operator++();
        }

        ChannelIterator() : channel(nullptr) {}
        ChannelIterator(const ChannelIterator&) = default;

    private:
        INPUTCHANNEL* channel;

    public:
        using difference_type   = long long;
        using value_type        = decltype(channel->pop());
        using pointer           = value_type*;
        using reference         = value_type&&;
        using const_reference = const value_type&;
        using iterator_category = std::input_iterator_tag;

    private:
        std::shared_ptr<value_type> element;

    public:
        ChannelIterator& operator++() {
            try {
                element = std::make_shared<value_type>(channel->pop());
            } catch (const ChannelClosed& err) {
                channel = nullptr;
            }
            return *this;
        };

        void operator++(int){
            this->operator++();
        }

        bool operator==(const ChannelIterator& other) const {
            return this->channel == other.channel;
        }

        bool operator!=(const ChannelIterator& other) const {
            return this->channel != other.channel;
        }

        reference operator*() {
            return std::move(*element);
        }

        const reference operator*() const {
            return *element;
        }

    };
} }
