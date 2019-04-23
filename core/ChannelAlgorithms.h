//
// Created by dchansen on 4/23/19.
//

#pragma once
#include "Channel.h"

namespace Gadgetron { namespace Core {

    template <class CHANNEL, class PRED> class TakeWhileChannel {

    public:
        TakeWhileChannel(CHANNEL& channel, PRED pred) : channel{ channel }, pred{ pred } {}

        auto pop() {
            if (done)
                throw ChannelClosed();
            auto message = channel.pop();
            done         = !pred(message);
            return std::move(message);
        }

    private:
        CHANNEL& channel;
        PRED pred;
        bool done = false;
    };

    template <class CHANNEL, class PRED>
    inline ChannelIterator<TakeWhileChannel<CHANNEL, PRED>> begin(TakeWhileChannel<CHANNEL, PRED>& channel) {
        return ChannelIterator<TakeWhileChannel<CHANNEL, PRED>>(&channel);
    }

    template <class CHANNEL, class PRED>
    inline ChannelIterator<TakeWhileChannel<CHANNEL, PRED>> end(TakeWhileChannel<CHANNEL, PRED>& channel) {
        return ChannelIterator<TakeWhileChannel<CHANNEL, PRED>>(&channel);
    }

    template <class CHANNEL, class PRED> auto take_while(CHANNEL& channel, PRED pred) {
        return TakeWhileChannel<CHANNEL, PRED>(channel, pred);
    }

}}
