//
// Created by dchansen on 4/23/19.
//

#pragma once

#include "Channel.h"

namespace Gadgetron {
    namespace Core {
        namespace Algorithm {

            namespace channel_detail {
                template<class CHANNEL, class PRED>
                class TakeWhileChannel {

                public:
                    TakeWhileChannel(CHANNEL &channel, PRED pred) : channel{channel}, pred{std::move(pred)} {}

                    auto pop() {
                        if (done)
                            throw ChannelClosed();
                        auto message = channel.pop();
                        done = !pred(message);
                        return std::move(message);
                    }

                private:
                    CHANNEL &channel;
                    PRED pred;
                    bool done = false;
                };

                template<class CHANNEL, class PRED>
                inline ChannelIterator<TakeWhileChannel<CHANNEL, PRED>>
                begin(TakeWhileChannel<CHANNEL, PRED> &channel) {
                    return ChannelIterator<TakeWhileChannel<CHANNEL, PRED>>(&channel);
                }

                template<class CHANNEL, class PRED>
                inline ChannelIterator<TakeWhileChannel<CHANNEL, PRED>> end(TakeWhileChannel<CHANNEL, PRED> &channel) {
                    return ChannelIterator<TakeWhileChannel<CHANNEL, PRED>>();
                }


                template<class CHANNEL, class PRED>
                class BufferChannel {
                public:
                    BufferChannel(Channel &channel, PRED pred) : channel{channel}, pred{std::move(pred)} {}

                    auto pop() {
                        if (closed) throw ChannelClosed();

                        std::vector<decltype(channel.pop())> buffer{};
                        buffer.push_back(channel.pop());
                        try {
                            while (!pred(buffer.back())) buffer.push_back(channel.pop());
                        } catch (ChannelClosed &channel) {
                            closed = true;
                        }
                        return buffer;

                    }

                private:
                    CHANNEL &channel;
                    PRED pred;
                    bool closed = false;

                };
            }

            /**
             * Takes a channel, and produces a channel-view which iterates over the channel until the
             * predicate returns false or the channel is closed. Note that it is inclusive, so the message returning false is included
             * @tparam CHANNEL A type which supports .pop and throws a ChannelClosed when closed. Typically an InputChannel
             * @tparam PRED
             * @param channel
             * @param pred
             * @return A new channel which closes when the predicate has been false
             */
            template<class CHANNEL, class PRED>
            auto take_while(CHANNEL &channel, PRED pred) {
                return channel_detail::TakeWhileChannel<CHANNEL, PRED>(channel, pred);
            }

            template<class CHANNEL, class PRED>
            auto take_until(CHANNEL &channel, PRED pred) {
                auto npred = [pred](auto m) { return !pred(m); };
                return take_while(channel, npred);
            }


            /**
             * Takes a channel and produces a channel-view which gathers messages into views, based on when the
             * predicate returns true (inclusive).
             * @tparam CHANNEL
             * @tparam PRED
             * @param channel
             * @param pred
             * @return
             */
            template<class CHANNEL, class PRED>
            auto buffer(CHANNEL &channel, PRED pred) {
                return channel_detail::BufferChannel<CHANNEL, PRED>(channel, pred);
            }

        }
    }
}
