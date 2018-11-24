#pragma once

#include "Channel.h"
#include <thread>
#include <memory>
#include <future>
#include "log.h"
#include <ismrmrd/xml.h>
#include "PropteryMixin.h"

namespace Gadgetron::Core {

    class Node {
    public:
        virtual ~Node() {};

    };


    class GadgetNode : public Node, public std::enable_shared_from_this<GadgetNode>, public PropertyMixin {
    public:
        using Input = std::shared_ptr<InputChannel<Message>>;
        using Output = std::shared_ptr<OutputChannel>;

        GadgetNode(std::tuple<Input, Output> channels) : input_channel(std::get<0>(channels),
                                                                       output_channel(std::get<1>(channels)) {
            auto self = shared_from_this();
            std::thread([self]() { self->start(); }
            ).detach();
        }

        virtual ~GadgetNode() {};


    protected:

        virtual void process(std::shared_ptr<InputChannel<Message>> in, std::shared_ptr<OutputChannel> out) = 0;

    private:
        void start() {

            try {
                this->process(input_channel, output_channel);

            }
            catch (const ChannelClosedError &e) {
                output_channel->close();
            }
            catch (const std::exception &e) {
                GERROR(e.what());
            }
        }

        Input input_channel;
        Output output_channel;


    };

    template<class T>
    class TypedGadgetNode : public GadgetNode {
        TypedGadgetNode(const ISMRMRD::IsmrmrdHeader &header) {

        }

        virtual void process(std::shared_ptr<InputChannel<Message>> in, std::shared_ptr<OutputChannel> out) {
            auto typed_input = TypedInputChannel<T>(in, out);
            this->process(typed_input, *out);
        }


        virtual void process(InputChannel<T> &in, OutputChannel &out) = 0;

    };


    class LegacyGadgetNode : public GadgetNode {

    };


    class MergeNode : public Node {

    public:

        MergeNode(std::shared_future<std::shared_ptr<OutputChannel>> output_future,
                  const std::vector<std::string> &channel_names) {
            for (auto &name : channel_names) {
                input_channels.emplace({name, std::make_shared<MessageChannel>()});
            }
        }


    private:
        std::map<std::string, std::shared_ptr<MessageChannel>> input_channels;

    };
}



