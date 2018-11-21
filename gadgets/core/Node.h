#pragma once

#include "Channel.h"
#include <thread>
#include <memory>
#include <future>
#include "log.h"
#include <ismrmrd/xml.h>

namespace Gadgetron::Core {

    class Node {
    public:
        virtual ~Node() {};

    };

    class UnaryNode : public Node, public std::enable_shared_from_this<UnaryNode> {
    public:

        UnaryNode(std::shared_future<std::shared_ptr<OutputChannel>> output_future) : queue(
                std::make_shared<MessageChannel>()) {
            auto self = shared_from_this();
            std::thread([self, output_future]() { self->start(output_future); }
            ).detach();
        }

        virtual ~UnaryNode() {};

    protected:

        virtual void process(std::shared_ptr<InputChannel<Message>> in, std::shared_ptr<OutputChannel> out) = 0;

    private:
        void start(std::shared_future<std::shared_ptr<OutputChannel>> output_future) {

            try {
                auto output_channel = output_future.get();
                this->process(queue, output_channel);

            }
            catch (const ChannelClosedError &e) {

            }
            catch (const std::exception &e) {
                GERROR(e.what());
            }


        }

        std::shared_ptr<MessageChannel> queue;
    };

    template<class T>
    class GadgetNode : public UnaryNode {
        GadgetNode(const ISMRMRD::IsmrmrdHeader &header) {

        }

        virtual void process(std::shared_ptr<InputChannel<Message>> in, std::shared_ptr<OutputChannel> out) {
            auto typed_input = InputMessageChannel<T>(in, out);
            this->process(typed_input, *out);
        }


        virtual void process(InputChannel<T> &in, OutputChannel &out) = 0;

    };


    class LegacyGadgetNode : public UnaryNode {

    };


    class MergeNode : public Node {

    public:

        MergeNode(std::shared_future<std::shared_ptr<OutputChannel>> output_future, const std::vector<std::string>& channel_names ){
            for (auto& name : channel_names){
                input_channels.emplace({name, boost::make_shared<MessageChannel>()});
            }
        }


    private:
        std::map<std::string,std::shared_ptr<MessageChannel>>  input_channels;

    };
}



