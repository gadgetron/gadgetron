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

    protected:

        virtual void process(std::shared_ptr<InputChannel<Message>> in, std::shared_ptr<OutputChannel> out) = 0;




    };


    class GadgetNode : public Node, public PropertyMixin {
    public:
            GadgetNode(const ISMRMRD::IsmrmrdHeader& header, const std::unordered_map<std::string,std::string>& properties) : PropertyMixin(properties){};
            virtual ~GadgetNode(){};
    };

    template<class T>
    class TypedGadgetNode : public GadgetNode {
        TypedGadgetNode(const ISMRMRD::IsmrmrdHeader &header) {

        }

        virtual void process(std::shared_ptr<InputChannel<Message>> in, std::shared_ptr<OutputChannel> out) override final  {
            auto typed_input = TypedInputChannel<T>(in, out);
            this->process(typed_input, *out);
        }


        virtual void process(InputChannel<T> &in, OutputChannel &out) = 0;

    };


/*
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
    */
}



