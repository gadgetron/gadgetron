#include "ChainBuilder.h"
#include "Gadget.h"
#include "AST.h"
#include "Node.h"
#include "Channel.h"
#include "Stream.h"
#include <ismrmrd/xml.h>

#include <boost/dll.hpp>
#include <boost/variant.hpp>


#include <memory>

namespace Gadgetron::Core {

    class StreamHandler { //TODO:: Name? Is this in fact a Chain? We need to keep state here, as the dlls will unload
                          // and calling Gadgets with unloaded dlls will yield undefined behaviour (demon bats out of your nose variety)
    private:
        using factory_type = std::shared_ptr<Node>(
                std::tuple<std::shared_ptr<InputChannel<Message>>,
                std::shared_ptr<OutputChannel>>,
        const ISMRMRD::IsmrmrdHeader header,
        const std::unordered_map<std::string, std::string>& parameters );

        using Channels = std::tuple<std::shared_ptr<InputChannel<Message>>, std::shared_ptr<OutputChannel>>;

        std::shared_ptr<Node> make_node(const AST::Parallel &parallel_ast,
                                        Channels channels) {
            //TODO: we might want a different implementation here
            throw std::runtime_error("Not implementated yet");
        }

        std::shared_ptr<Node> make_node(const AST::Gadget &gadget_ast,
                                          Channels channels) {
            auto library = boost::dll::shared_library(
                    gadget_ast.dll, boost::dll::load_mode::append_decorations);
            auto creator = library.get<factory_type>("make_" + gadget_ast.classname);
            libraries.push_back(library);

            return creator(channels, header, gadget_ast.properties);
        }

        std::shared_ptr<Stream> make_node(const AST::Stream &stream_ast,
                                          std::shared_ptr<InputChannel < Message>> channel) {
            auto input_channel = channel;
            auto output_channel = std::make_shared<MessageChannel>();


            for (auto &node_ast : stream_ast.nodes) {
                auto node = boost::apply_visitor(
                        [this,&input_channel,&output_channel](auto &ast) {
                            return this->make_node(ast, Channels{input_channel,output_channel}); },
                        node_ast);
                input_channel = output_channel;
                output_channel = std::make_shared<MessageChannel>();
            }

            return std::make_shared<Stream>(stream_ast.name, input_channel);
        }

        const ISMRMRD::IsmrmrdHeader header;
        std::vector<boost::dll::shared_library> libraries;
    };

}  // namespace Gadgetron::Core
