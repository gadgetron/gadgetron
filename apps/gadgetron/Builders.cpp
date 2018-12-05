

#include "log.h"

#include "Builders.h"

namespace Gadgetron::Server::Builders {

    using namespace Gadgetron::Core;
    using namespace boost::dll;

    Builder::Builder(const Config &config, const Context::Paths &paths) : config(config), paths(paths) {}

    boost::dll::shared_library Builder::load_library(const std::string &shared_library_name) {
        return boost::dll::shared_library(
                this->make_library_path(shared_library_name),
                boost::dll::load_mode::search_system_folders |
                boost::dll::load_mode::append_decorations
        );
    }

    boost::filesystem::path Builder::make_library_path(const std::string &shared_library_name) {
        return paths.gadgetron_home / "lib" / ("lib" +shared_library_name+ ".so");
    }

    ReaderBuilder::ReaderBuilder(const Config &config, const Context::Paths &paths) : Builder(config, paths) {}

    void ReaderBuilder::process(
            std::function<void(uint16_t, std::unique_ptr<Reader>, boost::dll::shared_library)> on_reader) {

        for (auto &reader_config : config.readers) {

            auto library = this->load_library(reader_config.dll);
            auto factory = library.get_alias<std::unique_ptr<Reader>(void)>(
                    "reader_factory_export_" + reader_config.classname);
            auto reader = factory();

            uint16_t port = reader_config.port.value_or(reader->port());

            on_reader(port, std::move(reader), library);
        }
    }

    WriterBuilder::WriterBuilder(const Config &config, const Context::Paths &paths) : Builder(config, paths) {}

    void WriterBuilder::process(std::function<void(std::unique_ptr<Writer>, boost::dll::shared_library)> on_reader) {
        throw std::runtime_error("Not implemented yet.");
    }

    StreamBuilder::StreamBuilder(
            const Config &config,
            const Context &context,
            std::shared_ptr<MessageChannel> input,
            std::shared_ptr<MessageChannel> output
    ) : Builder(config, context.paths), context(context), input(std::move(input)), output(std::move(output)) {}

    void StreamBuilder::process(std::function<void(std::unique_ptr<Node>, boost::dll::shared_library)> on_node) {

        for (auto &node_config : config.stream.nodes) {
            boost::apply_visitor(
                    [&](auto &node) {
                        this->load_node(node, on_node);
                    },
                    node_config
            );
        }
    }

    using gadget_factory = std::unique_ptr<Node>(const Context &, const std::unordered_map<std::string, std::string> &);

    void StreamBuilder::load_node(const Config::Gadget &gadget_config, node_callback on_node) {
        GDEBUG_STREAM("Loading Gadget node from: " << gadget_config.dll << std::endl);

        auto library = this->load_library(gadget_config.dll);
        auto factory = library.get_alias<gadget_factory>("gadget_factory_export_" + gadget_config.classname);

        auto node = factory(context, gadget_config.properties);

        on_node(std::move(node), library);
    }

    void StreamBuilder::load_node(const Config::Parallel &parallel_config, node_callback on_node) {
        GDEBUG_STREAM("Loading Parallel node." << std::endl);
    }

    void StreamBuilder::load_node(const Config::Distributed &parallel_config, node_callback on_node) {
        GDEBUG_STREAM("Loading Distributed node." << std::endl);
    }
}