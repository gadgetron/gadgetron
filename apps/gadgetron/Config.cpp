#include <pugixml.hpp>

#include <boost/optional.hpp>
#include <boost/parameter/name.hpp>
#include <boost/range/algorithm/transform.hpp>

#include "log.h"

#include "Config.h"

using namespace Gadgetron::Server;

namespace {

    Config::Reader parse_reader(const pugi::xml_node &reader_node) {

        std::string port_str = reader_node.child_value("port");

        boost::optional<uint16_t> port = boost::none;
        if (!port_str.empty())
            port = static_cast<uint16_t>(std::stoi(port_str));

        return Config::Reader{reader_node.child_value("dll"),
                              reader_node.child_value("classname"),
                              port};
    }

    std::vector<Config::Reader> parse_readers(const pugi::xml_node &reader_root) {
        std::vector<Config::Reader> readers;
        for (const auto &node : reader_root.children("reader")) {
            readers.push_back(parse_reader(node));
        }
        return readers;
    }

    Config::Writer parse_writer(const pugi::xml_node &writer_node) {
        return Config::Writer{writer_node.child_value("dll"),
                              writer_node.child_value("classname")};
    }

    std::vector<Config::Writer> parse_writers(const pugi::xml_node &writer_root) {
        std::vector<Config::Writer> writers;
        for (const auto &node : writer_root.children("reader")) {
            writers.push_back(parse_writer(node));
        }
        return writers;
    }


    std::unordered_map<std::string, std::string>
    parse_properties(const pugi::xml_node &root) {

        std::unordered_map<std::string, std::string> map;
        for (const auto &node : root.children("property")) {
            map.emplace(node.child_value("name"), node.child_value("value"));
        }

        return map;
    }

    Config::Gadget parse_gadget(const pugi::xml_node &gadget_node) {
        return Config::Gadget{gadget_node.child_value("name"),
                              gadget_node.child_value("dll"),
                              gadget_node.child_value("classname"),
                              parse_properties(gadget_node)};

    }

    namespace Legacy {

        std::vector<Config::Gadget> parse_gadgets(const pugi::xml_node& gadget_node){
            std::vector<Config::Gadget> gadgets;
            for (const auto& node : gadget_node.children("gadget")){
                gadgets.push_back(parse_gadget(node));
            }
            return gadgets;

        }

        Config::Stream parse_stream(const pugi::xml_node& stream_node){
            std::vector<Config::Node> nodes;
            boost::transform(parse_gadgets(stream_node),std::back_inserter(nodes),
                             [](auto gadget){return Config::Node(gadget);});

            return Config::Stream{"main",nodes};
        }

        Config parse(const pugi::xml_document &config) {

            pugi::xml_node root = config.child("gadgetronStreamConfiguration");

            if (!root) {
                throw std::runtime_error("gadgetronStreamConfiguration element not found in configuration file");
            }

            return Config{parse_readers(root),parse_writers(root),parse_stream(root)};
        }
    }

    namespace V2 {

        Config::Stream parse_stream(const pugi::xml_node& stream_node ); //Forward declaration. Eww.

        //NOTE: Branchnode, mergenode and gadget are all kind of the same. Should it be the same code?
        // Conceptually they're very different.

        Config::MergeNode parse_mergenode(const pugi::xml_node& merge_node){
            return Config::MergeNode{merge_node.child_value("name"), merge_node.child_value("dll"),
                                     merge_node.child_value("classname"), parse_properties(merge_node)};
        }

        Config::BranchNode parse_branchnode(const pugi::xml_node& branch_node){
            return Config::BranchNode{branch_node.child_value("name"), branch_node.child_value("dll"),
                                      branch_node.child_value("classname"), parse_properties(branch_node)};
        }

        Config::Parallel parse_parallel(const pugi::xml_node& parallel_node){

            auto branchnode = parse_branchnode(parallel_node.child("branchnode"));
            auto mergenode = parse_mergenode(parallel_node.child("mergenode"));

            std::vector<Config::Stream> streams;
            for (const auto& stream_node : parallel_node.children("stream")){
                streams.push_back(parse_stream(stream_node));
            }
            return Config::Parallel{branchnode,mergenode,streams};
        }

        static const std::unordered_map<std::string,std::function<Config::Node(const pugi::xml_node&)>>
                node_parsers = {{"gadget",[](const pugi::xml_node& n){return parse_gadget(n);}},
                                {"parallel",[](const pugi::xml_node& n){return parse_parallel(n);}}};

        Config::Stream parse_stream(const pugi::xml_node& stream_node ){
            std::vector<Config::Node> nodes;
            for (auto& node : stream_node.children() ){
                nodes.push_back(node_parsers.at(node.name())(node));
            }
            return Config::Stream{stream_node.attribute("name").value(),nodes};
        }

        Config parse(const pugi::xml_document& config){

            auto root = config.child("configuration");

            auto readers = parse_readers(root.child("readers"));
            auto writers = parse_writers(root.child("writers"));
            auto stream = parse_stream(root.child("stream"));
            return Config{readers, writers, stream};
        }

    }

    std::function<Config(const pugi::xml_document&)> select_config_parser(const pugi::xml_document &raw_config) {

        if (raw_config.child("gadgetronStreamConfiguration")){
            return [](const pugi::xml_document& doc){return Legacy::parse(doc);};
        } else {
            return [](const pugi::xml_document &doc){return V2::parse(doc);};
        }
    }
}

namespace Gadgetron::Server {

    Config parse_config(std::istream &stream) {

        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load(stream);

        if (result.status != pugi::status_ok) {
            GERROR("Loading config file failed with following error: %s (%d)\n", result.description(), result.status);
            throw std::runtime_error(result.description());
        }

        auto  parser = select_config_parser(doc);
        return parser(doc);
    }
}

