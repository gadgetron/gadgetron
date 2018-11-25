#include "ConfigParsers.h"

#include "AST.h"
#include <iostream>
#include <memory>
#include <pugixml.hpp>
#include <boost/optional.hpp>
#include <boost/range.hpp>
#include <log.h>
#include <boost/parameter/name.hpp>

using namespace Gadgetron;
namespace {
    AST::Reader parse_reader(const pugi::xml_node &reader_node) {

        std::string port_str = reader_node.child_value("port");

        boost::optional<uint16_t> port = boost::none;
        if (!port_str.empty())
            port = static_cast<uint16_t>(std::stoi(port_str));

        return AST::Reader{reader_node.child_value("dll"),
                           reader_node.child_value("classname"),
                           port};
    }

    std::vector<AST::Reader> parse_readers(const pugi::xml_node &reader_root) {
        std::vector<AST::Reader> readers;
        for (const auto &node : reader_root.children("reader")) {
            readers.push_back(parse_reader(node));
        }
        return readers;
    }

    AST::Writer parse_writer(const pugi::xml_node &writer_node) {

        std::string port_str = writer_node.child_value("port");

        boost::optional<uint16_t> port = boost::none;
        if (!port_str.empty())
            port = static_cast<uint16_t>(std::stoi(port_str));

        return AST::Writer{writer_node.child_value("dll"),
                           writer_node.child_value("classname"),
                           port};
    }

    std::vector<AST::Writer> parse_writers(const pugi::xml_node &writer_root) {
        std::vector<AST::Writer> writers;
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

    AST::Gadget parse_gadget(const pugi::xml_node &gadget_node) {
        return AST::Gadget{gadget_node.child_value("name"),
                           gadget_node.child_value("dll"),
                           gadget_node.child_value("classname"),
                           parse_properties(gadget_node)};

    }




    namespace Legacy {


        std::vector<AST::Gadget> parse_gadgets(const pugi::xml_node& gadget_node){
            std::vector<AST::Gadget> gadgets;
            for (const auto& node : gadget_node.children("gadget")){
                gadgets.push_back(parse_gadget(node));
            }
            return gadgets;

        }

        AST::Stream parse_stream(const pugi::xml_node& stream_node){
            std::vector<AST::Node> nodes;
            boost::transform([](auto& a){return a;},parse_gadgets(stream_node),std::back_inserter(nodes));

            return AST::Stream{"main",nodes};
        }

        AST::Chain parse(const pugi::xml_document &config) {

            pugi::xml_node root = config.document_element().child("gadgetronStreamConfiguration");

            if (!root) {
                throw std::runtime_error("gadgetronStreamConfiguration element not found in configuration file");
            }

            return AST::Chain{parse_readers(root),parse_writers(root),parse_stream(root)};

        }

    }

    namespace V2 {

        AST::Stream parse_stream(const pugi::xml_node& stream_node ); //Forward declaration. Eww.

        //NOTE: Branchnode, mergenode and gadget are all kind of the same. Should it be the same code?
        // Conceptually they're very different.

        AST::MergeNode parse_mergenode(const pugi::xml_node& merge_node){
            return AST::MergeNode{merge_node.child_value("name"), merge_node.child_value("dll"),
                                  merge_node.child_value("classname"), parse_properties(merge_node)};
        }

        AST::BranchNode parse_branchnode(const pugi::xml_node& branch_node){
            return AST::BranchNode{branch_node.child_value("name"), branch_node.child_value("dll"),
                                  branch_node.child_value("classname"), parse_properties(branch_node)};
        }

        AST::Parallel parse_parallel(const pugi::xml_node& parallel_node){

            auto branchnode = parse_branchnode(parallel_node.child("branchnode"));
            auto mergenode = parse_mergenode(parallel_node.child("mergenode"));

            std::vector<AST::Stream> streams;
            for (const auto& stream_node : parallel_node.children("stream")){
                streams.push_back(parse_stream(stream_node));
            }
            return AST::Parallel{branchnode,mergenode,streams};
        }

        static const std::unordered_map<std::string,std::function<AST::Node(const pugi::xml_node&)>>
        node_parsers = {{"gadget",[](const pugi::xml_node& n){return parse_gadget(n);}},
                        {"parallel",[](const pugi::xml_node& n){return parse_parallel(n);}}};

        AST::Stream parse_stream(const pugi::xml_node& stream_node ){
            std::vector<AST::Node> nodes;
            for (auto& node : stream_node.children() ){
                nodes.push_back(node_parsers.at(node.name())(node));
            }
            return AST::Stream{stream_node.attribute("name").value(),nodes};
        }

        AST::Chain parse(const pugi::xml_document& config){

            auto root = config.child("configuration");

            auto readers = parse_readers(root.child("readers"));
            auto writers = parse_writers(root.child("writers"));
            auto stream = parse_stream(root.child("stream"));
            return AST::Chain{readers,writers,stream};
        }

    }

    std::function<AST::Chain(const pugi::xml_document&)> select_config_parser(const pugi::xml_document &raw_config) {

        if (raw_config.child("gadgetronStreamConfiguration")){
            return [](const pugi::xml_document& doc){return Legacy::parse(doc);};
        } else {
            return [](const pugi::xml_document &doc){return V2::parse(doc);};
        }

    }
}

namespace Gadgetron {

    AST::Chain parse_stream_configuration(std::istream &stream) {

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