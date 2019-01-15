#include <pugixml.hpp>

#include <set>
#include <map>
#include <memory>
#include <string>

#include <boost/optional.hpp>
#include <boost/parameter/name.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/range/algorithm/find_if.hpp>
#include <numeric>

#include "log.h"

#include "Config.h"
#include "Types.h"

using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Core;

namespace {


    class ConfigNodeError : public std::runtime_error {
    public:
        ConfigNodeError(const std::string &message, const pugi::xml_node &node) : std::runtime_error(
                make_message(message, node)) {}

    private:
        static std::string make_message(const std::string &message, const pugi::xml_node &node) {
            std::stringstream stream;
            stream << message << " ";
            node.print(stream, "", pugi::format_raw);
            return stream.str();
        }
    };

    struct Property {
        std::string name;
        std::string value;
    };
    using Location = std::string;
    using Value = std::string;
    using PropertyMap = std::map<Location, Value>;


    class LegacySource {
    public:
        static std::string name(const pugi::xml_node &node) {
            return node.child_value("name");
        }

        static std::string value(const pugi::xml_node &node) {
            return node.child_value("value");
        }

        static bool accepts(const pugi::xml_node &node) {
            return node.child("name") && node.child("value");
        }
    };

    class V2Source {
    public:

        static std::string name(const pugi::xml_node &node) {
            return node.attribute("name").value();
        }

        static std::string value(const pugi::xml_node &node) {
            return node.attribute("value").value();
        }

        static bool accepts(const pugi::xml_node &node) {
            return node.attribute("name") && node.attribute("value");
        }
    };


    bool is_reference(const std::string &value) {
        return value.find('@') != std::string::npos;
    };


    template<class Source>
    optional<Property> make_property(const pugi::xml_node &node) {
        if (!Source::accepts(node)) return boost::none;
        return Property{Source::name(node), Source::value(node)};
    }


    template<class... Rulesets>
    Property parse_property(const pugi::xml_node &node) {
        std::vector<optional<Property>> potentials = {make_property<Rulesets>(node)...};
        auto to_bool = [](auto& potential) {return bool(potential);};

        auto n_valid = boost::count_if(potentials, to_bool );
        if (n_valid < 1) { throw ConfigNodeError("Unable to parse property", node); };
        if (n_valid > 1) { throw ConfigNodeError("Ambigous property parse", node); };
        return **boost::find_if(potentials,to_bool);
    }


    template<class... Rulesets>
    class Parser {

    protected:

        Parser(const pugi::xml_document &doc) : referenceable_properties(create_referenceable_properties(doc)) {

        }

        std::unordered_map<std::string, std::string>
        parse_properties(const pugi::xml_node &gadget_node) {

            std::unordered_map<std::string, std::string> properties;

            for (auto &node : gadget_node.children("property")) {
                auto property = parse_property<Rulesets...>(node);
                properties[property.name] = dereference(property.value);
            }

            return properties;
        }

        std::string dereference(const std::string &value_string) {
            if (is_reference(value_string)) return referenceable_properties[value_string];
            return value_string;
        }

        static Config::Reader parse_reader(const pugi::xml_node &reader_node) {

            std::string slot_str = reader_node.child_value("slot");

            boost::optional<uint16_t> slot = boost::none;
            if (!slot_str.empty())
                slot = static_cast<uint16_t>(std::stoi(slot_str));

            return Config::Reader{reader_node.child_value("dll"),
                                  reader_node.child_value("classname"),
                                  slot};
        }

        static std::vector<Config::Reader> parse_readers(const pugi::xml_node &reader_root) {
            std::vector<Config::Reader> readers{};
            for (const auto &node : reader_root.children("reader")) {
                readers.push_back(parse_reader(node));
            }
            return readers;
        }

        static Config::Writer parse_writer(const pugi::xml_node &writer_node) {
            return Config::Writer{writer_node.child_value("dll"),
                                  writer_node.child_value("classname")};
        }

        static std::vector<Config::Writer> parse_writers(const pugi::xml_node &writer_root) {
            std::vector<Config::Writer> writers{};
            for (const auto &node : writer_root.children("writer")) {
                writers.push_back(parse_writer(node));
            }
            return writers;
        }

        Config::Gadget parse_gadget(const pugi::xml_node &gadget_node) {
            return Config::Gadget{gadget_node.child_value("name"),
                                  gadget_node.child_value("dll"),
                                  gadget_node.child_value("classname"),
                                  parse_properties(gadget_node)};
        }

    private:
        PropertyMap referenceable_properties;

        static std::string dereference_key(PropertyMap &map, const std::string &key) {
            try {
                std::string value = map.at(key);
                if (!is_reference(value)) return value;
                map.erase(key);
                auto val = dereference_key(map, value);
                map[key] = val;
                return val;
            } catch (const std::out_of_range&) {
                throw std::runtime_error("Cycle detected in Gadget xml properties");
            }
        }

        static PropertyMap dereference_map(const PropertyMap &propertyMap) {
            PropertyMap resultMap = propertyMap;
            for (auto property : resultMap) {
                resultMap[property.first] = dereference_key(resultMap, property.first);
            }
            return resultMap;
        }

        static PropertyMap create_referenceable_properties(const pugi::xml_node &root) {

            PropertyMap propertyMap;

            for (auto selector : root.select_nodes("//*[child::name and child::property]")) {
                auto node = selector.node();
                auto parent_name = node.child_value("name");

                for (auto p_node : node.children("property")) {
                    auto property = parse_property<Rulesets...>(p_node);
                    auto location = property.name + "@" + parent_name;

                    propertyMap[location] = property.value;
                }
            }

            return dereference_map(propertyMap);
        }
    };

    class Legacy : public Parser<LegacySource> {
    public:

        static Config parse(const pugi::xml_document &config) {

            auto parser = Legacy(config);
            pugi::xml_node root = config.child("gadgetronStreamConfiguration");

            return Config{
                    parse_readers(root),
                    parse_writers(root),
                    parser.parse_stream(root)
            };
        }


    private:

        Legacy(const pugi::xml_document &config) : Parser<LegacySource>(config) {}

        std::vector<Config::Gadget> parse_gadgets(const pugi::xml_node &gadget_node) {
            std::vector<Config::Gadget> gadgets{};
            for (const auto &node : gadget_node.children("gadget")) {
                gadgets.push_back(parse_gadget(node));
            }
            return gadgets;
        }

        Config::Stream parse_stream(const pugi::xml_node &stream_node) {
            std::vector<Config::Node> nodes;
            boost::transform(
                    parse_gadgets(stream_node),
                    std::back_inserter(nodes),
                    [](auto gadget) {
                        return Config::Node(gadget);
                    }
            );

            return Config::Stream{"main", nodes};
        }
    };


    class V2 : public Parser<V2Source, LegacySource> {
    public:


        static Config parse(const pugi::xml_document &config) {

            auto parser = V2(config);
            auto root = config.child("configuration");

            return Config{
                    parse_readers(root.child("readers")),
                    parse_writers(root.child("writers")),
                    parser.parse_stream(root.child("stream"))
            };
        }

    private:

        V2(const pugi::xml_document &doc) : Parser<V2Source, LegacySource>(doc) {
            node_parsers["gadget"] = [&](const pugi::xml_node &n) { return this->parse_gadget(n); };
            node_parsers["parallel"] = [&](const pugi::xml_node &n) { return this->parse_parallel(n); };
            node_parsers["distributed"] = [&](const pugi::xml_node &n) { return this->parse_distributed(n); };

        }

        std::unordered_map<std::string, std::function<Config::Node(const pugi::xml_node &)>> node_parsers;

        Config::Merge parse_mergenode(const pugi::xml_node &merge_node) {
            return Config::Merge{merge_node.child_value("name"), merge_node.child_value("dll"),
                                 merge_node.child_value("classname"), parse_properties(merge_node)};
        }

        Config::Branch parse_branchnode(const pugi::xml_node &branch_node) {
            return Config::Branch{branch_node.child_value("name"), branch_node.child_value("dll"),
                                  branch_node.child_value("classname"), parse_properties(branch_node)};
        }

        Config::Parallel parse_parallel(const pugi::xml_node &parallel_node) {

            auto branch = parse_branchnode(parallel_node.child("branch"));
            auto merge = parse_mergenode(parallel_node.child("merge"));

            std::vector<Config::Stream> streams{};
            for (const auto &stream_node : parallel_node.children("stream")) {
                streams.push_back(parse_stream(stream_node));
            }
            return Config::Parallel{branch, merge, streams};
        }

        Config::Distributed parse_distributed(const pugi::xml_node& distributed_node){
            auto branch = parse_branchnode(distributed_node.child("branch"));
            auto merge = parse_mergenode(distributed_node.child("merge"));

            std::vector<Config::Stream> streams{};
            for (const auto &stream_node : distributed_node.children("stream")) {
                streams.push_back(parse_stream(stream_node));
            }

            throw std::runtime_error("Not really implemented yet");
        }

        Config::Stream parse_stream(const pugi::xml_node &stream_node) {
            std::vector<Config::Node> nodes;
            for (auto &node : stream_node.children()) {
                nodes.push_back(node_parsers.at(node.name())(node));
            }
            return Config::Stream{stream_node.attribute("key").value(), nodes};
        }
    };
}

namespace Gadgetron::Server::Connection {

    Config parse_config(std::istream &stream) {

        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load(stream);

        if (result.status != pugi::status_ok) {
            GERROR("Loading config file failed with following error: %s (%d)\n", result.description(), result.status);
            throw std::runtime_error(result.description());
        }

        if (doc.child("gadgetronStreamConfiguration"))
            return Legacy::parse(doc);

        return V2::parse(doc);
    }
}
