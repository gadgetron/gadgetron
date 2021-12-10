#include <pugixml.hpp>

#include <set>
#include <map>
#include <list>
#include <memory>
#include <string>

#include <boost/range/algorithm/transform.hpp>
#include <boost/range/algorithm/count_if.hpp>
#include <boost/range/algorithm/find_if.hpp>

#include "log.h"

#include "Config.h"
#include "Types.h"

using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Core;

namespace Gadgetron::Server::Connection {
    struct Config::External::Configuration {
        pugi::xml_document document;

        explicit Configuration(const pugi::xml_node &configuration_node) : document() {
            document.append_copy(configuration_node);
        }
    };
}

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

    template<class ConfigNode>
    constexpr const char *xml_name();

    template<>
    constexpr const char *xml_name<Config::Reader>() { return "reader"; }

    template<>
    constexpr const char *xml_name<Config::Writer>() { return "writer"; }

    template<>
    constexpr const char *xml_name<Config::Gadget>() { return "gadget"; }

    template<>
    constexpr const char *xml_name<Config::Branch>() { return "branch"; }

    template<>
    constexpr const char *xml_name<Config::Merge>() { return "merge"; }

    template<>
    constexpr const char *xml_name<Config::Distributor>() { return "distributor"; }

    struct XMLSerializer {

        template<class ConfigNode>
        static pugi::xml_node add_name(const ConfigNode &config, pugi::xml_node &node) {
            auto name = node.append_child("name");
            name.text().set(Config::name(config).c_str());
            return name;
        }

        template<class ConfigNode>
        static pugi::xml_node add_basenode(const ConfigNode &configNode, pugi::xml_node &node) {
            auto child_node = node.append_child(xml_name<ConfigNode>());
            auto dll = child_node.append_child("dll");
            dll.append_child(pugi::node_pcdata).set_value(configNode.dll.c_str());
            auto classname = child_node.append_child("classname");
            classname.append_child(pugi::node_pcdata).set_value(configNode.classname.c_str());
            return child_node;
        }

        static pugi::xml_node add_readers(const std::vector<Config::Reader> &readers, pugi::xml_node &node) {
            auto readers_node = node.append_child("readers");
            for (auto &reader : readers) add_basenode(reader, readers_node);
            return readers_node;
        }

        static pugi::xml_node add_writers(const std::vector<Config::Writer> &writers, pugi::xml_node &node) {
            auto writers_node = node.append_child("writers");
            for (auto &writer : writers) add_basenode(writer, writers_node);
            return writers_node;
        }

        static void add_property(const std::pair<std::string, std::string> &property, pugi::xml_node &node) {
            auto property_node = node.append_child("property");
            property_node.append_attribute("name").set_value(property.first.c_str());
            property_node.append_attribute("value").set_value(property.second.c_str());
        }

        template<class ConfigNode>
        static pugi::xml_node add_node(const ConfigNode &configNode, pugi::xml_node &node) {
            auto gadget_node = add_basenode(configNode, node);
            add_name(configNode, gadget_node);
            for (auto &property : configNode.properties) add_property(property, gadget_node);
            return gadget_node;
        }

        static pugi::xml_node add_node(const Config::Parallel &parallel, pugi::xml_node &node) {
            auto parallel_node = node.append_child("parallel");
            add_node(parallel.merge, parallel_node);
            add_node(parallel.branch, parallel_node);
            for (auto &stream : parallel.streams) {
                add_node(stream, parallel_node);
            }
            return parallel_node;
        }

        static pugi::xml_node add_node(const Config::Distributed &distributed, pugi::xml_node &node) {
            auto distributed_node = node.append_child("distributed");
            add_readers(distributed.readers, distributed_node);
            add_writers(distributed.writers, distributed_node);
            add_node(distributed.distributor, distributed_node);
            add_node(distributed.stream, distributed_node);

            return distributed_node;
        }

        static pugi::xml_node add_node(const Config::Execute &execute, pugi::xml_node &node) {
            auto execute_node = node.append_child("execute");
            execute_node.append_attribute("name").set_value(execute.name.c_str());
            execute_node.append_attribute("type").set_value(execute.type.c_str());
            return execute_node;
        }

        static pugi::xml_node add_node(const Config::Connect &connect, pugi::xml_node &node) {
            auto connect_node = node.append_child("connect");
            connect_node.append_attribute("port").set_value(connect.port.c_str());
            return connect_node;
        }

        static pugi::xml_node add_node(const Config::External &external, pugi::xml_node &node) {

            auto external_node = node.append_child("external");

            add_readers(external.readers, external_node);
            add_writers(external.writers, external_node);
            visit(
                    [&](auto action) { add_node(action, external_node); },
                    external.action
            );
            external_node.append_copy(external.configuration->document);

            return external_node;
        }

        static pugi::xml_node add_node(const Config::Stream &stream, pugi::xml_node &node) {
            auto stream_node = node.append_child("stream");
            stream_node.append_attribute("key").set_value(stream.key.c_str());
            for (auto n : stream.nodes) {
                visit([&stream_node](auto &typed_node) { add_node(typed_node, stream_node); }, n);
            }
            return stream_node;
        }

        static pugi::xml_node add_node(const Config::ParallelProcess& parallelProcess, pugi::xml_node & node){
            auto parallel_node = node.append_child("parallelprocess");
            parallel_node.append_attribute("workers").set_value((long long unsigned int)parallelProcess.workers);
            add_node(parallelProcess.stream, parallel_node);
            return parallel_node;
        }

        static pugi::xml_node add_node(const Config::PureStream& stream, pugi::xml_node& node){
            auto purestream_node = node.append_child("purestream");
            for (auto& gadget : stream.gadgets){
                add_node(gadget,purestream_node);
            }
            return purestream_node;
        }

        static pugi::xml_node add_node(const Config::PureDistributed& distributed, pugi::xml_node& node){
            auto puredistributed_node = node.append_child("puredistributed");
            add_readers(distributed.readers,puredistributed_node);
            add_writers(distributed.writers,puredistributed_node);
            add_node(distributed.stream,puredistributed_node);
            return puredistributed_node;
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
    }


    template<class Source>
    optional<Property> make_property(const pugi::xml_node &node) {
        if (!Source::accepts(node)) return none;
        return Property{Source::name(node), Source::value(node)};
    }


    template<class... Sources>
    Property parse_property(const pugi::xml_node &node) {
        std::vector<optional<Property>> potentials = {make_property<Sources>(node)...};
        auto to_bool = [](auto &potential) { return bool(potential); };

        auto n_valid = boost::count_if(potentials, to_bool);
        if (n_valid < 1) { throw ConfigNodeError("Unable to parse property", node); }
        if (n_valid > 1) { throw ConfigNodeError("Ambiguous property parse", node); }
        return **boost::find_if(potentials, to_bool);
    }


    template<class... Sources>
    class Parser {

    protected:

        explicit Parser(const pugi::xml_document &doc) : referenceable_properties(
                create_referenceable_properties(doc)) {}

        std::unordered_map<std::string, std::string>
        parse_properties(const pugi::xml_node &gadget_node) {

            std::unordered_map<std::string, std::string> properties;

            for (auto &node : gadget_node.children("property")) {
                auto property = parse_property<Sources...>(node);
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

            optional<uint16_t> slot = none;
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

        template<class NODE>
        NODE parse_node(const pugi::xml_node &gadget_node) {
            return NODE{gadget_node.child_value("name"),
                        gadget_node.child_value("dll"),
                        gadget_node.child_value("classname"),
                        parse_properties(gadget_node)};
        }

    private:
        PropertyMap referenceable_properties;

        static std::string value_of(const PropertyMap &map, const std::string &key, std::set<std::string> visited) {

            if (visited.count(key))
                throw std::runtime_error("Cyclic reference detected in Gadget xml property: " + key);

            auto value = map.at(key);

            if (is_reference(value)) {
                visited.insert(key);
                return value_of(map, value, visited);
            }

            return value;
        }

        static PropertyMap dereference_map(const PropertyMap &propertyMap) {
            PropertyMap resultMap = propertyMap;
            for (const auto &property : resultMap) {
                resultMap[property.first] = value_of(resultMap, property.first, std::set<std::string>());
            }
            return resultMap;
        }

        static PropertyMap create_referenceable_properties(const pugi::xml_node &root) {

            PropertyMap propertyMap;

            for (auto selector : root.select_nodes("//*[child::name and child::property]")) {
                auto node = selector.node();
                auto parent_name = node.child_value("name");

                for (auto p_node : node.children("property")) {
                    auto property = parse_property<Sources...>(p_node);
                    auto location = property.name + "@" + parent_name;

                    propertyMap[location] = property.value;
                }
            }

            return dereference_map(propertyMap);
        }
    };

    bool is_legacy_python_gadget(const Config::Gadget &gadget) {
        return gadget.dll == "gadgetron_python" && gadget.classname == "PythonGadget";
    }

    Config::Node transform_legacy_python_gadget(Config::Gadget gadget) {
        GDEBUG_STREAM("Legacy Python Gadget detected: " << gadget.name)

        pugi::xml_document document;
        auto configuration = document.append_child("configuration");

        for (auto &property : gadget.properties) {
            GDEBUG_STREAM("Appending property to configuration: " << property.first << ": " << property.second)
            XMLSerializer::add_property(property, configuration);
        }

        return Config::External{
            Config::Execute{
                    gadget.properties.at("python_module"),
                    "python",
                    gadget.properties.at("python_class")
            },
            std::make_shared<Config::External::Configuration>(configuration),
            std::vector<Config::Reader>(),
            std::vector<Config::Writer>()
        };
    }

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

        static bool accepts(const pugi::xml_document &config) {
            return config.child("gadgetronStreamConfiguration");
        }

    private:
        explicit Legacy(const pugi::xml_document &config) : Parser<LegacySource>(config) {}

        const std::list<std::pair<std::function<bool(const Config::Gadget &)>,
                                  std::function<Config::Node(Config::Gadget)>>> node_transformations{
            std::make_pair(is_legacy_python_gadget, transform_legacy_python_gadget),
            std::make_pair([](auto _) { return true; }, [=](auto c) { return Config::Node(c); })
        };

        Config::Node apply_transformation(Config::Gadget gadget) {

            auto pair = *std::find_if(
                    node_transformations.begin(),
                    node_transformations.end(),
                    [&](auto p) { return std::get<0>(p)(gadget); }
            );

            return std::get<1>(pair)(gadget);
        }

        std::vector<Config::Gadget> parse_gadgets(const pugi::xml_node &gadget_node) {
            std::vector<Config::Gadget> gadgets{};
            for (const auto &node : gadget_node.children("gadget")) {
                gadgets.push_back(parse_node<Config::Gadget>(node));
            }
            return gadgets;
        }

        Config::Stream parse_stream(const pugi::xml_node &stream_node) {
            std::vector<Config::Node> nodes;
            boost::transform(
                    parse_gadgets(stream_node),
                    std::back_inserter(nodes),
                    [&](auto gadget) { return apply_transformation(gadget); }
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

        static bool accepts(const pugi::xml_document &config) {
            auto configuration = config.child("configuration");
            return std::string(configuration.child_value("version")) == "2";
        }

    private:

        explicit V2(const pugi::xml_document &doc) : Parser<V2Source, LegacySource>(doc) {
            node_parsers["gadget"] = [&](const pugi::xml_node &n) { return this->parse_node<Config::Gadget>(n); };
            node_parsers["parallel"] = [&](const pugi::xml_node &n) { return this->parse_parallel(n); };
            node_parsers["external"] = [&](const pugi::xml_node &n) { return this->parse_external(n); };
            node_parsers["distributed"] = [&](const pugi::xml_node &n) { return this->parse_distributed(n); };
            node_parsers["parallelprocess"] = [&](const pugi::xml_node &n) { return this->parse_parallelprocess(n); };
            node_parsers["puredistributed"] = [&](const pugi::xml_node &n) { return this->parse_puredistributed(n); };
        }

        std::unordered_map<std::string, std::function<Config::Node(const pugi::xml_node &)>> node_parsers;

        Config::Parallel parse_parallel(const pugi::xml_node &parallel_node) {

            auto branch = parse_node<Config::Branch>(parallel_node.child("branch"));
            auto merge = parse_node<Config::Merge>(parallel_node.child("merge"));

            std::vector<Config::Stream> streams{};
            for (const auto &stream_node : parallel_node.children("stream")) {
                streams.push_back(parse_stream(stream_node));
            }
            return Config::Parallel{branch, merge, streams};
        }

        Config::External parse_external(const pugi::xml_node &external_node) {

            return Config::External{
                parse_action(external_node),
                parse_action_configuration(external_node),
                parse_readers(external_node.child("readers")),
                parse_writers(external_node.child("writers"))
            };
        }

        Config::Distributed parse_distributed(const pugi::xml_node &distributed_node) {
            auto distributor = parse_node<Config::Distributor>(distributed_node.child("distributor"));
            auto stream = parse_stream(distributed_node.child("stream"));
            auto readers = parse_readers(distributed_node.child("readers"));
            auto writers = parse_writers(distributed_node.child("writers"));
            return {readers,writers,distributor,stream};
        }

        Config::Stream parse_stream(const pugi::xml_node &stream_node) {
            std::vector<Config::Node> nodes;
            for (auto &node : stream_node.children()) {
                nodes.push_back(node_parsers.at(node.name())(node));
            }
            return Config::Stream{stream_node.attribute("key").value(), nodes};
        }

        Config::PureStream parse_purestream(const pugi::xml_node &purestream_node){
            std::vector<Config::Gadget> gadgets;
            boost::transform(purestream_node.children(), std::back_inserter(gadgets),
                [&,this](auto node) { return this->parse_node<Config::Gadget>(node); });
            return {gadgets};
        }

        Config::ParallelProcess parse_parallelprocess(const pugi::xml_node& parallelprocess_node)
        {
            size_t workers = std::stoul(parallelprocess_node.attribute("workers").value());
            return Config::ParallelProcess{workers,parse_purestream(parallelprocess_node.child("purestream"))};
        }

        Config::PureDistributed parse_puredistributed(const pugi::xml_node& puredistributedprocess_node){

            auto purestream = parse_purestream(puredistributedprocess_node.child("purestream"));
            auto readers = parse_readers(puredistributedprocess_node.child("readers"));
            auto writers = parse_writers(puredistributedprocess_node.child("writers"));
            return {readers,writers,purestream};
        }

        static optional<std::string> parse_target(std::string s) {
            if (s.empty()) return none;
            return s;
        }

        static Config::Execute parse_execute(const pugi::xml_node &execute_node) {
            return Config::Execute {
                execute_node.attribute("name").value(),
                execute_node.attribute("type").value(),
                parse_target(execute_node.attribute("target").value())
            };
        }

        static std::string address_or_localhost(const std::string &s) {
            return s.empty() ? "localhost" : s;
        }

        static Config::Connect parse_connect(const pugi::xml_node &connect_node) {
            return Config::Connect {
                address_or_localhost(connect_node.attribute("address").value()),
                connect_node.attribute("port").value()
            };
        }

        std::map<std::string, std::function<Config::Action(const pugi::xml_node &)>> action_parsers {
                {"execute", parse_execute},
                {"connect", parse_connect}
        };

        Config::Action parse_action(const pugi::xml_node &external_node) {

            for (pugi::xml_node child : external_node) {
                if (action_parsers.count(child.name())) {
                    return action_parsers.at(child.name())(child);
                }
            }

            throw ConfigNodeError("Unable to parse valid action for external node", external_node);
        }

        std::shared_ptr<Config::External::Configuration>
        parse_action_configuration(const pugi::xml_node &external_node) {
            return std::make_shared<Config::External::Configuration>(external_node.child("configuration"));
        }
    };
}

namespace Gadgetron::Server::Connection {

    Config parse_config(std::istream &stream) {

        auto parsers = {
                std::make_pair(Legacy::accepts, Legacy::parse),
                std::make_pair(V2::accepts, V2::parse)
        };

        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load(stream);

        if (result.status != pugi::status_ok) {
            GERROR("Loading config file failed with following error: %s (%d)\n", result.description(), result.status);
            throw std::runtime_error(result.description());
        }

        auto parser = std::find_if(parsers.begin(), parsers.end(), [&](auto pair) { return std::get<0>(pair)(doc); });
        if (parser == parsers.end()) throw std::runtime_error("No parser accepted config file.");
        return std::get<1>(*parser)(doc);
    }

    std::string serialize_config(const Config &config) {
        pugi::xml_document doc{};
        auto config_node = doc.append_child("configuration");
        config_node.append_child("version").text().set(2);
        XMLSerializer::add_readers(config.readers, config_node);
        XMLSerializer::add_writers(config.writers, config_node);
        XMLSerializer::add_node(config.stream, config_node);

        std::stringstream stream;
        doc.save(stream);
        return stream.str();
    }

    std::string serialize_config(const Config::External& external_config) {
        std::stringstream stream;
        external_config.configuration->document.save(stream);
        return stream.str();
    }
}
