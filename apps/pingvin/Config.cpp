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

using namespace Gadgetron::Main;

namespace {

    class ConfigNodeError : public std::runtime_error {
    public:
        ConfigNodeError(const std::string &message, const pugi::xml_node &node) : std::runtime_error(
                make_message(message, node)) {}

    private:
        static std::string make_message(const std::string &message, const pugi::xml_node &node) {
            std::stringstream stream;
            stream << message << "\n";
            node.print(stream);
            return stream.str();
        }
    };

    template<class ConfigNode>
    constexpr const char *xml_name();

    template<>
    constexpr const char *xml_name<Config::Gadget>() { return "gadget"; }

    template<>
    constexpr const char *xml_name<Config::Branch>() { return "branch"; }

    template<>
    constexpr const char *xml_name<Config::Merge>() { return "merge"; }

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
    std::optional<Property> make_property(const pugi::xml_node &node) {
        if (!Source::accepts(node)) return std::nullopt;
        return Property{Source::name(node), Source::value(node)};
    }


    template<class... Sources>
    Property parse_property(const pugi::xml_node &node) {
        std::vector<std::optional<Property>> potentials = {make_property<Sources>(node)...};
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

        template<class NODE>
        NODE parse_node(const pugi::xml_node &gadget_node) {
            std::string name(gadget_node.child_value("name"));
            std::string dll(gadget_node.child_value("dll"));
            std::string classname(gadget_node.child_value("classname"));

            if (dll.empty()) {
                throw ConfigNodeError("Missing dll for gadget " + name, gadget_node);
            }
            if (classname.empty()) {
                throw ConfigNodeError("Missing classname for gadget " + name, gadget_node);
            }

            return NODE{name, dll, classname, parse_properties(gadget_node)};
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

    class Legacy : public Parser<LegacySource> {
    public:

        static Config parse(const pugi::xml_document &config) {

            auto parser = Legacy(config);
            pugi::xml_node root = config.child("gadgetronStreamConfiguration");

            return Config{ parser.parse_stream(root) };
        }

        static bool accepts(const pugi::xml_document &config) {
            return config.child("gadgetronStreamConfiguration");
        }

    private:
        explicit Legacy(const pugi::xml_document &config) : Parser<LegacySource>(config) {}

        const std::list<std::pair<std::function<bool(const Config::Gadget &)>,
                                  std::function<Config::Node(Config::Gadget)>>> node_transformations{
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

            return Config{ parser.parse_stream(root.child("stream")) };
        }

        static bool accepts(const pugi::xml_document &config) {
            auto configuration = config.child("configuration");
            return std::string(configuration.child_value("version")) == "2";
        }

    private:

        explicit V2(const pugi::xml_document &doc) : Parser<V2Source, LegacySource>(doc) {
            node_parsers["gadget"] = [&](const pugi::xml_node &n) { return this->parse_node<Config::Gadget>(n); };
            node_parsers["parallel"] = [&](const pugi::xml_node &n) { return this->parse_parallel(n); };
            node_parsers["parallelprocess"] = [&](const pugi::xml_node &n) { return this->parse_parallelprocess(n); };
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

    };
}

namespace Gadgetron::Main {

    Config Config::parse(std::istream &stream) {

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

    std::string Config::serialize(const Config &config) {
        pugi::xml_document doc{};
        auto config_node = doc.append_child("configuration");
        config_node.append_child("version").text().set(2);
        XMLSerializer::add_node(config.stream, config_node);

        std::stringstream stream;
        doc.save(stream);
        return stream.str();
    }
}
