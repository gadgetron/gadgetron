#include <pugixml.hpp>

#include <set>
#include <map>
#include <memory>
#include <string>

#include <boost/optional.hpp>
#include <boost/parameter/name.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <numeric>

#include "log.h"

#include "Config.h"

using namespace Gadgetron::Server::Connection;

namespace {

    void throw_node_error(const std::string &message, const pugi::xml_node &node) {

        std::stringstream stream;

        stream << message << " ";
        node.print(stream, "", pugi::format_raw);

        throw std::runtime_error(stream.str());
    }

    class Property {
    public:
        virtual std::string name() = 0;
        virtual std::string value() = 0;
        virtual std::string value(std::set<std::string> &visited) {
            return value();
        };
        virtual ~Property() = default;
    };

    class PropertyBuilder {
    public:
        virtual bool accepts(const pugi::xml_node &node) const = 0;
        virtual std::unique_ptr<Property> build(const pugi::xml_node &node) const = 0;
        virtual ~PropertyBuilder() = default;
    };

    using Location = std::string;
    using PropertyMap = std::map<Location, std::unique_ptr<Property>>;
    using PropertyBuilders = std::vector<std::unique_ptr<PropertyBuilder>>;

    class Source {
    public:
        virtual std::string name(const pugi::xml_node &node) = 0;
        virtual std::string value(const pugi::xml_node &node) = 0;

        virtual bool accepts(const pugi::xml_node &node) = 0;
        virtual bool is_reference(const pugi::xml_node &node) {
            return value(node).find('@') != std::string::npos;
        };

        virtual ~Source() = default;
    };

    class LegacySource : public Source {
    public:
        std::string name(const pugi::xml_node &node) override {
            return node.child_value("name");
        }

        std::string value(const pugi::xml_node &node) override {
            return node.child_value("value");
        }

        bool accepts(const pugi::xml_node &node) override {
            return node.child("name") && node.child("value");
        }
    };

    class V2Source : public Source{
    public:
        std::string name(const pugi::xml_node &node) override {
            return node.attribute("name").value();
        }

        std::string value(const pugi::xml_node &node) override {
            return node.attribute("value").value();
        }

        bool accepts(const pugi::xml_node &node) override {
            return node.attribute("name") && node.attribute("value");
        }
    };


    class ValueBuilder : public PropertyBuilder {
    public:
        explicit ValueBuilder(Source &source) : source(source) {}

        bool accepts(const pugi::xml_node &node) const override {
            return source.accepts(node) && !source.is_reference(node);
        }

        std::unique_ptr<Property> build(const pugi::xml_node &node) const override {
            return std::make_unique<ValueProperty>(
                source.name(node),
                source.value(node)
            );
        }

    private:
        class ValueProperty : public Property {
        public:
            ValueProperty(std::string name, std::string value)
            : name_(std::move(name)), value_(std::move(value)) {}

            std::string name() override { return name_; }
            std::string value() override { return value_; }
        private:
            std::string name_, value_;
        };

        Source &source;
    };

    class ReferenceBuilder : public PropertyBuilder {
    public:
        ReferenceBuilder(Source &source, PropertyMap &properties)
        : source(source), properties(properties) {}

        bool accepts(const pugi::xml_node &node) const override {
            return source.accepts(node) && source.is_reference(node);
        }

        std::unique_ptr<Property> build(const pugi::xml_node &node) const override {
            return std::make_unique<ReferenceProperty>(
                source.name(node),
                source.value(node),
                properties
            );
        }

    private:
        class ReferenceProperty : public Property {
        public:
            ReferenceProperty(std::string name, std::string reference, PropertyMap &properties)
            : name_(std::move(name)), reference_(std::move(reference)), properties(properties) {}

            std::string name() override { return name_; }
            std::string value() override {
                std::set<std::string> visited;
                return value(visited);
            }

            std::string value(std::set<std::string> &visited) override {
                if (visited.count(reference_)) {
                    throw std::runtime_error("Cyclical property reference: " + reference_);
                }

                visited.insert(reference_);
                return properties.at(reference_)->value(visited);
            }

        private:
            std::string name_, reference_;
            PropertyMap &properties;
        };

        Source &source;
        PropertyMap &properties;
    };

    class Parser {
    public:
        virtual Config parse(const pugi::xml_document &) = 0;

        virtual ~Parser() = default;

    protected:
        PropertyMap referenceable_properties;
        PropertyBuilders property_builders;

        virtual void assemble_referenceable_properties(const pugi::xml_node &root) {

            for (auto selector : root.select_nodes("//*[child::name and child::property]")) {
                auto node = selector.node();
                auto parent_name = node.child_value("name");

                for (auto p_node : node.children("property")) {
                    auto property = parse_property(p_node);
                    auto location = property->name() + "@" + parent_name;

                    referenceable_properties[location] = std::move(property);
                }
            }
        }

        std::unique_ptr<Property> parse_property(const pugi::xml_node &node) {

            std::vector<std::unique_ptr<Property>> properties;

            for (auto &builder : property_builders) {
                if (builder->accepts(node)) {
                    properties.push_back(builder->build(node));
                }
            }

            if (0 == properties.size()) {
                throw_node_error("Unable to parse property:", node);
            }

            if (2 <= properties.size()) {
                throw_node_error("Ambiguous property parse:", node);
            }

            return std::move(properties[0]);
        }

        std::unordered_map<std::string, std::string>
        parse_properties(const pugi::xml_node &gadget_node) {

            std::unordered_map<std::string, std::string> properties;

            for (auto &node : gadget_node.children("property")) {
                auto property = parse_property(node);
                properties[property->name()] = property->value();
            }

            return properties;
        }

        Config::Reader parse_reader(const pugi::xml_node &reader_node) {

            std::string slot_str = reader_node.child_value("slot");

            boost::optional<uint16_t> slot = boost::none;
            if (!slot_str.empty())
                slot = static_cast<uint16_t>(std::stoi(slot_str));

            return Config::Reader{reader_node.child_value("dll"),
                                  reader_node.child_value("classname"),
                                  slot};
        }

        std::vector<Config::Reader> parse_readers(const pugi::xml_node &reader_root) {
            std::vector<Config::Reader> readers{};
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
    };

    class Legacy : public Parser {
    public:
        Legacy() {
            property_builders.push_back(std::make_unique<ValueBuilder>(legacy_source));
            property_builders.push_back(std::make_unique<ReferenceBuilder>(legacy_source, referenceable_properties));
        }

        Config parse(const pugi::xml_document &config) override {

            assemble_referenceable_properties(config);

            pugi::xml_node root = config.child("gadgetronStreamConfiguration");

            return Config{
                parse_readers(root),
                parse_writers(root),
                parse_stream(root)
            };
        }

    private:
        LegacySource legacy_source;

        std::vector<Config::Gadget> parse_gadgets(const pugi::xml_node& gadget_node){
            std::vector<Config::Gadget> gadgets{};
            for (const auto& node : gadget_node.children("gadget")){
                gadgets.push_back(parse_gadget(node));
            }
            return gadgets;
        }

        Config::Stream parse_stream(const pugi::xml_node& stream_node){
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


    class V2 : public Parser {
    public:
        V2() {
            node_parsers["gadget"] = [&](const pugi::xml_node& n){ return this->parse_gadget(n); };
            node_parsers["parallel"] = [&](const pugi::xml_node& n){ return this->parse_parallel(n); };

            property_builders.push_back(std::make_unique<ValueBuilder>(v2_source));
            property_builders.push_back(std::make_unique<ReferenceBuilder>(v2_source, referenceable_properties));
            property_builders.push_back(std::make_unique<ValueBuilder>(legacy_source));
            property_builders.push_back(std::make_unique<ReferenceBuilder>(legacy_source, referenceable_properties));
        }

        Config parse(const pugi::xml_document& config) override {

            assemble_referenceable_properties(config);

            auto root = config.child("configuration");

            return Config{
                parse_readers(root.child("readers")),
                parse_writers(root.child("writers")),
                parse_stream(root.child("stream"))
            };
        }
    private:
        std::unordered_map<std::string,std::function<Config::Node(const pugi::xml_node&)>> node_parsers;

        LegacySource legacy_source;
        V2Source     v2_source;

        Config::Merge parse_mergenode(const pugi::xml_node& merge_node){
            return Config::Merge{merge_node.child_value("name"), merge_node.child_value("dll"),
                                 merge_node.child_value("classname"), parse_properties(merge_node)};
        }

        Config::Branch parse_branchnode(const pugi::xml_node& branch_node){
            return Config::Branch{branch_node.child_value("name"), branch_node.child_value("dll"),
                                  branch_node.child_value("classname"), parse_properties(branch_node)};
        }

        Config::Parallel parse_parallel(const pugi::xml_node& parallel_node){

            auto branch = parse_branchnode(parallel_node.child("branch"));
            auto merge = parse_mergenode(parallel_node.child("merge"));

            std::vector<Config::Stream> streams{};
            for (const auto& stream_node : parallel_node.children("stream")){
                streams.push_back(parse_stream(stream_node));
            }
            return Config::Parallel{branch, merge, streams};
        }

        Config::Stream parse_stream(const pugi::xml_node& stream_node ){
            std::vector<Config::Node> nodes;
            for (auto& node : stream_node.children() ){
                nodes.push_back(node_parsers.at(node.name())(node));
            }
            return Config::Stream{stream_node.attribute("key").value(), nodes};
        }
    };

   std::unique_ptr<Parser> select_config_parser(const pugi::xml_document &raw_config) {
        if (raw_config.child("gadgetronStreamConfiguration")){
            return std::make_unique<Legacy>();
        } else {
            return std::make_unique<V2>();
        }
    }
}

namespace Gadgetron::Server::Connection {

    Config parse_config(std::istream &stream) {

        pugi::xml_document doc;
        pugi::xml_parse_result result = doc.load(stream);

        if (result.status != pugi::status_ok) {
            GERROR("Loading config file failed with following error: %s (%d)\n", result.description(), result.status);
            throw std::runtime_error(result.description());
        }

        auto parser = select_config_parser(doc);
        return parser->parse(doc);
    }
}
