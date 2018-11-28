#pragma once

#include <string>
#include <unordered_map>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

namespace Gadgetron::Server {

    struct Config {

        struct Gadget;
        struct Stream;
        struct Parallel;
        using Node = boost::variant<Gadget, Parallel>;

        struct Gadget {
            std::string name, dll, classname;
            std::unordered_map<std::string, std::string> properties;
        };

        struct Reader {
            std::string dll, classname;
            boost::optional<uint16_t> port;
        };

        struct Writer {
            std::string dll, classname;
        };

        struct Stream {
            std::string name;
            std::vector<Node> nodes;
        };

        struct BranchNode : Gadget {

        };

        struct MergeNode : Gadget {

        };

        struct Parallel {
            BranchNode branch;
            MergeNode merge;
            std::vector<Stream> streams;
        };

        std::vector<Reader> readers;
        std::vector<Writer> writers;
        Stream stream;
    };

    Config parse_config(std::istream &stream);
}