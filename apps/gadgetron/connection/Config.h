#pragma once

#include <string>
#include <unordered_map>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

namespace Gadgetron::Server::Connection {

    struct Config {

        struct Gadget;
        struct Parallel;
        struct Distributed;
        using Node = boost::variant<Gadget, Parallel, Distributed>;

        struct Reader {
            std::string dll, classname;
            boost::optional<uint16_t> slot;
        };

        struct Writer {
            std::string dll, classname;
        };

        struct Stream {
            std::string key;
            std::vector<Node> nodes;
        };

        struct Gadget {
            std::string name, dll, classname;
            std::unordered_map<std::string, std::string> properties;
        };

        struct Branch : Gadget {};
        struct Merge : Gadget {};

        struct Parallel {
            Branch branch;
            Merge merge;
            std::vector<Stream> streams;
        };

        struct Distributor : Gadget {};


        struct Distributed {
            std::vector<Reader> readers;
            std::vector<Writer> writers;
            Distributor distributor;
            Stream stream;
        };

        std::vector<Reader> readers;
        std::vector<Writer> writers;
        Stream stream;
    };

    Config parse_config(std::istream &stream);
    std::string serialize_config(const Config& config);
}