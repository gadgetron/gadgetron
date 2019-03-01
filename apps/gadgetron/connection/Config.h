#pragma once

#include <string>
#include <unordered_map>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

namespace Gadgetron::Server::Connection {

    struct Config {

        struct Gadget;
        struct External;
        struct Parallel;
        struct Distributed;
        struct ParallelProcess;
        struct PureDistributed;
        using Node = boost::variant<Gadget, External, Parallel, Distributed, ParallelProcess, PureDistributed>;

        template<class CONFIG>
        static std::string name(CONFIG config) {
            return config.name.empty() ? config.classname : config.name;
        }

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

        struct PureStream{
            std::vector<Gadget> gadgets;
        };

        struct Gadget {
            std::string name, dll, classname;
            std::unordered_map<std::string, std::string> properties;
        };

        struct Module { std::string name, type; };
        struct Connect { std::string port; };
        using Action = boost::variant<Module, Connect>;

        struct External {
            Action action;

            struct Configuration;
            std::shared_ptr<Configuration> configuration;

            std::vector<Reader> readers;
            std::vector<Writer> writers;
        };

        struct Branch : Gadget {};
        struct Merge : Gadget {};

        struct Parallel {
            Branch branch;
            Merge merge;
            std::vector<Stream> streams;
        };

        struct PureDistributed {
            std::vector<Reader> readers;
            std::vector<Writer> writers;
            PureStream stream;
        };

        struct ParallelProcess {

            size_t workers = 0;
            PureStream stream;
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
    std::string serialize_external_config(const Config::External& external_config);
}