#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "Types.h"

namespace Gadgetron::Server::Connection {

    struct Config {

        struct Gadget;
        struct External;
        struct Parallel;
        struct Distributed;
        struct ParallelProcess;
        struct PureDistributed;
        using Node = Core::variant<Gadget, External, Parallel, Distributed, ParallelProcess, PureDistributed>;

        template<class CONFIG>
        static std::string name(CONFIG config) {
            return config.name.empty() ? config.classname : config.name;
        }

        struct Reader {
            std::string dll, classname;
            Core::optional<uint16_t> slot;
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
            Gadget(std::string name, std::string dll, std::string classname, std::unordered_map<std::string, std::string> properties):
            name(std::move(name)), dll(std::move(dll)), classname(std::move(classname)), properties(std::move(properties))
            {

            }
        };

        struct Execute {
            std::string name, type;
            Core::optional<std::string> target;
        };

        struct Connect {
            std::string address, port;
        };

        using Action = Core::variant<Execute, Connect>;

        struct External {
            Action action;

            struct Configuration;
            std::shared_ptr<Configuration> configuration;

            std::vector<Reader> readers;
            std::vector<Writer> writers;
        };

        struct Branch : Gadget { using Gadget::Gadget;};
        struct Merge : Gadget { using Gadget::Gadget;};

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

        struct Distributor : Gadget { using Gadget::Gadget;};

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
    std::string serialize_config(const Config::External& external_config);
}
