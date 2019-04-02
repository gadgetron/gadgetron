#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

namespace Gadgetron::Server::Connection {

    struct Config {

        struct Gadget;
        struct Parallel;
        struct Distributed;
        struct ParallelProcess;
        struct PureDistributed;
        using Node = boost::variant<Gadget, Parallel, Distributed, ParallelProcess, PureDistributed>;

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
            Gadget(std::string name, std::string dll, std::string classname, std::unordered_map<std::string, std::string> properties):
            name(std::move(name)), dll(std::move(dll)), classname(std::move(classname)), properties(std::move(properties))
            {

            }
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
}