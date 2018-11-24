#pragma once

#include <boost/variant.hpp>
#include <boost/optional.hpp>
namespace Gadgetron::AST {

    struct Gadget;
    struct Stream;
    struct Parallel;
    using Node = boost::variant<Gadget, Parallel>;

    struct Gadget {
        std::string name;
        std::string dll;
        std::string classname;

        std::unordered_map<std::string, std::string> properties;

    };

    struct Reader {
        std::string dll;
        std::string classname;
        boost::optional<uint16_t> port;
    };

    struct Writer {
        std::string dll;
        std::string classname;
        boost::optional<uint16_t> port;
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



    struct Chain { //Not sure I like this name
        std::vector<Reader> readers;
        std::vector<Writer> writers;
        Stream stream;

    };
}