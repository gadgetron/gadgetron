#pragma once

namespace Gadgetron::Core::AST{

    struct Gadget {
        std::string name;
        std::string dll;
        std::string classname;

        std::unordered_map<std::string,std::string> properties;

    };

    struct Reader {
        std::string dll;
        std::string classname;
        boost::optional<int16_t> port;
    };

    struct Writer {
        std::string dll;
        std::string classname;
        boost::optional<int16_t> port;
    };

    struct Stream {
        std::string name;
        std::vector<Gadget> gadgets;
    };

    struct BranchNode : Gadget {};
    struct MergeNode : Gadget {};

    struct Parallel {
        BranchNode branch;
        MergeNode merge;
        std::vector<Stream> streams;
    };
}