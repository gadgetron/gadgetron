#include "External.h"

namespace {

}

namespace Gadgetron::Server::Connection::Stream {

    External::External(
            const Config::External &config,
            const Core::Context &context,
            Loader &loader
    ) {



    }

    void External::process(
            Core::InputChannel input,
            Core::OutputChannel output,
            ErrorHandler &error_handler
    ) {

    }

    const std::string &External::name() {
        static const std::string name = "external";
        return name;
    }
}