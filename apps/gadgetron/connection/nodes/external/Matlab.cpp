#include "Matlab.h"

#include <list>
#include <boost/process.hpp>

#include "connection/config/Config.h"

#include "log.h"

namespace Gadgetron::Server::Connection::Nodes {

    boost::process::child start_matlab_module(
        const Config::Execute &execute,
        unsigned short port,
        const Gadgetron::Core::StreamContext &context
    ) {
        boost::process::child module(
                boost::process::search_path("matlab"),
                boost::process::args={"-batch", "gadgetron.external.main"},
                boost::process::env["GADGETRON_EXTERNAL_PORT"] = std::to_string(port),
                boost::process::env["GADGETRON_EXTERNAL_MODULE"] = execute.name,
                boost::process::env["GADGETRON_STORAGE_ADDRESS"] = context.storage_address
        );

        GINFO_STREAM("Started external MATLAB module (pid: " << module.id() << ").");

        return std::move(module);
    }

    bool matlab_available() noexcept {
        try {
            return !boost::process::system(
                    boost::process::search_path("matlab"),
                    boost::process::args={"-batch", "gadgetron.external.test_available"},
                    boost::process::std_out > boost::process::null,
                    boost::process::std_err > boost::process::null
            );
        }
        catch (...) {
            return false;
        }
    }
}
