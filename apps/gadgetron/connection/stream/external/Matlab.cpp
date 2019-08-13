#include "Matlab.h"

#include <list>
#include <boost/process.hpp>
#include <boost/optional.hpp>

#include "connection/Config.h"

#include "log.h"

namespace Gadgetron::Server::Connection::Stream {

    boost::process::child start_matlab_module(const Config::Execute &execute, unsigned short port, const Gadgetron::Core::Context &context) {

        auto matlab_path = (context.paths.gadgetron_home / "share" / "gadgetron" / "matlab").string();

        boost::process::child module(
                boost::process::search_path("matlab"),
                boost::process::args={"-batch", "gadgetron.external.main"},
                boost::process::env["MATLABPATH"]+={matlab_path},
                boost::process::env["GADGETRON_EXTERNAL_PORT"] = std::to_string(port),
                boost::process::env["GADGETRON_EXTERNAL_MODULE"] = execute.name
        );

        GINFO_STREAM("Started external MATLAB module (pid: " << module.id() << ").");

        return std::move(module);
    }
}
