#include "Python.h"

#include <boost/process.hpp>

#include "log.h"

#define MODE "GADGETRON_EXTERNAL_MODE"
#define PORT "GADGETRON_EXTERNAL_PORT"
#define PATH "PYTHONPATH"

namespace Gadgetron::Server::Connection::Stream {

    void start_python_module(std::string module, unsigned short port, const Gadgetron::Core::Context &context) {

        auto python_path = (context.paths.gadgetron_home / "share" / "gadgetron" / "python").string();

        GINFO_STREAM("Starting Python module '" << module << "'");

        boost::process::spawn(
                boost::process::search_path("python3"),
                "-m", "gadgetron", module,
                boost::process::env[MODE]="active",
                boost::process::env[PORT]=std::to_string(port),
                boost::process::env[PATH]+={python_path}
        );
    }
}
