#include "Python.h"

#include <list>
#include <boost/process.hpp>
#include <boost/optional.hpp>

#include "connection/Config.h"

#include "log.h"

namespace Gadgetron::Server::Connection::Stream {

    void start_python_module(const Config::Execute &execute, unsigned short port, const Gadgetron::Core::Context &context) {

        auto python_path = (context.paths.gadgetron_home / "share" / "gadgetron" / "python").string();

        std::list<std::string> args{
                "-m", "gadgetron",
                std::to_string(port),
                execute.name
        };

        if(execute.target) args.push_back(execute.target.get());

        boost::process::child module(
                boost::process::search_path("python3"),
                boost::process::args=args,
                boost::process::env["PYTHONPATH"]+={python_path}
        );

        module.detach();
        GINFO_STREAM("Started external Python module (pid: " << module.id() << ").");
    }
}
