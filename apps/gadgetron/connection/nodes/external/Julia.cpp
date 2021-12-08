#include "Python.h"

#include <list>
#include "Process.h"

#include "connection/config/Config.h"

#include <boost/asio.hpp>

#include "log.h"
#include <regex>

namespace Gadgetron::Server::Connection::Nodes {

    using namespace Gadgetron::Core;


    boost::process::child start_julia_module(
        const Config::Execute &execute,
        unsigned short port,
        const StreamContext &context
    ) {
        if (!execute.target) throw std::invalid_argument("Target must be specified for Julia modules");

        std::list<std::string> args{
            "-e", "import Gadgetron; Gadgetron.External.main()",
            std::to_string(port),
            execute.name, *execute.target
        };
        namespace bp = boost::process;
        //Workaround for bug in Boost process
        auto env = boost::this_process::environment();

        env.set("GADGETRON_STORAGE_ADDRESS",context.storage_address);

        auto module = Process::child(
                boost::process::search_path("julia"),
                boost::process::args = args,
                env,
                boost::process::limit_handles,
                boost::process::std_out > stdout,
                boost::process::std_err > stderr
        );

        GINFO_STREAM("Started external Julia module (pid: " << module.id() << ").");
        return std::move(module);
    }

    bool julia_available() noexcept {
        try {
            return !Process::system(
                    boost::process::search_path("julia"),
                    boost::process::args={"-e", "using Gadgetron"},
                    boost::process::std_out > boost::process::null,
                    boost::process::std_err > boost::process::null
            );
        }
        catch (...) {
            return false;
        }
    }
}
