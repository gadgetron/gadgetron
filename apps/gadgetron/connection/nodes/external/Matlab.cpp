#include "Matlab.h"

#include <list>
#include "Process.h"
#include "connection/config/Config.h"

#include "log.h"

namespace Gadgetron::Server::Connection::Nodes {

    boost::process::child start_matlab_module(
        const Config::Execute &execute,
        unsigned short port,
        const Gadgetron::Core::StreamContext &context
    ) {

        auto env = boost::this_process::environment();
        env.set("GADGETRON_EXTERNAL_PORT",std::to_string(port));
        env.set("GADGETRON_EXTERNAL_MODULE",execute.name);
        env.set("GADGETRON_STORAGE_ADDRESS", context.storage_address);

        auto module = Process::child(
                boost::process::search_path("matlab"),
                boost::process::args={"-batch", "gadgetron.external.main"},
                env,
                boost::process::limit_handles,
                boost::process::std_out > stdout,
                boost::process::std_err > stderr
        );

        GINFO_STREAM("Started external MATLAB module (pid: " << module.id() << ").");

        return std::move(module);
    }

    bool matlab_available() noexcept {
        try {
            return !Process::system(
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
