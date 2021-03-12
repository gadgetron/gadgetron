#include "Python.h"

#include <list>
#include <boost/process.hpp>
#include <boost/optional.hpp>

#include "connection/Config.h"

#include <boost/asio.hpp>

#include "log.h"
#include <regex> 
namespace Gadgetron::Server::Connection::Stream {

    namespace {

        std::tuple<int, int, int> parse_version(const std::string& versionstring) {

            auto version_rules = std::regex(R"(.*?([0-9]+)\.([0-9]+)\.([0-9]+))");
            std::smatch matches;
            auto found = std::regex_search(versionstring, matches, version_rules);

            if (!found)
                return {0, 0, 0};

            int major_version = std::stoi(matches[1].str());
            int minor_version = std::stoi(matches[2].str());
            int patch_version = std::stoi(matches[3].str());

            return {major_version, minor_version, patch_version};
        }

        bool is_valid_python3(const std::string& pythonname){
        try {
            boost::asio::io_service ios;
            std::future<std::string> output_stream;

            auto c = boost::process::child(
                    boost::process::search_path(pythonname),
                    boost::process::args={"--version"},
                    boost::process::std_in.close(),
                    boost::process::std_out > output_stream,
                    boost::process::std_err > boost::process::null,
                    ios
            );
            ios.run();

            auto output = output_stream.get();

            
            auto [major, minor, patch] = parse_version(output);

            if ((major == 3) && (minor > 5))
                return true;


       }
        catch (...){
            return false;
        }
       return false;


    };

    const std::string& get_python_executable() {
        static std::string python_executable_name = "";
        if (python_executable_name.empty()) {
            auto possible_names = std::vector<std::string>{"python", "python3"};
            
            for (auto& name : possible_names) {
                if (is_valid_python3(name)) {
                    python_executable_name = name;
                    return python_executable_name;
                } 
            }
            python_executable_name = " ";
        }

        return python_executable_name;
    }
    }

    boost::process::child start_python_module(const Config::Execute &execute, unsigned short port, const Gadgetron::Core::Context &context) {

        auto python_path = (context.paths.gadgetron_home / "share" / "gadgetron" / "python").string();

        std::list<std::string> args{
                "-m", "gadgetron",
                std::to_string(port),
                execute.name
        };

        if(execute.target) args.push_back(execute.target.value());

        boost::process::child module(
                boost::process::search_path(get_python_executable()),
                boost::process::args=args,
                boost::process::env["PYTHONPATH"]+={python_path}
        );

        GINFO_STREAM("Started external Python module (pid: " << module.id() << ").");
        return std::move(module);
    }

    bool python_available() noexcept {
        try {
            return !boost::process::system(
                    boost::process::search_path(get_python_executable()),
                    boost::process::args={"-m", "gadgetron"},
                    boost::process::std_out > boost::process::null,
                    boost::process::std_err > boost::process::null
            );
        }
        catch (...) {
            return false;
        }
    }
}
