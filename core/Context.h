#pragma once
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <mrd/types.h>

namespace Gadgetron::Core {

    struct Context {
        using Header = mrd::Header;

        struct Paths {
            boost::filesystem::path gadgetron_home;
        };

        Header header;
        Paths  paths;
        std::map<std::string, std::string> parameters;
    };

    struct StreamContext : Context {
        using Args = boost::program_options::variables_map;

        StreamContext(
            mrd::Header header,
            const Paths paths,
            const Args args
        ) : Context{
                std::move(header),
                paths,
                GetParameters(args)
            },
            args{args} {}


        Args args;

        private:
        static std::map<std::string, std::string> GetParameters(const boost::program_options::variables_map& args) {
            std::map<std::string, std::string> parameters;
            if (args.count("parameter")) {
                auto params = args["parameter"].as<std::vector<std::pair<std::string, std::string>>>();
                for (auto &arg : params) {
                    parameters[arg.first] = arg.second;
                }
            }
            return parameters;
        }

    };
}
