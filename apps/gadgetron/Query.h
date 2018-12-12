#pragma once

#include <map>
#include <future>
#include <sstream>

namespace Gadgetron::Server::Query {

    class Handler {
    public:
        ~Handler() = default;

        virtual bool accepts(const std::string &query) = 0;
        virtual std::string handle(const std::string &query) = 0;
    };

    class GadgetronHandler : public Handler {
    public:
        GadgetronHandler();
        bool accepts(const std::string &query) override;
        std::string handle(const std::string &query) override;
    private:
        std::map<std::string, std::function<std::string()>> handlers;
    };

    class ISMRMRDHandler : public Handler {
        bool accepts(const std::string &query) override;
        std::string handle(const std::string &query) override;
    };

    class ConfigHandler : public Handler {
    public:
        explicit ConfigHandler(std::future<std::stringstream> &&raw_config_future);

        bool accepts(const std::string &query) override;
        std::string handle(const std::string &query) override;
    private:
        std::future<std::stringstream> config;
    };
}

