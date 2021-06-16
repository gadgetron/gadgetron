#pragma once

#include <memory>

#include "Serialization.h"
#include "Configuration.h"

#include "Channel.h"

namespace Gadgetron::Server::Connection::Nodes {

    class RemoteError : public std::runtime_error {
    public:
        explicit RemoteError(std::list<std::string> messages);
    private:
        const std::list<std::string> messages;
    };

    struct Local  {};
    struct Remote { std::string address, port; };

    using  Address = Core::variant<Remote, Local>;

    std::ostream & operator<<(std::ostream &, const Local &);
    std::ostream & operator<<(std::ostream &, const Remote &);

    std::string to_string(const Local &);
    std::string to_string(const Remote &);

    std::unique_ptr<std::iostream> listen(std::string port);
    std::unique_ptr<std::iostream> connect(const std::string &address, const std::string &port);
    std::unique_ptr<std::iostream> connect(const Remote &remote);
    std::unique_ptr<std::iostream> connect(const Address &address, std::shared_ptr<Configuration> configuration);
}