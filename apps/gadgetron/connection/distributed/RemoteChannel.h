#pragma once

#include <map>
#include <vector>
#include <atomic>
#include <iostream>
#include <boost/asio/ip/tcp.hpp>

#include "Channel.h"
#include "Reader.h"
#include "Writer.h"

namespace Gadgetron::Server::Distributed {

    class RemoteChannel : public Core::Channel {
    public:


        RemoteChannel(
                const Address& address,
                const std::string& xml_config,
                const ISMRMRD::IsmrmrdHeader& header,
                const std::map<uint16_t,std::unique_ptr<Core::Reader>>& readers,
                const std::vector<std::unique_ptr<Core::Writer>>& writers
        );

        Core::Message pop() override;
        Core::optional<Core::Message> try_pop() override {
            return Core::none;
        };

        void push_message(Core::Message ptr) override;
        void close() override;

    private:
        void handle_close();
        void save_error(const std::string& error_message);

        std::unique_ptr<boost::asio::ip::tcp::iostream> stream;
        const Address address;
        const std::map<uint16_t,std::unique_ptr<Core::Reader>>& readers;
        const std::vector<std::unique_ptr<Core::Writer>>& writers;

        std::map<uint16_t, std::function<void(std::istream&)>> info_handlers;

        bool closed_input = false;
        bool closed_output = false;

        std::mutex closed_mutex;

        std::vector<std::string> error_messages;
    };

    class RemoteError : public std::runtime_error {
    public:
        RemoteError(const Address& address, const std::vector<std::string>& errors);
    };

}

