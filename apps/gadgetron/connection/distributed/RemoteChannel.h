#pragma once

#include "Channel.h"
#include "Reader.h"
#include "Writer.h"
#include <map>
#include <vector>
#include <atomic>
#include <iostream>

namespace Gadgetron::Server::Distributed{

struct Address {
    std::string ip;
    std::string port;
};

class RemoteChannel : public Core::Channel {
    public:

        RemoteChannel(const Address& address, const std::string& xml_config, const ISMRMRD::IsmrmrdHeader& header, const std::map<uint16_t,std::unique_ptr<Core::Reader>>& readers,const std::vector<std::unique_ptr<Core::Writer>>& writers );


       Core::Message pop() override;

        void push_message(Core::Message ptr) override;

        void close() override;

    private:

        void handle_close();
        void save_error(const std::string& error_message);

        std::unique_ptr<std::iostream> stream;
        const Address address;
        const std::map<uint16_t,std::unique_ptr<Core::Reader>>& readers;
        const std::vector<std::unique_ptr<Core::Writer>>& writers;

        std::map<uint16_t, std::function<void(std::istream&)>> info_handlers;

        bool closed = false;
        std::mutex closed_mutex;

        std::vector<std::string> error_messages;


    };

class RemoteError : public std::runtime_error {
public:
    RemoteError(const Address& address, const std::vector<std::string>& errors);
};

}

