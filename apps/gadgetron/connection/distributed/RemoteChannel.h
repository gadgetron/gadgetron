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

        RemoteChannel(const Address& address, const std::string& xml_config, const std::map<uint16_t,std::shared_ptr<Core::Reader>>& readers,const std::vector<std::shared_ptr<Core::Writer>>& writers );


        std::unique_ptr<Core::Message> pop() override;

        void push_message(std::unique_ptr<Core::Message> &&ptr) override;

        void close() override;

    private:

        std::unique_ptr<std::iostream> stream;
        const std::map<uint16_t,std::shared_ptr<Core::Reader>> readers;
        const std::vector<std::shared_ptr<Core::Writer>>& writers;

        std::map<uint16_t, std::function<void(std::ostream&)>> error_readers;

        bool closed = false;
        std::mutex closed_mutex;


    };

}

