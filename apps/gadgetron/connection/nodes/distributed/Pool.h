#pragma once

#include <mutex>
#include <vector>
#include <memory>
#include <future>
#include <algorithm>

#include "connection/nodes/common/External.h"

#include "Worker.h"

#include "Message.h"


namespace Gadgetron::Server::Connection::Nodes {

    class Pool {
    public:
        explicit Pool(std::list<std::unique_ptr<Worker>> workers);
        std::future<Core::Message> push(Core::Message message);

    private:
        std::shared_ptr<std::list<std::unique_ptr<Worker>>> workers;
    };
}