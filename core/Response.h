#pragma once

#include <string>

namespace Gadgetron::Core {
    struct Response {
        Response(uint64_t correlation_id, const std::string &string);

        uint64_t correlation_id;
        std::string response;
    };
}
