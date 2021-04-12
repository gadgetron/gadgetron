#pragma once

namespace Gadgetron::Sessions {

    enum class MessageID : uint64_t {
        INFO = 0,
        STORE = 1,
        QUERY = 2,
        FETCH = 3,
        ERROR = 4,
        OK = 5
    };

    struct FrameHeader {
        MessageID id;
        uint64_t message_length;
    };


    struct Info {

    };


}