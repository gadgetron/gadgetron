#pragma once

#include <map>
#include <istream>
#include <functional>

#include "Channel.h"

namespace Gadgetron::Server::Connection::Handlers {

    enum message_id : uint16_t {
        FILENAME    = 1,
        CONFIG      = 2,
        HEADER      = 3,
        CLOSE       = 4,
        TEXT        = 5,
        QUERY       = 6,
        RESPONSE    = 7,
        ERROR       = 8
    };

    class Handler {
    public:
        virtual void handle(std::istream &stream) = 0;
        virtual ~Handler() = default;
    };

    class QueryHandler : public Handler {
    public:
        explicit QueryHandler(Gadgetron::Core::OutputChannel &channel);
        QueryHandler(
                Gadgetron::Core::OutputChannel &channel,
                std::map<std::string, std::string> additional_answers
        );

        void handle(std::istream &stream) override;

        std::map<std::string, std::string> answers;
        Gadgetron::Core::OutputChannel &channel;
    };

    class ErrorProducingHandler : public Handler {
    public:
        explicit ErrorProducingHandler(std::string message);

        void handle(std::istream &stream) override;

    private:
        std::string message;
    };

    class CloseHandler : public Handler {
    public:
        explicit CloseHandler(std::function<void()> callback);
        void handle(std::istream &stream) override;

    private:
        std::function<void()> callback;
    };


}

