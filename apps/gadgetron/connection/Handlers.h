#pragma once

#include <map>
#include <istream>
#include <functional>

#include "Channel.h"

namespace Gadgetron::Server::Connection::Handlers {


    class Handler {
    public:
        virtual void handle(std::istream &stream, Gadgetron::Core::OutputChannel &channel) = 0;

        virtual ~Handler() = default;
    };

    class QueryHandler : public Handler {
    public:
        QueryHandler();

        explicit QueryHandler(
                std::map<std::string, std::string> additional_answers
        );

        void handle(std::istream &stream, Gadgetron::Core::OutputChannel &channel) override;

        std::map<std::string, std::string> answers;
    };

    class ErrorProducingHandler : public Handler {
    public:
        explicit ErrorProducingHandler(std::string message);

        void handle(std::istream &stream, Gadgetron::Core::OutputChannel &channel) override;

    private:
        std::string message;
    };

    class CloseHandler : public Handler {
    public:
        explicit CloseHandler(std::function<void()> callback);

        void handle(std::istream &stream, Gadgetron::Core::OutputChannel &channel) override;

    private:
        std::function<void()> callback;
    };


}

