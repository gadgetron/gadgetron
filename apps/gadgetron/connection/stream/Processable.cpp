#include "Processable.h"


std::thread Gadgetron::Server::Connection::Stream::Processable::process_async(std::shared_ptr<Processable> processable, Core::InputChannel input,
                                           Core::OutputChannel output, const ErrorHandler &errorHandler) {
        return std::thread([](auto processable, auto input, auto output, auto error_handler) {
            error_handler.handle([processable, &error_handler](auto input, auto output) {
                                     processable->process(std::move(input), std::move(output), error_handler);
                                 },
                                 std::move(input), std::move(output));
        }, processable, std::move(input), std::move(output), ErrorHandler{errorHandler,processable->name()});

    }
