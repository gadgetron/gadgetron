#include "Processable.h"


std::thread Gadgetron::Main::Processable::process_async(
    std::shared_ptr<Processable> processable,
    Core::GenericInputChannel input,
    Core::OutputChannel output,
    const ErrorHandler &error_handler
) {
    ErrorHandler nested_handler{error_handler, processable->name()};

    return nested_handler.run(
        [=](auto input, auto output, auto error_handler) {
          processable->process(std::move(input), std::move(output), error_handler);
        },
        std::move(input),
        std::move(output),
        nested_handler
    );
}
