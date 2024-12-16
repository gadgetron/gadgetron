#pragma once

#include <thread>
#include <memory>
#include <functional>


#include "io/primitives.h"
#include "Channel.h"
#include "Context.h"

namespace Gadgetron::Main {

    class ErrorReporter {
    public:
        virtual void operator()(const std::string&, const std::string&) = 0;
        virtual ~ErrorReporter() = default;
    };

    class ErrorHandler {
    public:

        ErrorHandler( ErrorReporter &callback, std::string location)
                : location{
                std::move(location)}, push_error{callback} {}

        ErrorHandler( const ErrorHandler &parent,const std::string &new_location) : push_error{parent.push_error}, location{
                add_location(parent.location, new_location)} {

        }

        ErrorHandler(const ErrorHandler &parent) = default;


#if defined(NDEBUG)

        template<class F, class...ARGS>
        void handle(F fn, ARGS &&... args) {
            // When debugging, it is useful to have all exceptions bubble up to the
            // debugger. To enable this, we sabotage the error handler on debug builds.
                try {
                    fn(std::forward<ARGS>(args)...);
                }
                catch (const Core::ChannelClosed &e) {
                    // Ignored.
                }
                catch (const std::exception &e) {
                    push_error(location, e.what());
                }
                catch (...) {
                    push_error(location, "Unknown error.");
                }
        }
#else

        template<class F, class...ARGS>
        void handle(F fn, ARGS &&... args) {
            try {
                fn(std::forward<ARGS>(args)...);
            }
            catch (const Core::ChannelClosed &e) {
                // Ignored.
            }
        }

#endif

        template<class F, class... ARGS>
        std::thread run(F fn, ARGS &&... args) {
            return std::thread(
                    []( auto handler, auto fn, auto &&... iargs) {
                        handler.handle(fn, std::forward<ARGS>(iargs)...);
                    },
                    *this,
                    std::forward<F>(fn),
                    std::forward<ARGS>(args)...
            );
        }

    private:

        static std::string add_location(const std::string &old_location, const std::string &new_location) {
            return old_location + "/" + new_location;
        }

        const std::string location;
        ErrorReporter &push_error;
    };
}