#pragma once

#include "Writer.h"
#include "Response.h"

namespace Gadgetron::Server::Connection::Writers {

    class ResponseWriter : public Core::TypedWriter<Core::Response> {
    public:
        void serialize(std::ostream &, std::unique_ptr<Core::Response>) override;
    };

    class TextWriter : public Core::TypedWriter<std::string> {
    public:
        void serialize(std::ostream &, std::unique_ptr<std::string>) override;
    };

}



