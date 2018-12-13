#pragma once

#include "Response.h"
#include "Writer.h"

namespace Gadgetron::Core::Writers {
    class ResponseWriter : public TypedWriter<Response> {
    public:
        void serialize(std::ostream &stream, std::unique_ptr<Response> ) override;
    };
}


