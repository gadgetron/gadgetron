#pragma once

#include "Writer.h"

namespace Gadgetron::Core::Writers {

    class TextWriter : public TypedWriter<std::string>
    {
    protected:
        void serialize(std::ostream &stream, const std::string& str) override;
    };
}