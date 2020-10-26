#pragma once

#include "Writer.h"

namespace Gadgetron::Core::Writers
{
    /// \brief This writer does not do anything: it allows to use a pipeline without writing any output.
    class FakeWriter : public Writer
    {
    public:
        bool accepts(const Message& message) override;
        void write(std::ostream &stream, Message message) override;
    };
}

