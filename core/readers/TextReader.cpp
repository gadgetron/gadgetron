
#include "io/primitives.h"
#include "MessageID.h"
#include "TextReader.h"

namespace Gadgetron::Core::Readers
{
    Core::Message TextReader::read(std::istream &stream)
    {
        std::string str = IO::read_string_from_stream<uint32_t>(stream);
        return Message(str);
    }

    uint16_t TextReader::slot() {
        return TEXT;
    }

    GADGETRON_READER_EXPORT(TextReader)
}