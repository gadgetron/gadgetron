#include "TextWriter.h"
#include "io/primitives.h"
#include "MessageID.h"

namespace Gadgetron::Core::Writers
{
    void TextWriter::serialize(std::ostream &stream, const std::string& str)
    {
        IO::write(stream, TEXT);
        IO::write_string_to_stream<uint32_t>(stream, str);
    }

    GADGETRON_WRITER_EXPORT(TextWriter);
}


