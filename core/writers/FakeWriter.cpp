#include "FakeWriter.h"

namespace Gadgetron::Core::Writers
{
bool FakeWriter::accepts(const Message& /*message*/)
{
    return false;
}

void FakeWriter::write(std::ostream& /*stream*/, Message /*message*/)
{

}

GADGETRON_WRITER_EXPORT(FakeWriter)
}
