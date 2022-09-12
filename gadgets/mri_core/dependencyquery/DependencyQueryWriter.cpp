#include "readers/GadgetIsmrmrdReader.h"
#include "DependencyQueryWriter.h"
#include "GadgetContainerMessage.h"

#include "io/primitives.h"

namespace Gadgetron{

    void DependencyQueryWriter::serialize(std::ostream &stream, const DependencyQuery::Dependency& dependencies) {
        using namespace Core;

        IO::write(stream,GADGET_MESSAGE_DEPENDENCY_QUERY);


        auto& meta = dependencies.dependencies;

        std::stringstream meta_stream;
        ISMRMRD::serialize(meta,meta_stream);
        auto meta_string = meta_stream.str();

        IO::write(stream,meta_string.size());
        stream.write(meta_string.data(),meta_string.size());

    }


    GADGETRON_WRITER_EXPORT(DependencyQueryWriter)
}
