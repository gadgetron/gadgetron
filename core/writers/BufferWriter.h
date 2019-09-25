#pragma once

#include <mri_core_data.h>
#include "Writer.h"

namespace Gadgetron::Core::Writers {

class BufferWriter : public TypedWriter<IsmrmrdReconData> {
public:
    void serialize(std::ostream &stream, const IsmrmrdReconData & args) override;
};
}


