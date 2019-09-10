#pragma once
#include "Writer.h"
#include "mri_core_data.h"


namespace Gadgetron::Core::Writers {
    class IsmrmrdImageArrayWriter : public TypedWriter<IsmrmrdImageArray> {
    protected:
        void serialize(std::ostream& stream, const IsmrmrdImageArray& args) override;
    };
}
