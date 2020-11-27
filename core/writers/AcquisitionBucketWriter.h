#pragma once

#include <ostream>

#include "Writer.h"
#include "mri_core_acquisition_bucket.h"


namespace Gadgetron::Core::Writers {

    class AcquisitionBucketWriter : public TypedWriter<AcquisitionBucket> {
    protected:
        void serialize(std::ostream& stream, const AcquisitionBucket& args) override;
    };
}
