#pragma once

#include <ostream>

#include "Writer.h"
#include "mri_core_acquisition_bucket.h"


namespace Gadgetron::Core::Writers {

    class AcquisitionBucketWriter : public TypedWriter<IsmrmrdAcquisitionBucket> {
    protected:
        void serialize(std::ostream& stream, const IsmrmrdAcquisitionBucket& args) override;
    };
}
