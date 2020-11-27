#ifndef DICOMIMAGEWRITER_H
#define DICOMIMAGEWRITER_H

#include <dcmtk/dcmdata/dctk.h>
#include "gadgetron_dicom_export.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd/ismrmrd.h"
#include "Writer.h"

namespace Gadgetron {

    class EXPORTGADGETSDICOM DicomImageWriter
        : public Core::TypedWriter<DcmFileFormat, Core::optional<std::string>, Core::optional<ISMRMRD::MetaContainer>> {
    protected:
        void serialize(std::ostream& stream, const DcmFileFormat&,
            const Core::optional<std::string>&,
            const Core::optional<ISMRMRD::MetaContainer>& args) override;

    protected:
    public:
    };

} /* namespace Gadgetron */

#endif
