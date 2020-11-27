#include <complex>
#include <fstream>
#include <io/primitives.h>
#include <time.h>

// Gadgetron includes
#include "DicomImageWriter.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"

// DCMTK includes
#define INCLUDE_CSTDLIB
#define INCLUDE_CSTDIO
#define INCLUDE_CSTRING
#define NOGDI
#include "dcmtk/config/osconfig.h"
#include "dcmtk/ofstd/ofstdinc.h"

#include "dcmtk/dcmdata/dcostrmb.h"
#include "dcmtk/dcmdata/dctk.h"


namespace Gadgetron {
    void DicomImageWriter::serialize(std::ostream& stream, const DcmFileFormat& dcmInput,
        const Core::optional<std::string>& dcm_filename_message,
        const Core::optional<ISMRMRD::MetaContainer>& dcm_meta_message) {
        using namespace Gadgetron::Core;


        auto dcmFile = DcmFileFormat(dcmInput); //Copy here, because DCMTK has decided that the const keyword is bad
        // Initialize transfer state of DcmDataset
        dcmFile.transferInit();

        // Calculate size of DcmFileFormat and create a SUFFICIENTLY sized buffer
        long buffer_length = dcmFile.calcElementLength(EXS_LittleEndianExplicit, EET_ExplicitLength) * 2;
        std::vector<char> bufferChar(buffer_length);
        char* buffer = &bufferChar[0];

        DcmOutputBufferStream out_stream(buffer, buffer_length);

        OFCondition status;

        status = dcmFile.write(out_stream, EXS_LittleEndianExplicit, EET_ExplicitLength, NULL);

        offile_off_t serialized_length = 0;

        char* serialized;
        out_stream.flushBuffer(reinterpret_cast<void*&>(serialized), serialized_length);

        // finalize transfer state of DcmDataset
        dcmFile.transferEnd();


        Core::IO::write(stream, GADGET_MESSAGE_DICOM_WITHNAME);

        uint32_t nbytes = (uint32_t)serialized_length;
        Core::IO::write(stream, nbytes);

        stream.write(serialized, serialized_length);

        // chech whether the image filename is attached
        if (dcm_filename_message) {
            unsigned long long len = dcm_filename_message->length();
            Core::IO::write_string_to_stream<unsigned long long>(stream, *dcm_filename_message);
        } else {
            Core::IO::write(stream, (unsigned long long)0);
        }

        if (dcm_meta_message) {
                std::stringstream str;
                ISMRMRD::serialize(*dcm_meta_message, str);
                std::string attribContent = str.str();
                Core::IO::write_string_to_stream(stream, attribContent);
        } else {
            Core::IO::write(stream, (unsigned long long)0);
        }

    }
    GADGETRON_WRITER_EXPORT(DicomImageWriter)
} /* namespace Gadgetron */
