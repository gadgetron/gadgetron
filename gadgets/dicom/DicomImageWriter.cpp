#include <complex>
#include <fstream>
#include <time.h>

// Gadgetron includes
#include "GadgetIsmrmrdReadWrite.h"
#include "DicomImageWriter.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"

// DCMTK includes
#include "dcmtk/config/osconfig.h"
#include "dcmtk/ofstd/ofstdinc.h"
#define INCLUDE_CSTDLIB
#define INCLUDE_CSTDIO
#define INCLUDE_CSTRING
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcostrmb.h"


namespace Gadgetron {

int DicomImageWriter::write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb)
{
    GadgetContainerMessage<DcmFileFormat>* dcm_file_message = AsContainerMessage<DcmFileFormat>(mb);
    if (!dcm_file_message) {
        ACE_DEBUG( (LM_ERROR, ACE_TEXT("(%P,%l), DicomImageWriter::write, invalid image message objects, 1\n")) );
        return -1;
    }

    DcmFileFormat *dcmFile = dcm_file_message->getObjectPtr();

/* BEGIN DEBUG

    OFString modality;
    DcmTagKey key(0x0008, 0x0060);
    OFCondition s = dcmFile->getDataset()->findAndGetOFString(key, modality);
    if (s.bad()) {
        GADGET_DEBUG1("Failed to set Modality\n");
        return GADGET_FAIL;
    }

    GADGET_DEBUG2("Verifying that DcmDataset is valid... Modality: %s\n", modality.c_str());

END DEBUG */

    //GADGET_DEBUG1("Initializing transfer state for DICOM file\n");
    // Initialize transfer state of DcmDataset
    dcmFile->transferInit();

    // Calculate size of DcmFileFormat and create a SUFFICIENTLY sized buffer
    long buffer_length = dcmFile->calcElementLength(EXS_LittleEndianExplicit, EET_ExplicitLength) * 2;
    char buffer[buffer_length];

    DcmOutputBufferStream out_stream(buffer, buffer_length);

    OFCondition status;

    status = dcmFile->write(out_stream, EXS_LittleEndianExplicit, EET_ExplicitLength, NULL);
    if (!status.good()) {
        GADGET_DEBUG2("Failed to write DcmFileFormat to DcmOutputStream(%s)\n", status.text());
        return GADGET_FAIL;
    }

    void *serialized = NULL;
    offile_off_t serialized_length = 0;
    out_stream.flushBuffer(serialized, serialized_length);

    // finalize transfer state of DcmDataset
    dcmFile->transferEnd();

    ssize_t send_cnt = 0;

    GadgetMessageIdentifier id;
    id.id = GADGET_MESSAGE_DICOM;
    //GADGET_DEBUG2("Sending GadgetMessageIdentifier %d\n", id.id);
    if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0) {
        ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send DICOM message identifier\n")));
        return -1;
    }
    //GADGET_DEBUG2("Sent GadgetMessageIdentifier %d\n", id.id);


    uint32_t nbytes = (uint32_t)serialized_length;
    //GADGET_DEBUG2("Sending bytes length %d\n", serialized_length);
    if ((send_cnt = sock->send_n (&nbytes, sizeof(nbytes))) <= 0) {
        ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send DICOM bytes length\n")));
        return -1;
    }
    //GADGET_DEBUG2("Sent bytes length %d\n", serialized_length);


    //GADGET_DEBUG1("Begin sending DICOM image bytes\n");
    if ((send_cnt = sock->send_n (serialized, serialized_length)) <= 0) {
        ACE_DEBUG ((LM_ERROR, ACE_TEXT ("(%P|%t) Unable to send DICOM bytes\n")));
        return -1;
    }
    //GADGET_DEBUG1("Finished sending DICOM image bytes\n");

    return 0;
}

GADGETRON_WRITER_FACTORY_DECLARE(DicomImageWriter)

} /* namespace Gadgetron */
