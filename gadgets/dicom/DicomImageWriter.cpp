#include <complex>
#include <fstream>
#include <time.h>

// Gadgetron includes
#include "GadgetIsmrmrdReadWrite.h"
#include "DicomImageWriter.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"

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
    if (!dcm_file_message)
    {
      GERROR("DicomImageWriter::write, invalid image message objects\n");
      return -1;
    }

    DcmFileFormat *dcmFile = dcm_file_message->getObjectPtr();

    // Initialize transfer state of DcmDataset
    dcmFile->transferInit();

    // Calculate size of DcmFileFormat and create a SUFFICIENTLY sized buffer
    long buffer_length = dcmFile->calcElementLength(EXS_LittleEndianExplicit, EET_ExplicitLength) * 2;
    std::vector<char> bufferChar(buffer_length);
    char* buffer = &bufferChar[0];

    DcmOutputBufferStream out_stream(buffer, buffer_length);

    OFCondition status;

    status = dcmFile->write(out_stream, EXS_LittleEndianExplicit, EET_ExplicitLength, NULL);
    if (!status.good()) {
      GERROR("Failed to write DcmFileFormat to DcmOutputStream(%s)\n", status.text());
      return GADGET_FAIL;
    }

    void *serialized = NULL;
    offile_off_t serialized_length = 0;
    out_stream.flushBuffer(serialized, serialized_length);

    // finalize transfer state of DcmDataset
    dcmFile->transferEnd();

    ssize_t send_cnt = 0;

    GadgetMessageIdentifier id;
    id.id = GADGET_MESSAGE_DICOM_WITHNAME;

    if ((send_cnt = sock->send_n (&id, sizeof(GadgetMessageIdentifier))) <= 0)
    {
      GERROR("Unable to send DICOM message identifier\n");
      return -1;
    }

    uint32_t nbytes = (uint32_t)serialized_length;
    if ((send_cnt = sock->send_n (&nbytes, sizeof(nbytes))) <= 0)
    {
      GERROR("Unable to send DICOM bytes length\n");
      return -1;
    }

    if ((send_cnt = sock->send_n (serialized, serialized_length)) <= 0)
    {
      GERROR("Unable to send DICOM bytes\n");
      return -1;
    }

    // chech whether the image filename is attached
    GadgetContainerMessage<std::string>* dcm_filename_message = AsContainerMessage<std::string>(mb->cont());
    if (dcm_filename_message)
    {
        unsigned long long len = dcm_filename_message->getObjectPtr()->length();
        if ((send_cnt = sock->send_n (&len, sizeof(unsigned long long))) <= 0)
        {
          GERROR("Unable to send DICOM filename length\n");
          return -1;
        }

        const char* filename = dcm_filename_message->getObjectPtr()->c_str();
        if ((send_cnt = sock->send_n (filename, len)) <= 0)
        {
          GERROR("Unable to send DICOM filename\n");
          return -1;
        }

        GadgetContainerMessage<ISMRMRD::MetaContainer>* dcm_meta_message = AsContainerMessage<ISMRMRD::MetaContainer>(dcm_filename_message->cont());

        char* buf = NULL;
        len = 0;

        if (dcm_meta_message)
        {
            try
            {
                std::stringstream str;
                ISMRMRD::serialize( *dcm_meta_message->getObjectPtr(), str);
                std::string attribContent = str.str();
                len = attribContent.length();

                buf = new char[len];
                GADGET_CHECK_THROW(buf != NULL);

                memcpy(buf, attribContent.c_str(), len);
            }
            catch(...)
            {
              GERROR("Unable to serialize dicom image meta attributes \n");
              return -1;
            }

            if ((send_cnt = sock->send_n(&len, sizeof(unsigned long long))) <= 0)
            {
                GERROR("Unable to send dicom image meta attributes length\n");
                if (buf != NULL) delete[] buf;
                return -1;
            }

            if ( (send_cnt = sock->send_n (buf, len)) <= 0 )
            {
              GERROR("Unable to send dicom image meta attributes\n");
              if ( buf != NULL ) delete [] buf;
              return -1;
            }

            if ( buf != NULL ) delete [] buf;
        }
        else
        {
            if ((send_cnt = sock->send_n(&len, sizeof(unsigned long long))) <= 0)
            {
                GERROR("Unable to send dicom image meta attributes length\n");
                return -1;
            }
        }
    }

    return 0;
}

GADGETRON_WRITER_FACTORY_DECLARE(DicomImageWriter)

} /* namespace Gadgetron */
