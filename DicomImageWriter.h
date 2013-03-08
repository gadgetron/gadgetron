#ifndef DICOMIMAGEWRITER_H
#define DICOMIMAGEWRITER_H

#include "gadgetroncore_export.h"
#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"


class EXPORTGADGETSCORE DicomImageWriter : public GadgetMessageWriter
{
 public:
  virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);

  GADGETRON_WRITER_DECLARE(DicomImageWriter);
};

/*
class EXPORTGADGETSCORE DicomImageWriterUSHORT : public DicomImageWriter<ACE_UINT16>
{
 public:
  GADGETRON_WRITER_DECLARE(DicomImageWriterUSHORT);
};

class EXPORTGADGETSCORE DicomImageWriterFLOAT : public DicomImageWriter<float>
{
 public:
  GADGETRON_WRITER_DECLARE(DicomImageWriterFLOAT);
};

class EXPORTGADGETSCORE DicomImageWriterCPLX : public DicomImageWriter< std::complex<float> >
{
 public:
  GADGETRON_WRITER_DECLARE(DicomImageWriterCPLX);
};
*/

#endif
