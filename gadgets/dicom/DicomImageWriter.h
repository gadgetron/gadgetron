#ifndef DICOMIMAGEWRITER_H
#define DICOMIMAGEWRITER_H

#include "gadgetron_dicom_export.h"
#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd/ismrmrd.h"


namespace Gadgetron {

class EXPORTGADGETSDICOM DicomImageWriter : public GadgetMessageWriter
{
 public:
  virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);

  GADGETRON_WRITER_DECLARE(DicomImageWriter);
};

class EXPORTGADGETSDICOM DicomImageAttribWriter : public GadgetMessageWriter
{
 public:
  virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);

  GADGETRON_WRITER_DECLARE(DicomImageAttribWriter);
};

} /* namespace Gadgetron */

#endif
