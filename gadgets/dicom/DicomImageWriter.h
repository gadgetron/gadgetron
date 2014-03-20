#ifndef DICOMIMAGEWRITER_H
#define DICOMIMAGEWRITER_H

#include "gadgetron_dicom_export.h"
#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"


namespace Gadgetron {

class EXPORTGADGETSDICOM DicomImageWriter : public GadgetMessageWriter
{
 public:
  virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);
};

} /* namespace Gadgetron */

#endif
