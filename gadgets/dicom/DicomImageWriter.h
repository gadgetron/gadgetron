#ifndef DICOMIMAGEWRITER_H
#define DICOMIMAGEWRITER_H

#include "gadgetron_dicom_export.h"
// DCMTK includes
#include "dcmtk/config/osconfig.h"
#include "dcmtk/ofstd/ofstdinc.h"
#define INCLUDE_CSTDLIB
#define INCLUDE_CSTDIO
#define INCLUDE_CSTRING
#include "dcmtk/dcmdata/dctk.h"
#include "dcmtk/dcmdata/dcostrmb.h"

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
