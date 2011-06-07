#ifndef MRIIMAGEWRITER_H
#define MRIIMAGEWRITER_H

#include "gadgetron_export.h"
#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"

class EXPORTGADGETSCORE MRIImageWriter : public GadgetMessageWriter
{
 public:
  GADGETRON_WRITER_DECLARE(GadgetMessageWriter);

  virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);
};

#endif
