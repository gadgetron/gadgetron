#ifndef MRIACQUISITIONREADER_H
#define MRIACQUISITIONREADER_H

#include "gadgetroncore_export.h"
#include "GadgetMessageInterface.h"

class EXPORTGADGETSCORE MRIAcquisitionReader : public GadgetMessageReader
{

 public:
  GADGETRON_READER_DECLARE(MRIAcquisitionReader);
  

  virtual ACE_Message_Block* read(ACE_SOCK_Stream* socket);

};

#endif //MRIACQUISITIONREADER_H
