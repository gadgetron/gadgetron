#ifndef MRIACQUISITIONREADER_H
#define MRIACQUISITIONREADER_H

#include "GadgetMessageInterface.h"

class MRIAcquisitionReader : public GadgetMessageReader
{

 public:
  GADGETRON_READER_DECLARE(MRIAcquisitionReader);
  

  virtual ACE_Message_Block* read(ACE_SOCK_Stream* socket);

};

#endif //MRIACQUISITIONREADER_H
