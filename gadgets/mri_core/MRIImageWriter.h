#ifndef MRIIMAGEWRITER_H
#define MRIIMAGEWRITER_H

#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{

  template<typename T> class MRIImageWriter : public GadgetMessageWriter
  {
  public:
    virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);
  };

  class EXPORTGADGETSMRICORE MRIImageWriterUSHORT : public MRIImageWriter<ACE_UINT16>
  {
  public:
    GADGETRON_WRITER_DECLARE(GadgetMessageWriterUSHORT);
  };

  class EXPORTGADGETSMRICORE MRIImageWriterFLOAT : public MRIImageWriter<float>
  {
  public:
    GADGETRON_WRITER_DECLARE(GadgetMessageWriterFLOAT);
  };

  class EXPORTGADGETSMRICORE MRIImageWriterCPLX : public MRIImageWriter< std::complex<float> >
  {
  public:
    GADGETRON_WRITER_DECLARE(GadgetMessageWriterCPLX);
  };
}
#endif
