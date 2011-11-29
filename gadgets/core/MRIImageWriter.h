#ifndef MRIIMAGEWRITER_H
#define MRIIMAGEWRITER_H

#include "gadgetroncore_export.h"
#include "GadgetMessageInterface.h"
#include "GadgetMRIHeaders.h"

#include <complex>

template<typename T> class MRIImageWriter : public GadgetMessageWriter
{
 public:

  virtual int write(ACE_SOCK_Stream* sock, ACE_Message_Block* mb);
};

class EXPORTGADGETSCORE MRIImageWriterUSHORT : public MRIImageWriter<ACE_UINT16>
{
 public:
  GADGETRON_WRITER_DECLARE(GadgetMessageWriterUSHORT);
};

class EXPORTGADGETSCORE MRIImageWriterFLOAT : public MRIImageWriter<float>
{
 public:
  GADGETRON_WRITER_DECLARE(GadgetMessageWriterFLOAT);
};

class EXPORTGADGETSCORE MRIImageWriterCPLX : public MRIImageWriter< std::complex<float> >
{
 public:
  GADGETRON_WRITER_DECLARE(GadgetMessageWriterCPLX);
};

#endif
