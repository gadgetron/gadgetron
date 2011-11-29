#ifndef IMAGEWRITERGADGET_H
#define IMAGEWRITERGADGET_H

#include <complex>

#include "gadgetroncore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"


template <typename T> class ImageWriterGadget :
public Gadget2<GadgetMessageImage, hoNDArray< T > >
{
 public:

  ImageWriterGadget()
    : calls_(0)
    {}

 protected:
  virtual int process( GadgetContainerMessage< GadgetMessageImage>* m1,
		       GadgetContainerMessage< hoNDArray< T > >* m2);

  long calls_;

};

class EXPORTGADGETSCORE ImageWriterGadgetUSHORT :
public ImageWriterGadget<ACE_UINT16>
{
 public:
  GADGET_DECLARE(ImageWriterGadgetUSHORT)
};

class EXPORTGADGETSCORE ImageWriterGadgetFLOAT :
public ImageWriterGadget<float>
{
 public:
  GADGET_DECLARE(ImageWriterGadgetFLOAT)
};

class EXPORTGADGETSCORE ImageWriterGadgetCPLX :
public ImageWriterGadget< std::complex<float> >
{
 public:
  GADGET_DECLARE(ImageWriterGadgetCPLX)
};

#endif //IMAGEWRITERGADGET_H
