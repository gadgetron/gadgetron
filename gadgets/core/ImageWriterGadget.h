#ifndef IMAGEWRITERGADGET_H
#define IMAGEWRITERGADGET_H

#include <complex>

#include "gadgetroncore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"


class EXPORTGADGETSCORE ImageWriterGadget :
public Gadget2<GadgetMessageImage, hoNDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(ImageWriterGadget)

  ImageWriterGadget()
    : calls_(0)
    {}

 protected:
  virtual int process( GadgetContainerMessage< GadgetMessageImage>* m1,
		       GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

  long calls_;

};

#endif //IMAGEWRITERGADGET_H
