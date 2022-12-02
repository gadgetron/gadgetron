#ifndef IMAGEWRITERGADGET_H
#define IMAGEWRITERGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_debugging_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{
  
  template <typename T> class ImageWriterGadget :
  public Gadget2<ISMRMRD::ImageHeader, hoNDArray< T > >
  {
    public:
      
    ImageWriterGadget()
      : calls_(0)
	{}
      
    protected:
      virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
			   GadgetContainerMessage< hoNDArray< T > >* m2);
      
      long calls_;      
  };
  
  class EXPORTGADGETSDEBUGGING ImageWriterGadgetUSHORT :
  public ImageWriterGadget<uint16_t >
  {
  public:
    GADGET_DECLARE(ImageWriterGadgetUSHORT)
  };

  class EXPORTGADGETSDEBUGGING ImageWriterGadgetFLOAT :
  public ImageWriterGadget<float>
  {
  public:
    GADGET_DECLARE(ImageWriterGadgetFLOAT)
  };

  class EXPORTGADGETSDEBUGGING ImageWriterGadgetCPLX :
  public ImageWriterGadget< std::complex<float> >
  {
  public:
    GADGET_DECLARE(ImageWriterGadgetCPLX)
  };
}
#endif //IMAGEWRITERGADGET_H
