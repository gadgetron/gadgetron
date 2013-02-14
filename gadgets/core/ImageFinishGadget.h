#ifndef IMAGEFINISHGADGET_H
#define IMAGEFINISHGADGET_H

#include "gadgetroncore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"
#include "GadgetStreamController.h"

#include <complex>
namespace Gadgetron{
template <typename T>
class EXPORTGADGETSCORE ImageFinishGadget : 
public Gadget2<ISMRMRD::ImageHeader,hoNDArray< T > >
{
 protected:
  virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, 
                      GadgetContainerMessage< hoNDArray< T > >* m2);
};

class EXPORTGADGETSCORE ImageFinishGadgetUSHORT :
public ImageFinishGadget<ACE_UINT16>
{
 public:
  GADGET_DECLARE(ImageFinishGadgetUSHORT);
};

class EXPORTGADGETSCORE ImageFinishGadgetFLOAT :
public ImageFinishGadget<float>
{
 public:
  GADGET_DECLARE(ImageFinishGadgetFLOAT);
};

class EXPORTGADGETSCORE ImageFinishGadgetCPLX :
public ImageFinishGadget< std::complex<float> >
{
 public:
  GADGET_DECLARE(ImageFinishGadgetCPLX);
};
}
#endif //IMAGEFINISHGADGET_H
