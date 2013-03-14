#ifndef ACCUMULATORGADGET_H
#define ACCUMULATORGADGET_H

#include <complex>

#include "gadgetroncore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd.h"
namespace Gadgetron{
class EXPORTGADGETSCORE AccumulatorGadget : 
public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
{
  
 public:
  GADGET_DECLARE(AccumulatorGadget);

  AccumulatorGadget();
  ~AccumulatorGadget();

 protected:
  virtual int process_config(ACE_Message_Block* mb);
  virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1,
		      GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);

  hoNDArray< std::complex<float> >* buffer_;
  std::vector<unsigned int> dimensions_;
  std::vector<float> field_of_view_;
  unsigned int slices_;

  int image_counter_;
  int image_series_;

};
}
#endif //ACCUMULATORGADGET_H
