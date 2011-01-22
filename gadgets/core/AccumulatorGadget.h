#ifndef ACCUMULATORGADGET_H
#define ACCUMULATORGADGET_H

#include <complex>

#include "Gadget.h"
#include "gadgetheaders.h"
#include "NDArray.h"

class AccumulatorGadget : 
public Gadget2< GadgetMessageAcquisition, NDArray< std::complex<float> > >
{
  
 public:
  GADGET_DECLARE(AccumulatorGadget);

  AccumulatorGadget();
  ~AccumulatorGadget();

 protected:
  virtual int process_config(ACE_Message_Block* mb);
  virtual int process(GadgetContainerMessage< GadgetMessageAcquisition >* m1,
		      GadgetContainerMessage< NDArray< std::complex<float> > > * m2);

  NDArray< std::complex<float> >* buffer_;
  std::vector<int> dimensions_;

};

#endif //ACCUMULATORGADGET_H
