#ifndef FlowPhaseSubtractionGadget_H
#define FlowPhaseSubtractionGadget_H

#include <complex>
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "spiral_gadgets_export.h"
#include "ismrmrd.h"


class EXPORTGADGETSSPIRAL FlowPhaseSubtractionGadget :
public Gadget2< GadgetMessageImage, hoNDArray< std::complex<float> > >
{
  
 public:
  GADGET_DECLARE(FlowPhaseSubtractionGadget);

  FlowPhaseSubtractionGadget();
  virtual ~FlowPhaseSubtractionGadget();

 protected:
  virtual int process_config(ACE_Message_Block* mb);
  virtual int process(GadgetContainerMessage< GadgetMessageImage >* m1,
		      GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);

 private:
  unsigned int sets_;
  ACE_Message_Queue<ACE_MT_SYNCH> buffer_;

};

#endif //FlowPhaseSubtractionGadget_H
