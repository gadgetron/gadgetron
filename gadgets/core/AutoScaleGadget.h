#ifndef AUTOSCALEGADGET_H_
#define AUTOSCALEGADGET_H_

#include "gadgetron_core_export.h"
#include "Gadget.h"
#include "ismrmrd.h"
#include "hoNDArray.h"

namespace Gadgetron{

  class EXPORTGADGETSCORE AutoScaleGadget:
    public Gadget2<ISMRMRD::ImageHeader,hoNDArray< float > >
  {
  public:
    GADGET_DECLARE(AutoScaleGadget);

    AutoScaleGadget();
    virtual ~AutoScaleGadget();

  protected:
    virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
			GadgetContainerMessage< hoNDArray< float > >* m2);

    unsigned int histogram_bins_;
    std::vector<unsigned int> histogram_;
    float current_scale_;
    float max_value_;
  };
}

#endif /* AUTOSCALEGADGET_H_ */
