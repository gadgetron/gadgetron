#ifndef AUTOSCALEGADGET_H_
#define AUTOSCALEGADGET_H_

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>

namespace Gadgetron{

  class EXPORTGADGETSMRICORE AutoScaleGadget:
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
    std::vector<size_t> histogram_;
    float current_scale_;
    float max_value_;
  };
}

#endif /* AUTOSCALEGADGET_H_ */
