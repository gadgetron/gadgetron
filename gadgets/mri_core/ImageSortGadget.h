#ifndef IMAGESORTGADGET_H
#define IMAGESORTGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron{

  struct ImageEntry
  {
    int index_;
    GadgetContainerMessage<ISMRMRD::ImageHeader>* mb_;
  };

  bool image_entry_compare(const ImageEntry& i, const ImageEntry& j)
  {
    return (i.index_<j.index_);
  }

  class EXPORTGADGETSMRICORE ImageSortGadget : public Gadget1 < ISMRMRD::ImageHeader >
  {
  public:
    GADGET_DECLARE(ImageSortGadget);

  protected:
    GADGET_PROPERTY_LIMITS(sorting_dimension, std::string, "Dimension that data will be sorted by", "slice",
			   GadgetPropertyLimitsEnumeration, 
			   "average",
			   "slice",
			   "contrast",
			   "phase",
			   "repetition",
			   "set");

    virtual int close(unsigned long flags);
    virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);
    int index(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1);
    
    std::vector<ImageEntry> images_;
    
  };
}

#endif //IMAGESORTGADGET_H
