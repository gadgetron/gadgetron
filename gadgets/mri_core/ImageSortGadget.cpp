#include "ImageSortGadget.h"

namespace Gadgetron{

  int ImageSortGadget::index(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
  {
    std::string sorting_dimension_local = sorting_dimension.value();
    
    if (sorting_dimension_local.size() == 0) {
      return -1;
    } else if (sorting_dimension_local.compare("average") == 0) {
      return m1->getObjectPtr()->average;
    } else if (sorting_dimension_local.compare("slice") == 0) {
      return m1->getObjectPtr()->slice;
    } else if (sorting_dimension_local.compare("contrast") == 0) {
      return m1->getObjectPtr()->contrast;
    } else if (sorting_dimension_local.compare("phase") == 0) {
      return m1->getObjectPtr()->phase;
    } else if (sorting_dimension_local.compare("repetition") == 0) {
      return m1->getObjectPtr()->repetition;
    } else if (sorting_dimension_local.compare("set") == 0) {
      return m1->getObjectPtr()->set;
    } else {
      return -1;
    }

    return -1;
  }

  int ImageSortGadget::close(unsigned long flags)
  {
    GDEBUG("++++++ close call with %d images\n", images_.size());
    if (images_.size()) {
      
      std::sort(images_.begin(),images_.end(), image_entry_compare);
      
      for (auto it = images_.begin(); it != images_.end(); it++) {
	if (this->next()->putq(it->mb_) == -1) {
	  it->mb_->release();
	  GERROR("Error passing data on to next gadget\n");
	  return GADGET_FAIL;
	}
      }
      
      images_.clear();
    }
    return GADGET_OK;
  }
  
  int ImageSortGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1)
  {
    if (index(m1) < 0) {
      if (this->next()->putq(m1) == -1) {
	m1->release();
	GERROR("Error passing data on to next gadget\n");
	return GADGET_FAIL;
      }
    }

    ImageEntry i;
    i.index_ = index(m1);
    i.mb_ = m1;

    images_.push_back(i);

    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(ImageSortGadget);
}
