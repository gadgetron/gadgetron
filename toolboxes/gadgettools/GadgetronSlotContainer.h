#ifndef GADGETRONSLOTCONTAINER_H_
#define GADGETRONSLOTCONTAINER_H_

#include <algorithm>
#include <vector>

template <typename T> class GadgetronSlotContainer {

public:
	GadgetronSlotContainer() {}

	virtual ~GadgetronSlotContainer()
	{
		clear();
	}

	T* find(unsigned int slot) {
	    T* ret = 0;
	    for (unsigned int i = 0; i < slots_.size(); i++) {
	    	if (slots_[i] == slot) {
	    		ret = items_[i];
	    		break;
	    	}
	    }
	    return ret;
	  }

	  int insert ( unsigned short slot, T* item) {
		  if (this->find(slot)) {
			  return -1;
		  } else {
			  slots_.push_back(slot);
			  items_.push_back(item);
		  }
		  return 0;
	  }

	  int clear()
	  {
		  for (unsigned int i = 0; i < items_.size(); i++) {
			  if (items_[i]) delete items_[i];
		  }
		  slots_.clear();
		  items_.clear();
		  return 0;
	  }

protected:
	std::vector<unsigned int> slots_;
	std::vector<T*> items_;
};

#endif /* GADGETRONSLOTCONTAINER_H_ */
