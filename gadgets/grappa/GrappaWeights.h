#pragma once 

#include "gadgetron_grappa_export.h"
#include "hoNDArray.h"

#include <ace/Synch.h>
#include <complex>

namespace Gadgetron{

template <class T> class EXPORTGADGETSGRAPPA GrappaWeights
{
 public:
  GrappaWeights()
  	  : weights_are_valid_(false)
  	  , cond_(cond_mutex_)
  	  {

  	  }
  virtual ~GrappaWeights() {}
  
  int update(hoNDArray< std::complex<T> >* new_weights);

  int apply(hoNDArray< std::complex<T> >* data_in,
	    hoNDArray< std::complex<T> >* data_out, 
	    T scale = 1.0);

 private:
  ACE_Thread_Mutex mutex_;
  bool weights_are_valid_;

  ACE_Thread_Mutex cond_mutex_;
  ACE_Condition_Thread_Mutex cond_;
  hoNDArray< std::complex<T> > weights_;

};
}
