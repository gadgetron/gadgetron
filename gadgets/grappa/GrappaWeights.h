#pragma once 

#include <ace/Synch.h>

#include <complex>

#include "hoNDArray.h"

template <class T> class GrappaWeights
{
 public:
  GrappaWeights() {}
  virtual ~GrappaWeights() {}
  
  int update(hoNDArray< std::complex<T> >* new_weights);

  int apply(hoNDArray< std::complex<T> >* data_in,
	    hoNDArray< std::complex<T> >* data_out);

 protected:
  ACE_Thread_Mutex mutex_;
  hoNDArray< std::complex<T> > weights_;

};
