#pragma once 

#include <ace/Synch.h>

#include <complex>

#include "gadgetron_export.h"
#include "hoNDArray.h"

template <class T> class EXPORTGADGETSGRAPPA GrappaWeights
{
 public:
  GrappaWeights() {}
  virtual ~GrappaWeights() {}
  
  int update(hoNDArray< std::complex<T> >* new_weights);

  int apply(hoNDArray< std::complex<T> >* data_in,
	    hoNDArray< std::complex<T> >* data_out, 
	    T scale = 1.0);

 private:
  ACE_Thread_Mutex mutex_;
  hoNDArray< std::complex<T> > weights_;

};
