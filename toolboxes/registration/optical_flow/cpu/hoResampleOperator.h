#pragma once

#include "resampleOperator.h"
#include "hoNDArray.h"
#include "gpureg_export.h"

template <class REAL, class T, unsigned int D>
class EXPORTGPUREG hoResampleOperator : public resampleOperator< REAL, hoNDArray<REAL>, hoNDArray<T> >
{
  
 public:
  
  hoResampleOperator() : resampleOperator< REAL, hoNDArray<REAL>, hoNDArray<T> >() {}
  virtual ~hoResampleOperator() {}
  
  // Todo: probably move this function to the Cuda (cu) branch
  virtual bool mult_MH_preprocess() 
  { 
    if(this->offsets_.get()){
      this->preprocessed_ = true; 
      return true;
    }
    else{
      this->preprocessed_ = false; 
      return false;
    }
  }
};
