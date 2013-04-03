#pragma once

#include "linearOperator.h"

template <class ARRAY_TYPE_REAL, class ARRAY_TYPE_ELEMENT> 
class resampleOperator : public linearOperator<ARRAY_TYPE_ELEMENT>
{
public:
  
  resampleOperator() : linearOperator<ARRAY_TYPE_ELEMENT>() { preprocessed_ = false; }
  virtual ~resampleOperator() {}
  
  // Expected format: the vector field dimension should be the slowest varying
  //

  virtual void set_displacement_field( boost::shared_ptr<ARRAY_TYPE_REAL> offsets )
  {
    offsets_ = offsets;
  }
  
  virtual boost::shared_ptr<ARRAY_TYPE_REAL> get_displacement_field()
  {
    return offsets_;
  }

  virtual unsigned int get_number_of_displacement_vectors() 
  {
    if( !offsets_.get() ) return 0;
    offsets_->get_number_of_elements()/offsets_->get_size(offsets_->get_number_of_dimensions()-1);
  }
    
protected:
  boost::shared_ptr<ARRAY_TYPE_REAL> offsets_;
};
