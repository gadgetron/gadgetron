#ifndef MATRIXOPERATOR_H
#define MATRIXOPERATOR_H
#pragma once

#include "vector_td_utilities.h"
#include "solvers_export.h"

template <class REAL, class ARRAY_TYPE> class matrixOperator
{
 public:

  matrixOperator() { weight_ = REAL(1); }
  virtual ~matrixOperator() {}

  inline void set_weight( REAL weight ){ weight_ = weight; }
  inline REAL get_weight(){ return weight_; }

  inline void set_domain_dimensions( std::vector<unsigned int> dims ) { domain_dims_ = dims; }  
  inline void set_codomain_dimensions( std::vector<unsigned int> dims ) { codomain_dims_ = dims; }
  
  inline std::vector<unsigned int> get_domain_dimensions() { return domain_dims_; }
  inline std::vector<unsigned int> get_codomain_dimensions() { return codomain_dims_; }

  virtual int mult_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;
  virtual int mult_MH( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;
  virtual int mult_MH_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;
  
  void* operator new (size_t bytes) { return ::new char[bytes]; }
  void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
  void * operator new(size_t s, void * p) { return p; }

private:
  REAL weight_;
  std::vector<unsigned int> domain_dims_;
  std::vector<unsigned int> codomain_dims_; 
};

#endif //MATRIXOPERATOR_H
