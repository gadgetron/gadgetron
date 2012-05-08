#pragma once

#include "vector_td_utilities.h"
#include "solvers_export.h"

#include <boost/shared_ptr.hpp>

template <class REAL, class ARRAY_TYPE> class linearOperator
{
 public:

  linearOperator() { weight_ = REAL(1); }
  virtual ~linearOperator() {}

  virtual void set_weight( REAL weight ){ weight_ = weight; }
  virtual REAL get_weight(){ return weight_; }

  virtual bool set_domain_dimensions( std::vector<unsigned int> *dims ) 
  { 
    if( dims == 0x0 ) return false;
    domain_dims_ = *dims; 
    return true;
  }  

  virtual bool set_codomain_dimensions( std::vector<unsigned int> *dims ) 
  { 
    if( dims == 0x0 ) return false;
    codomain_dims_ = *dims; 
    return true;
  }
  
  virtual boost::shared_ptr< std::vector<unsigned int> > get_domain_dimensions() 
  { 
    std::vector<unsigned int> *dims = new std::vector<unsigned int>();
    *dims = domain_dims_; 
    return boost::shared_ptr< std::vector<unsigned int> >(dims);
  }

  virtual boost::shared_ptr< std::vector<unsigned int> > get_codomain_dimensions() 
  { 
    std::vector<unsigned int> *dims = new std::vector<unsigned int>();
    *dims = codomain_dims_; 
    return boost::shared_ptr< std::vector<unsigned int> >(dims);
  }

  virtual int mult_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;
  virtual int mult_MH( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;
  virtual int mult_MH_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;
  
  virtual boost::shared_ptr< linearOperator< REAL, ARRAY_TYPE > > clone() = 0;

  void* operator new (size_t bytes) { return ::new char[bytes]; }
  void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
  void * operator new(size_t s, void * p) { return p; }

protected:

  // The template below is useful for implementing the pure virtual 'clone' method 
  // To be used in _all_ classes that can be instantiated (i.e. non-abstract classes).
  //
  
  template <class T>
  boost::shared_ptr<T> clone( T *orig )
  {
    boost::shared_ptr<T> copy( new T() );
    if( !copy.get() ) return boost::shared_ptr<T>();
    *copy = *orig;
    return copy;
  } 
  
private:
  REAL weight_;
  std::vector<unsigned int> domain_dims_;
  std::vector<unsigned int> codomain_dims_; 
};
