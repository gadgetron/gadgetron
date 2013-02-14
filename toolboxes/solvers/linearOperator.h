#pragma once

#include "vector_td_utilities.h"
#include "solvers_export.h"

#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include "complext.h"

namespace Gadgetron{
template < class ARRAY_TYPE> class linearOperator
{
private:
		typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
	  typedef typename realType<ELEMENT_TYPE>::type REAL;
 public:

  linearOperator() { weight_ = REAL(1); }
  linearOperator(std::vector<unsigned int> *dims) { weight_ = REAL(1);
  	  set_domain_dimensions(dims);
  	  set_codomain_dimensions(dims);
  }
  linearOperator(std::vector<unsigned int> *dims, std::vector<unsigned int> *codims) {
	  weight_ = REAL(1);
	  set_domain_dimensions(dims);
	  set_codomain_dimensions(codims);
  }
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

  virtual void mult_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;
  virtual void mult_MH( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false) = 0;

  virtual void mult_MH_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false )
  {    
    if( codomain_dims_.size() == 0 ){
      throw std::runtime_error("Error: linearOperator::mult_MH_M : codomain dimensions not set");

    }

    ARRAY_TYPE tmp;
    tmp.create(&codomain_dims_);
    mult_M( in, &tmp, false );
    mult_MH( &tmp, out, accumulate );

  }
  
  virtual boost::shared_ptr< linearOperator<ARRAY_TYPE > > clone() = 0;

  void* operator new (size_t bytes) { return ::new char[bytes]; }
  void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
  void * operator new(size_t s, void * p) { return p; }

protected:

  // The template below is useful for implementing the pure virtual 'clone' method 
  //
  
  template <class T>
  boost::shared_ptr<T> clone( T *orig )
  {
    boost::shared_ptr<T> copy( new T() );
    if( !copy.get() ) return boost::shared_ptr<T>();
    *copy = *orig;
    return copy;
  } 
  
protected:
  REAL weight_;
  std::vector<unsigned int> domain_dims_;
  std::vector<unsigned int> codomain_dims_; 
};
}
