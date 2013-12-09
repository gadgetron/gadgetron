/** \file generalOperator.h
    \brief Base class for all operators on which we can compute a gradient.
*/

#pragma once


#include "complext.h"

#include <boost/shared_ptr.hpp>
#include <vector>
#include <stdexcept>

namespace Gadgetron{

  template <class ARRAY> class generalOperator
  {
   public:

    typedef typename ARRAY::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;

    generalOperator() : weight_(REAL(1)){};

    generalOperator(std::vector<unsigned long long> *dims) : weight_(REAL(1)){
      set_domain_dimensions(dims);
    }
    
    virtual ~generalOperator(){};

    virtual void gradient(ARRAY* in, ARRAY* out, bool accumulate = false ) = 0;
    
    virtual void set_domain_dimensions( std::vector<unsigned long long> *dims )
    {
      if( dims == 0x0 ) throw std::runtime_error("Null pointer provided");
      domain_dims_ = *dims;  
    }
    
    virtual boost::shared_ptr< std::vector<unsigned long long> > get_domain_dimensions()
    {
      std::vector<unsigned long long> *dims = new std::vector<unsigned long long>();
      *dims = domain_dims_;
      return boost::shared_ptr< std::vector<unsigned long long> >(dims);
    }
    
    virtual void set_weight( REAL weight ){ weight_ = weight; }
    virtual REAL get_weight(){ return weight_; }
    
    void* operator new (size_t bytes) { return ::new char[bytes]; }
    void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
    void * operator new(size_t s, void * p) { return p; }
    
  protected:
    REAL weight_;
    std::vector<unsigned long long> domain_dims_;
  };  
}
