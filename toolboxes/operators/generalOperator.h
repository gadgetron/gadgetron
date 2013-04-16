/** \file generalOperator.h
    \brief Base class for all operators on which we can compute a gradient.
*/

#pragma once

#include "GadgetronException.h"
#include "complext.h"

#include <boost/shared_ptr.hpp>
#include <vector>

namespace Gadgetron{

  template <class ARRAY> class generalOperator
  {
   public:

    typedef typename ARRAY::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;

    generalOperator() : weight_(REAL(1)){};

    generalOperator(std::vector<unsigned int> *dims) : weight_(REAL(1)){
      set_domain_dimensions(dims);
    }
    
    virtual ~generalOperator(){};

    virtual void gradient(ARRAY* in, ARRAY* out, bool accumulate = false ) = 0;
    
    virtual void set_domain_dimensions( std::vector<unsigned int> *dims )
    {
      if( dims == 0x0 ) BOOST_THROW_EXCEPTION(runtime_error("Null pointer provided"));
      domain_dims_ = *dims;  
    }
    
    virtual boost::shared_ptr< std::vector<unsigned int> > get_domain_dimensions()
    {
      std::vector<unsigned int> *dims = new std::vector<unsigned int>();
      *dims = domain_dims_;
      return boost::shared_ptr< std::vector<unsigned int> >(dims);
    }
    
    virtual void set_weight( REAL weight ){ weight_ = weight; }
    virtual REAL get_weight(){ return weight_; }
    
    void* operator new (size_t bytes) { return ::new char[bytes]; }
    void operator delete (void *ptr) { delete [] static_cast <char *> (ptr); } 
    void * operator new(size_t s, void * p) { return p; }
    
  protected:
    REAL weight_;
    std::vector<unsigned int> domain_dims_;
  };  
}
