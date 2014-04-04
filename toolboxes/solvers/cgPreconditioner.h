/** \file cgPreconditioner.h
    \brief Base class for preconditioners for the cgSolver class.
*/

#ifndef CGPRECONDITIONER_H
#define CGPRECONDITIONER_H
#pragma once

#include <boost/shared_ptr.hpp>

namespace Gadgetron{

  template <class ARRAY_TYPE> class cgPreconditioner
  {
  public:
    
    cgPreconditioner() {}
    virtual ~cgPreconditioner() {}
    
    virtual void set_weights( boost::shared_ptr<ARRAY_TYPE> w ){
      weights_ = w;
    }

    virtual void apply( ARRAY_TYPE *in, ARRAY_TYPE *out )
    {
      if( !weights_.get() ){
	throw std::runtime_error( "cgPreconditioner::apply(): weights not set");
      }
      
      if ( !in || !out || in->get_number_of_elements() != out->get_number_of_elements()) {
	throw std::runtime_error("cgPreconditioner::apply(): input and output dimensions mismatch");
      }
      
      if (in->get_number_of_elements() % weights_->get_number_of_elements()) {
	throw std::runtime_error( "cgPreconditioner::apply(): unexpected dimensionality of computed weights" );
      }
      *out = *in;
      *out *= *weights_;
    };

  protected:
    boost::shared_ptr<ARRAY_TYPE> weights_;    
  };
}

#endif //CGPRECONDITIONER_H
