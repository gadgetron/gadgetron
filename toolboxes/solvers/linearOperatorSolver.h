/** \file linearOperatorSolver.h
    \brief Base class for all of Gadgetron's solvers operating on linear operators.
*/

#pragma once

#include "solver.h"
#include "linearOperator.h"

#include <vector>
#include <iostream>
#include <stdexcept>

namespace Gadgetron{

  template <class ARRAY_TYPE> class linearOperatorSolver : public solver<ARRAY_TYPE, ARRAY_TYPE>
  {
    
  public:

    // Constructor
    linearOperatorSolver() : solver<ARRAY_TYPE,ARRAY_TYPE>() {}
  
    // Destructor
    virtual ~linearOperatorSolver() {}

    // Add encoding operator to solver (only one allowed)
    virtual void set_encoding_operator( boost::shared_ptr< linearOperator<ARRAY_TYPE> > op)
    {
      if( !op.get() ){
        throw std::runtime_error( "Error: linearOperatorSolver::set_encoding_operator : NULL operator provided" );
      }     
      encoding_operator_ = op;    
    }
  
    virtual boost::shared_ptr< linearOperator<ARRAY_TYPE> >
    get_encoding_operator()
    {
      return encoding_operator_;
    }  
  
    // Add linear operator to solver (in addition to the encoding operator)
    virtual void add_regularization_operator( boost::shared_ptr< linearOperator< ARRAY_TYPE> > op)
    {
      if( !op.get() ){
        throw std::runtime_error( "Error: linearOperatorSolver::add_regularization_operator : NULL operator provided" );
      }    
      regularization_operators_.push_back(op);
    }
  
    virtual boost::shared_ptr< linearOperator< ARRAY_TYPE> >
    get_regularization_operator( unsigned int i )
    {
      if( i >= get_number_of_regularization_operators() ){
        throw std::runtime_error( "Error: linearOperatorSolver::get_regularization_operator : index out of range" );
      }    
      return regularization_operators_[i];
    }  
  
    virtual unsigned int get_number_of_regularization_operators()
    {
      return regularization_operators_.size();
    }
    
  protected:
  
    // Single encoding operator
    boost::shared_ptr< linearOperator<ARRAY_TYPE> > encoding_operator_;
  
    // Vector of linear regularization operators
    std::vector< boost::shared_ptr< linearOperator<ARRAY_TYPE> > > regularization_operators_;
    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;
  };
}
