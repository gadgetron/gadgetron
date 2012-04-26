#pragma once

#include "solver.h"
#include "linearOperator.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class linearSolver
  : public solver<ARRAY_TYPE, ARRAY_TYPE>
{

public:

  // Constructor
  linearSolver() : solver<ARRAY_TYPE,ARRAY_TYPE>() {}
  
  // Destructor
  virtual ~linearSolver() {}

  // Add encoding operator to solver (only one allowed)
  virtual bool set_encoding_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > op)
  {
    if( !op.get() ){
      this->solver_error( "Error: linearSolver::set_encoding_operator : NULL operator provided" );
      return false;
    } 
    
    encoding_operator_ = op;
    
    return true;
  }
  
  virtual boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > 
  get_encoding_operator()
  {
    return encoding_operator_;
  }  
  
  // Add linear operator to solver (in addition to the encoding operator)
  virtual bool add_regularization_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > op)
  {
    if( !op.get() ){
      this->solver_error( "Error: linearSolver::add_regularization_operator : NULL operator provided" );
      return false;
    }
    
    regularization_operators_.push_back(op);
    
    return true;
  }
  
  virtual boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > 
  get_regularization_operator( unsigned int i )
  {
    if( i >= get_number_of_regularization_operators() ){
      this->solver_error( "Error: linearSolver::get_regularization_operator : index out of range" );
      return boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> >();
    }
    
    return regularization_operators_[i];
  }  

  virtual unsigned int get_number_of_regularization_operators()
  {
    return regularization_operators_.size();
  }
    
protected:
  
  // Single encoding operator
  boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > encoding_operator_;
  
  // Vector of linear regularization operators
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > > regularization_operators_;
};
