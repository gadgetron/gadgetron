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
      this->solver_error( "Error: linearSolver::add_matrix_operator : NULL operator provided" );
      return false;
    } 
    
    encoding_operator_ = op;
    
    return true;
  }
  
  // Add linear operator to solver (in addition to the encoding operator)
  virtual bool add_regularization_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > op)
  {
    if( !op.get() ){
      this->solver_error( "Error: linearSolver::add_matrix_operator : NULL operator provided" );
      return false;
    }
    
    regularization_operators_.push_back(op);
    
    return true;
  }
  
protected:
  
  // Single encoding operator
  boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > encoding_operator_;
  
  // Vector of linear regularization operators
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > > regularization_operators_;
};
