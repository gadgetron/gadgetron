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

  // Add matrix operator to the solver
  // ---------------------------------
  // The latter two arguments are used during subsequent calls to 'solve' and 'solve_from_rhs'
  //
  // When using the 'solve_from_rhs' interface, 'contributes_to_rhs' and 'rhs_data' are ignored.
  // When using the 'solve' interface
  // - 'contributions_to_rhs' indicates if this operator contributes to the right hand side (rhs):
  // - if true, the adjoint matrix operator (op.mult_MH) is computed on 'rhs_data' during the rhs computation
  // - if true and 'rhs_data' is 0x0, (op.mult_MH) is computed on the input data to the 'solve' method 

  inline bool add_linear_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > op)
    {
      if( !op.get() ){
        this->solver_error( "Error: linearSolver::add_matrix_operator : NULL operator provided" );
        return false;
      }

      operators_.push_back(op);

      return true;
    }

  inline bool add_encoding_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > op)
  {
    if( !op.get() ){
      this->solver_error( "Error: linearSolver::add_matrix_operator : NULL operator provided" );
      return false;
    } 
    
    encoding_op = op;

    return true;
  }

protected:
  
  // Vector of matrix operators
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > > operators_;
  boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > encoding_op;
  
  
};
