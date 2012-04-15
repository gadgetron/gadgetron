#pragma once

#include "solver.h"
#include "matrixOperator.h"

#include <vector>
#include <iostream>

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class matrixOpSolver 
  : public solver<ARRAY_TYPE, ARRAY_TYPE>
{

public:

  // Constructor
  matrixOpSolver() : solver<ARRAY_TYPE,ARRAY_TYPE>() {}
  
  // Destructor
  virtual ~matrixOpSolver() {}

  // Add matrix operator to the solver
  // ---------------------------------
  // The latter two arguments are used during subsequent calls to 'solve' and 'solve_from_rhs'
  //
  // When using the 'solve_from_rhs' interface, 'contributes_to_rhs' and 'rhs_data' are ignored.
  // When using the 'solve' interface
  // - 'contributions_to_rhs' indicates if this operator contributes to the right hand side (rhs):
  // - if true, the adjoint matrix operator (op.mult_MH) is computed on 'rhs_data' during the rhs computation
  // - if true and 'rhs_data' is 0x0, (op.mult_MH) is computed on the input data to the 'solve' method 

  inline bool add_matrix_operator( boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > op,
				   bool contributes_to_rhs = false, 
				   boost::shared_ptr<ARRAY_TYPE> rhs_data = boost::shared_ptr<ARRAY_TYPE>() ) 
  {
    if( !op.get() ){
      this->solver_error( "Error: matrixOpSolver::add_matrix_operator : NULL operator provided" );
      return false;
    } 
    
    operators_.push_back(op);
    indicators_.push_back(contributes_to_rhs);
    rhs_data_.push_back(rhs_data);
    
    return true;
  }

protected:
  
  // Vector of matrix operators
  std::vector< boost::shared_ptr< matrixOperator<REAL, ARRAY_TYPE> > > operators_;
  
  // Vector of boolean rhs indicators for the operators
  std::vector<bool> indicators_;
  
  // Vector of rhs priors for the operators
  std::vector< boost::shared_ptr<ARRAY_TYPE> > rhs_data_;  
};
