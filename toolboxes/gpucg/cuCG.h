#pragma once

#include "cuNDArray.h"
#include "cuCGMatrixOperator.h"
#include "cuCGPreconditioner.h"

#include <memory>
#include <vector>
#include <iostream>

#include <cublas_v2.h>

// Maximum number of matrix operators in each cg solver
static unsigned int const MAX_CG_MATRIX_OPERATORS = 32;

template<class REAL, class T> class cuCG
{
 public:

  enum cuCGOutputModes {
    OUTPUT_SILENT   = 0,
    OUTPUT_WARNINGS = 1,
    OUTPUT_VERBOSE  = 2,
    OUTPUT_MAX      = 3
  };

  cuCG()
    : num_operators_(0)
    , iterations_(10)
    , limit_(1e-3)
    , output_mode_(OUTPUT_SILENT)
  {
    if (cublasCreate(&cublas_handle_) != CUBLAS_STATUS_SUCCESS) {
      std::cerr << "cuCG unable to create cublas handle" << std::endl;
    }
  }

  virtual ~cuCG() {
    if (cublasDestroy(cublas_handle_) != CUBLAS_STATUS_SUCCESS) {
      std::cerr << "cuCG unable to create cublas handle" << std::endl;
    }
  }
  
  int add_matrix_operator( std::auto_ptr< cuCGMatrixOperator<T> > op, REAL weight )
  {
    operators_[num_operators_++] = op;
    op_weights_.push_back(weight);
    return 0;
  }

  int set_preconditioner( std::auto_ptr< cuCGPreconditioner<T> > precond ) {
    precond_ = precond;
    return 0;
  }

  void set_output_mode(int output_mode) {
    if (!(output_mode >= OUTPUT_MAX || output_mode < 0)) {
      output_mode_ = output_mode;
    }
  }

  void set_limit(REAL limit) {
    limit_ = limit;
  }

  void set_iterations(unsigned int iterations) {
    iterations_ = iterations;
  }
  
  cublasHandle_t get_cublas_handle() {
    return cublas_handle_; 
  }

  std::auto_ptr< cuNDArray<T> > solve(cuNDArray<T>* rhs);

 protected:
  // It is not safe to use auto_ptrs in std STL containers, so we store them in an array instead.
  std::auto_ptr< cuCGMatrixOperator<T> > operators_[MAX_CG_MATRIX_OPERATORS];
  unsigned int num_operators_;
  std::vector<REAL> op_weights_;
  std::auto_ptr< cuCGPreconditioner<T> > precond_;
  cublasHandle_t cublas_handle_;
  unsigned int iterations_;
  REAL limit_;
  int output_mode_;
};
