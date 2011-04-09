#ifndef CUCG_H
#define CUCG_H

#include <vector>
#include <iostream>

#include <cublas_v2.h>

#include "cuNDArray.h"
#include "cuCGMatrixOperator.h"
#include "cuCGPreconditioner.h"

template <class T> class cuCG
{
 public:
  cuCG() 
    : precond_(0) 
    , iterations_(10)
    , limit_(1e-3)
    , output_mode_(0)
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
  
  int add_matrix_operator(cuCGMatrixOperator<T>* op, double weight)
  {
    operators_.push_back(op);
    weights_.push_back(weight);
    return 0;
  }

  int set_preconditioner(cuCGPreconditioner<T>* precond) {
    precond_ = precond;
    return 0;
  }

  cuNDArray<T> solve(cuNDArray<T>* rhs);


 protected:
  std::vector< cuCGMatrixOperator<T>* > operators_;
  std::vector<double> weights_;
  cuCGPreconditioner<T>* precond_;
  cublasHandle_t cublas_handle_;
  unsigned int iterations_;
  double limit_;
  int output_mode_;

};

#endif //CUCG_H
