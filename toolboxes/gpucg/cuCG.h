#ifndef CUCG_H
#define CUCG_H

#include <vector>

#include "cuNDArray.h"
#include "cuCGMatrixOperator.h"

template <class T> class cuCG
{

 public:
  cuCG() {}
  virtual ~cuCG() {}
  
  int add_matrix_operator(cuCGMatrixOperator<T>* op, double weight)
  {
    operators_.push_back(op);
    weights_.push_back(weight);
    return 0;
  }

  cuNDArray<T> solve(cuNDArray<T>* rhs);


 protected:
  std::vector< cuCGMatrixOperator<T>* > operators_;
  std::vector<double> weights_;


};

#endif //CUCG_H
