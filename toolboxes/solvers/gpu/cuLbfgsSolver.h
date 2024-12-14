#pragma once

#include "lbfgsSolver.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"

#include <fstream>
#include "cuSolverUtils.h"

namespace Gadgetron{

  template <class T> class cuLbfgsSolver : public lbfgsSolver<cuNDArray<T> >
  {
  public:

    cuLbfgsSolver() : lbfgsSolver<cuNDArray<T> >() {}
    virtual ~cuLbfgsSolver() {}
/*
    virtual void iteration_callback(cuNDArray<T>* x ,int iteration,typename realType<T>::Type value){
  	  if (iteration == 0){
  		  std::ofstream textFile("residual.txt",std::ios::trunc);
  	  	  textFile << value << std::endl;
  	  } else{
  		  std::ofstream textFile("residual.txt",std::ios::app);
  		  textFile << value << std::endl;
  	  }

    };
    */
  };
}
