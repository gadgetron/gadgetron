#pragma once

#include "hoNDArray_math.h"
#include "hoCuNDArray_math.h"
#include "hoNDArray_fileio.h"
#include "complext.h"
#include "nlcgSolver.h"
#include "hoSolverUtils.h"

namespace Gadgetron{

template<class T> class hoCuNlcgSolver: public nlcgSolver<hoCuNDArray<T> >{
	typedef typename realType<T>::Type REAL;
public:
	hoCuNlcgSolver():nlcgSolver<hoCuNDArray<T> >(){

	}

	virtual ~hoCuNlcgSolver(){};

  virtual void iteration_callback(hoCuNDArray<T>* x,int i,REAL data_res,REAL reg_res){
	  /*
	  if (i == 0){
		  std::ofstream textFile("residual.txt",std::ios::trunc);
	  	  textFile << data_res << std::endl;
	  } else{
		  std::ofstream textFile("residual.txt",std::ios::app);
		  textFile << data_res << std::endl;
	  }
	  std::stringstream ss;
	  ss << "iteration-" << i << ".real";
	  write_nd_array(x,ss.str().c_str());*/
  };
};
}
