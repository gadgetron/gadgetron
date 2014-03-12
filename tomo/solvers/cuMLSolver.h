#pragma once

#include "mlSolver.h"

namespace Gadgetron{
template<class T> class cuMLSolver : public mlSolver<cuNDArray<T> >{
public:
	cuMLSolver() : mlSolver<cuNDArray<T> >() {}
	virtual ~cuMLSolver(){};

};
}
