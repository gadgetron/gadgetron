#pragma once

#include "sartSolver.h"

namespace Gadgetron{
template <class T> class cuSARTSolver : public sartSolver<cuNDArray<T> >{
public:
	cuSARTSolver() : sartSolver<cuNDArray<T> >() {}
	virtual ~cuSARTSolver(){};

};
}
