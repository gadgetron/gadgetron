#pragma once

#include "gpBbSolver.h"
#include "hoNDArray_math.h"
#include "real_utilities.h"
#include "vector_td_utilities.h"



namespace Gadgetron{

  template <class T> class hoGpBbSolver : public gpBbSolver< hoNDArray<T> >
  {  
  public:

    hoGpBbSolver() : gpBbSolver< hoNDArray<T> >() {};
    virtual ~hoGpBbSolver() {};
        

  };
}
