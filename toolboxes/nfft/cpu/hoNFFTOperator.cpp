#include "hoNDArray_math.h"
#include "hoNFFT.h"
#include "../NFFTOperator.hpp"

namespace Gadgetron{
    template class NFFTOperator<hoNDArray,float,1>;
    template class NFFTOperator<hoNDArray,float,2>;
    template class NFFTOperator<hoNDArray,float,3>;

    template class NFFTOperator<hoNDArray,double,1>;
    template class NFFTOperator<hoNDArray,double,2>;
    template class NFFTOperator<hoNDArray,double,3>;
}

