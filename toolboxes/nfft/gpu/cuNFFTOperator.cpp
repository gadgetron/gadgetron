#include "cuNDArray_math.h"
#include "cuNFFT.h"
#include "../NFFTOperator.hpp"

namespace Gadgetron{
    template class NFFTOperator<cuNDArray,float,1>;
    template class NFFTOperator<cuNDArray,float,2>;
    template class NFFTOperator<cuNDArray,float,3>;
    template class NFFTOperator<cuNDArray,float,4>;

    template class NFFTOperator<cuNDArray,double,1>;
    template class NFFTOperator<cuNDArray,double,2>;
    template class NFFTOperator<cuNDArray,double,3>;
    template class NFFTOperator<cuNDArray,double,4>;


}
