#include "gpunfft_export.h"
#include "cuNDArray_math.h"
#include "cuNFFT.h"
#include "../NFFTOperator.hpp"

namespace Gadgetron{
    template EXPORTGPUNFFT class NFFTOperator<cuNDArray,float,1>;
    template EXPORTGPUNFFT class NFFTOperator<cuNDArray,float,2>;
    template EXPORTGPUNFFT class NFFTOperator<cuNDArray,float,3>;
    template EXPORTGPUNFFT class NFFTOperator<cuNDArray,float,4>;

    template EXPORTGPUNFFT class NFFTOperator<cuNDArray,double,1>;
    template EXPORTGPUNFFT class NFFTOperator<cuNDArray,double,2>;
    template EXPORTGPUNFFT class NFFTOperator<cuNDArray,double,3>;
    template EXPORTGPUNFFT class NFFTOperator<cuNDArray,double,4>;


}
