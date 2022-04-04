/**
    \brief  Performs fft on complex image
    \author Original: Hui Xue
    \author PureGadget Conversion: Andrew Dupuis
    \test   Untested
*/

#pragma once

#include "PureGadget.h"
#include "Types.h"
#include "hoNDArray_math.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"
#include "hoNDFFT.h"
#include "mri_core_def.h"

namespace Gadgetron{

    using ComplexImage = Core::variant<Core::Image<std::complex<float>>, Core::Image<std::complex<double>>>;

    class ImageFFTGadget : public Core::PureGadget<ComplexImage, ComplexImage> {
    public:
        using Core::PureGadget<ComplexImage,ComplexImage>::PureGadget;
        ComplexImage process_function(ComplexImage image) const override; 
    };
}
