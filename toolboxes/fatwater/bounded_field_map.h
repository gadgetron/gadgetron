#pragma once
#include "fatwater.h"


namespace Gadgetron {
    namespace FatWater {
            void bounded_field_map(hoNDArray<float> &field_map,
                                        const hoNDArray<std::complex<float>> &input_data,
                                        const Parameters &parameters, float delta_field);
    }
}


