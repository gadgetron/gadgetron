#pragma once

#include <cpu/hoNDArray.h>
#include "fatwater.h"

namespace Gadgetron {

    namespace FatWater {
        void EXPORTFATWATER fat_water_fitting(hoNDArray<float> &field_map, hoNDArray<float> &r2star_map,
                               hoNDArray<std::complex<float>> &fractions,
                               const hoNDArray<std::complex<float>> &input_data,
                               const hoNDArray<float> &lambda_map, const Parameters &parameters );
    }

}