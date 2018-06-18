#pragma once



#include "hoNDArray.h"
namespace  Gadgetron {


    hoNDArray <uint16_t>
    update_field_map(const hoNDArray <uint16_t> &field_map_index, const hoNDArray <uint16_t> &proposed_field_map_index,
                     const hoNDArray<float> &residuals_map, const hoNDArray<float> &lambda_map);

}