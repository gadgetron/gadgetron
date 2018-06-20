//
// Created by dchansen on 6/16/18.
//
#pragma once

#include <complex>
#include "fatwater.h"
#include "hoNDArray.h"
namespace Gadgetron{
    namespace FatWater {

        void EXPORTFATWATER correct_frequency_shift(hoNDArray<std::complex<float>> &species_images, const Parameters &params);


    }
}
