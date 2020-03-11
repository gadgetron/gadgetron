//
// Created by dchansen on 3/4/20.
//

#pragma once

#include "hoNDArray.h"
namespace Gadgetron::T1 {

    struct T1_2param {
        hoNDArray<float> A;
        hoNDArray<float> T1;
    };

    struct T1_3param {
        hoNDArray<float> A;
        hoNDArray<float> B;
        hoNDArray<float> T1star;
    };

    /**
     * Fits a T1 map using the 2 parameter model
     * @param data Data of shape (X,Y,TI)
     * @param TI Inversion times
     * @return Magnitude (A) and T1 mapping
     */
    T1_2param fit_T1_2param(const hoNDArray<float>& data, const std::vector<float>& TI);

    /**
    * Fits a T1 map using the 2 parameter model
    * @param data Data of shape (X,Y,TI)
    * @param TI Inversion times
    * @return Magnitude (A), inverse magnitude (B) and T1 mapping
    */
    T1_3param fit_T1_3param(const hoNDArray<float>& data, const std::vector<float>& TI);


    T1_3param motion_compensated_t1_fit(const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI, unsigned int iterations=5);



}
