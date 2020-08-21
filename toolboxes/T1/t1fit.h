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


    struct registration_params {
        unsigned int iterations = 40;
        float regularization_sigma = 2.0f;
        float step_size = 2.0;
        float noise_sigma = 0.0f;
    }
    ;

    /**
     * Performs registration on a T1 dataset by iteratively creating synthetic data based 
     * on a two parameter T1 fit, and registering it to the original data using intensity based registration.
     * @param data Input data of shape (X,Y,TI)
     * @param TI inversion times
     * @param iterations Number of iterations to be used.
     * @return Deformation vector fields for bringing each dataset into a common reference frame
     **/ 
    hoNDArray<vector_td<float,2>> t1_registration(const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI, unsigned int iterations=5, registration_params params = {});
    hoNDArray<vector_td<float, 2>> multi_scale_t1_registration(   const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI, unsigned int levels=3, unsigned int iterations=5, registration_params params = {});

    hoNDArray<std::complex<float>> deform_groups(const hoNDArray<std::complex<float>>& data,const hoNDArray<vector_td<float,2>>& vector_field);

    hoNDArray<float> predict_signal(const T1_2param& params, const std::vector<float>& TI);

    hoNDArray<float> predict_signal(const T1_3param& params, const std::vector<float>& TI);

    hoNDArray<float> phase_correct(const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI);

    hoNDArray<vector_td<float,2>> register_groups_CMR(const hoNDArray<float>& phase_corrected_data,
                                                  const hoNDArray<float>& predicted );
    /**
     * Preforms registration on a T1 dataset by iteratively creating synthetic data. Uses the CMR motion correction code, and returns the deformed data
     * @param data Input data of shape (X,Y,TI)
     * @param TI inversion times
     * @param iterations Number of iterations between synthetic and registration
     * @return Motion corrected data
     */
    hoNDArray<vector_td<float,2>> t1_moco_cmr(const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI, unsigned int iterations);

}
