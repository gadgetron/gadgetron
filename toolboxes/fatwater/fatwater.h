#ifndef FATWATER_H
#define FATWATER_H

#include <vector>
#include <utility>

#include "fatwater_export.h"
#include "hoNDArray.h"

namespace Gadgetron
{
    namespace FatWater {

        struct EXPORTFATWATER Config {
            std::pair<float, float> frequency_range = {-500, 500};
            size_t number_of_frequency_samples = 200;

            std::pair<float, float> r2_range = {5, 500};
            size_t number_of_r2_samples = 5;
            size_t number_of_r2_fine_samples = 200;

            size_t number_of_iterations = 40;

            float lambda = 0.02;
            float lambda_extra = 0.01;
            bool do_gradient_descent = true;
            unsigned int downsamples = 0;


        };

        /**
           Amplitudes and frequences (in ppm)
         */
        struct EXPORTFATWATER ChemicalSpecies {
            std::string name;
            std::vector<std::pair<std::complex<float>, float> > amplitude_frequency_pairs;
        };


        struct EXPORTFATWATER Output {
            hoNDArray<std::complex<float>> images;
            hoNDArray<float> field_map;
            hoNDArray<float> r2star_map;
        };

        struct EXPORTFATWATER Parameters {

            float field_strength_T;
            bool precession_is_clockwise;
            float gyromagnetic_ratio_Mhz = 42.57747892;
            float sample_time_us;

            std::vector<float> echo_times_s;
            std::vector<ChemicalSpecies> species;

        };




        /**
           Main interface for water fat separation.

           data array is assumed to be a 7D array [X, Y, Z, CHA, N, S, LOC]

         */
        EXPORTFATWATER Output
        fatwater_separation(const hoNDArray <std::complex<float>> &data, Parameters p, Config config);
    }
}

#endif //FATWATER_H
