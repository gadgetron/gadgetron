//
// Created by dchansen on 6/16/18.
//
#include <numeric>
#include <boost/math/constants/constants.hpp>
#include "correct_frequency_shift.h"
#include "hoNDFFT.h"
#include "hoNDArray_elemwise.h"
namespace Gadgetron{
    namespace FatWater {

        using namespace std::complex_literals;
        static constexpr float PI = boost::math::constants::pi<float>();

        void correct_frequency_shift(hoNDArray<std::complex<float>> &species_images, const Parameters &parameters){

            uint16_t X = species_images.get_size(0);
            uint16_t Y = species_images.get_size(1);
            uint16_t Z = species_images.get_size(2);
            uint16_t CHA = species_images.get_size(3);
            uint16_t N = species_images.get_size(4);
            uint16_t S = species_images.get_size(5);
            uint16_t LOC = species_images.get_size(6);

            if (LOC > 1) throw std::invalid_argument("Phase shift does not support LOC dimension greater than 1");


            std::vector<size_t> sub_dimension = {X,Y,Z,CHA,N,1,1};
            auto data_ptr = species_images.get_data_ptr();


            for (int kspecies = 0; kspecies < parameters.species.size(); kspecies++){

                auto& species = parameters.species[kspecies];

                auto mean_frequency_offset = std::accumulate(
                        species.amplitude_frequency_pairs.begin(), species.amplitude_frequency_pairs.end(),0.0f,
                                              [](auto v, auto pair){ return v+pair.first.real()*pair.second;}
                                              ) /
                              std::accumulate(
                        species.amplitude_frequency_pairs.begin(),species.amplitude_frequency_pairs.end(), 0.0f,
                                            [](auto v, auto pair) { return v + pair.first.real();});

                if (std::abs(mean_frequency_offset) < 1.0 ) continue;



                hoNDArray<std::complex<float>> phase_ramp(species_images.get_size(0));

                for (int i = 0; i < phase_ramp.get_number_of_elements(); i++){
                    phase_ramp[i] = std::exp(2if*PI*parameters.sample_time_us*1e-6f*float(i)*mean_frequency_offset);
                }


                auto view = hoNDArray<std::complex<float>>(sub_dimension,data_ptr);
                hoNDFFT<float>::instance()->fft(&view,0);

                view *= phase_ramp;

                hoNDFFT<float>::instance()->ifft(&view,0);


                data_ptr += view.get_number_of_elements();




            }









        }
    }
}

