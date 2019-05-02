//
// Created by dchansen on 9/4/18.
//

#include "bounded_field_map.h"
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/minima.hpp>
#include <numeric>
#include <GadgetronTimer.h>

constexpr float PI = boost::math::constants::pi<float>();


namespace Gadgetron {
    namespace FatWater {
        namespace {
            template<unsigned int N>
            struct FieldMapModel {

                FieldMapModel(const Parameters &parameters,
                              const std::array<complext<float>, N> &data) : TEs_(parameters.echo_times_s) {

                    std::transform(data.begin(), data.end(), angles.begin(), [](auto c) { return arg(c); });


                    auto data_norm = std::accumulate(data.begin(), data.end(), 0.0,
                                                     [](auto acc, auto c) { return acc + norm(c); });
                    for (int i = 0; i < data.size(); i++) {
                        for (int j = 0; j < data.size(); j++) {
                            weights[j + i * N] = norm(data[i] * data[j]) / data_norm;
                        }
                    }
                }


                float operator()(float field_value) const {

                    float result = 0;
                    for (int i = 0; i < N; i++) {
                        for (int j = 0; j < N; j++) {
                            result += magnitude_internal(field_value, TEs_[i], TEs_[j], angles[i], angles[j],
                                                         weights[j + i * N]);
                        }
                    }
                    return result;
                }

                float magnitude_internal(float field_value, float time1, float time2, float angle1, float angle2,
                                         float weight) const {
                    assert(weight >= 0);
                    return weight * (1.0f - std::cos(field_value * (time1 - time2) + angle1 - angle2));

                }

                const std::vector<float> TEs_;
                std::array<float, N * N> weights;
                std::array<float, N> angles;

            };


            template<unsigned int N>
            void bounded_field_map_N(Gadgetron::hoNDArray<float> &field_map,
                                     const Gadgetron::hoNDArray<std::complex<float>> &input_data,
                                     const Gadgetron::FatWater::Parameters &parameters,
                                     float delta_field) {


                const size_t X = input_data.get_size(0);
                const size_t Y = input_data.get_size(1);
                const size_t Z = input_data.get_size(2);
                const size_t S = input_data.get_size(5);

#ifdef WIN32
    #pragma omp parallel for
#else
    #pragma omp parallel for collapse(2)
#endif
                for (int ky = 0; ky < Y; ky++) {
                    for (size_t kx = 0; kx < X; kx++) {

                        std::array<complext<float>, N> signal;

                        for (int k3 = 0; k3 < S; k3++) {
                            signal[k3] = input_data(kx, ky, 0, 0, 0, k3, 0);
                        }


                        auto model = FieldMapModel<N>(parameters, signal);

                        auto result_pair = boost::math::tools::brent_find_minima(model, field_map(kx,ky)-delta_field,field_map(kx,ky)+delta_field, 24);
                        field_map(kx, ky) = result_pair.first;
                    }
                }



            }
        }

        void bounded_field_map(Gadgetron::hoNDArray<float> &field_map,
                                 const Gadgetron::hoNDArray<std::complex<float>> &input_data,
                                 const Gadgetron::FatWater::Parameters &parameters,
                                 float delta_field
        ) {


            if (input_data.get_size(4) > 1) throw std::runtime_error("Only single repetition supported");

            switch (input_data.get_size(5)) {
                case 2:
                    bounded_field_map_N<2>(field_map, input_data, parameters, delta_field);
                    break;
                case 3:
                    bounded_field_map_N<3>(field_map, input_data, parameters, delta_field);
                    break;
                case 4:
                    bounded_field_map_N<4>(field_map, input_data, parameters, delta_field);
                    break;
                case 5:
                    bounded_field_map_N<5>(field_map, input_data, parameters, delta_field);
                    break;
                case 6:
                    bounded_field_map_N<6>(field_map, input_data, parameters, delta_field);
                    break;
                case 7:
                    bounded_field_map_N<7>(field_map, input_data, parameters, delta_field);
                    break;
                case 8:
                    bounded_field_map_N<8>(field_map, input_data, parameters, delta_field);
                    break;
                default:
                    throw std::runtime_error("Unsupported number of echoes");

            }
        }

    }
}

