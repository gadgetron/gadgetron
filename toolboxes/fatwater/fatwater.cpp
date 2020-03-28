#include "fatwater.h"

#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "hoArmadillo.h"

#include <boost/config.hpp>

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>

#include <boost/timer/timer.hpp>
#include <boost/iterator/function_input_iterator.hpp>
#include <iterator>
#include "graph_cut.h"

#include "hoNDArray_math.h"

#include <boost/random.hpp>
#include <boost/math/constants/constants.hpp>
#include <armadillo>
#include <random>
#include <cpu/hoNDArray_fileio.h>
#include <GadgetronTimer.h>
#include <complex>
#include <cpu/math/hoNDImage_util.h>

#ifdef max
#undef max
#endif // max
#ifdef min
#undef min
#endif // min

#include <algorithm>
#include "bounded_field_map.h"

using namespace boost;


namespace Gadgetron {
    namespace FatWater {
        using namespace std::complex_literals;
        static constexpr float PI = boost::math::constants::pi<float>();
        static std::mt19937 rng_state(4242);


        hoNDArray<std::vector<uint16_t>> find_local_minima(const hoNDArray<float> &residuals, float threshold = 0.00f) {


            auto threshold_signal = sum(residuals, 0);
            threshold_signal /= max(threshold_signal);
            sqrt_inplace(&threshold_signal);

            auto min_residuals = min(residuals, 0);
            auto max_residuals = max(residuals, 0);


            const auto Y = residuals.get_size(2);
            const auto X = residuals.get_size(1);
            hoNDArray<std::vector<uint16_t>> result(X, Y);
            const auto steps = residuals.get_size(0);
            size_t count = 0;
            for (size_t k2 = 0; k2 < Y; k2++) {
                for (size_t k1 = 0; k1 < X; k1++) {

                    std::vector<uint16_t> minima;
//                    if (threshold_signal(k1, k2) > threshold) {
                    for (size_t k0 = 1; k0 < steps - 1; k0++) {
                        if ((residuals(k0, k1, k2) < residuals(k0 - 1, k1, k2)) &&
                            (residuals(k0 + 1, k1, k2) >= residuals(k0, k1, k2))) {
                            minima.push_back(k0);
                        }

                    }
//                    }
                    if (minima.size() < 2) count++;
                    result(k1, k2) = std::move(minima);
                }
            }

            GDEBUG("COUNT DRACULA %i\n", count);
            return result;
        }

        hoNDArray<float> approx_second_derivative(const hoNDArray<float> &residuals,
                                                  const hoNDArray<std::vector<uint16_t>> &local_min_indices,
                                                  float step_size) {
            hoNDArray<float> second_deriv(local_min_indices.dimensions());

            const auto Y = second_deriv.get_size(1);
            const auto X = second_deriv.get_size(0);
            const auto nfields = residuals.get_size(0);

            for (uint16_t k2 = 0; k2 < Y; k2++) {
                for (uint16_t k1 = 0; k1 < X; k1++) {

                    auto minimum = std::min_element(&residuals(1, k1, k2), &residuals(nfields - 1, k1, k2)) -
                                   &residuals(0, k1, k2);

                    auto sd =
                            (residuals(minimum - 1, k1, k2) + residuals(minimum + 1, k1, k2) -
                             2 * residuals(minimum, k1, k2)) / (step_size * step_size);

                    second_deriv(k1, k2) = sd;

                }
            }

            return second_deriv;


        }


        hoNDArray<float>
        create_field_map(const hoNDArray<uint16_t> &field_map_index, const std::vector<float> &field_map_strengths) {
            const uint16_t max_val = field_map_strengths.size() - 1;
            hoNDArray<float> field_map(field_map_index.dimensions());
            std::transform(field_map_index.begin(), field_map_index.end(), field_map.begin(),
                           [&](uint16_t i) { return field_map_strengths[std::min(i, max_val)]; });
            return field_map;
        }


        hoNDArray<uint16_t> create_field_map_proposal1(const hoNDArray<uint16_t> &field_map_index,
                                                       const hoNDArray<std::vector<uint16_t>> &minima,
                                                       const hoNDArray<float> &residuals,
                                                       const std::vector<float> &field_map_strengths, float fat_freq,
                                                       float dF, float dTE) {

            const size_t elements = field_map_index.get_number_of_elements();
            hoNDArray<uint16_t> proposed_field_map_index(field_map_index.dimensions());
            const size_t field_maps = field_map_strengths.size();
            std::uniform_int_distribution<int> coinflip(0, 1);
            int jump;
            if (coinflip(rng_state)) {
                jump = round(std::abs(fat_freq / dF));

            } else {
                jump = round((1.0 / dTE - std::abs(fat_freq)) / dF);
            }


            for (size_t i = 0; i < elements; i++) {
                auto &mins = minima[i];
                auto fbi = field_map_index[i];
                auto fqmi = std::find_if(mins.begin(), mins.end(),
                                         [&](auto fqi) { return fqi > fbi + 20; }); //Find smallest
                proposed_field_map_index[i] = (fqmi == mins.end()) ? std::min<int>(fbi + jump, field_maps - 1) : *fqmi;
            }

            return proposed_field_map_index;


        }

        hoNDArray<uint16_t> create_field_map_proposal2(const hoNDArray<uint16_t> &field_map_index,
                                                       const hoNDArray<std::vector<uint16_t>> &minima,
                                                       const hoNDArray<float> &residuals,
                                                       const std::vector<float> &field_map_strengths, float fat_freq,
                                                       float dF, float dTE) {

            const size_t elements = field_map_index.get_number_of_elements();
            hoNDArray<uint16_t> proposed_field_map_index(field_map_index.dimensions());
            std::uniform_int_distribution<int> coinflip(0, 1);
            int jump;
            if (coinflip(rng_state)) {
                jump = round(std::abs(fat_freq / dF));
            } else {
                jump = round((1.0 / dTE - std::abs(fat_freq)) / dF);
            }

            for (size_t i = 0; i < elements; i++) {
                auto &mins = minima[i];
                int fbi = field_map_index[i];
                auto fqmi = std::find_if(mins.rbegin(), mins.rend(),
                                         [&](auto fqi) { return fqi < (fbi - 20); }); //Find smallest
                proposed_field_map_index[i] = (fqmi == mins.rend()) ? std::max<int>(fbi - jump, 0) : *fqmi;
            }

            return proposed_field_map_index;


        }

        hoNDArray<uint16_t>
        create_field_map_proposal_standard(const hoNDArray<uint16_t> &field_map_index, int sign,
                                           uint16_t max_field_value) {


            std::uniform_int_distribution<int> rng(1, 3);

            int step_size = sign * rng(rng_state);
//        int step_size = sign;
            hoNDArray<uint16_t> proposed_field_map_index(field_map_index.dimensions());
            std::transform(field_map_index.begin(), field_map_index.end(), proposed_field_map_index.begin(),
                           [&](uint16_t j) {
                               return uint16_t(std::min(std::max(j + step_size, 0), int(max_field_value)));
                           });

            return proposed_field_map_index;
        }

        hoNDArray<uint16_t>
        solve_MRF(const Config &config, const std::vector<float> &field_map_strengths,
                  const hoNDArray<float> &residual, const hoNDArray<std::vector<uint16_t>> &local_min_indices,
                  const hoNDArray<float> &second_deriv, float fat_freq, float dF, float dTE) {

            hoNDArray<uint16_t> fmIndex(local_min_indices.dimensions());

            std::uniform_int_distribution<int> coinflip(0, 2);
            fmIndex.fill(field_map_strengths.size() / 2);

            hoNDArray<uint16_t> fmIndex_update;
            for (int i = 0; i < config.number_of_iterations; i++) {
                if (coinflip(rng_state) == 0 || i < 15) {
                    if (!(i % 2)) {
                        fmIndex_update = create_field_map_proposal1(fmIndex, local_min_indices, residual,
                                                                    field_map_strengths, fat_freq, dF, dTE);
                    } else {
                        fmIndex_update = create_field_map_proposal2(fmIndex, local_min_indices, residual,
                                                                    field_map_strengths, fat_freq, dF, dTE);
                    }
                } else {
                    fmIndex_update = create_field_map_proposal_standard(fmIndex, std::pow(-1, i),
                                                                        field_map_strengths.size() - 1);
                }

                fmIndex = update_field_map(fmIndex, fmIndex_update, residual, second_deriv);
            }

            return fmIndex;
        }

        arma::Mat<std::complex<float>>
        calculate_psi_matrix(const std::vector<float> &echoTimes,
                             const arma::Mat<std::complex<float>> &phiMatrix,  float r2star) {
            arma::Mat<std::complex<float>> psiMatrix(phiMatrix.n_rows, phiMatrix.n_cols);
            for (int k1 = 0; k1 < phiMatrix.n_rows; k1++) {
                auto curModulation = exp(-r2star * (echoTimes[k1] - echoTimes[0]));
                for (int k2 = 0; k2 < phiMatrix.n_cols; k2++) {
                    psiMatrix(k1, k2) = phiMatrix(k1, k2) * curModulation;
                }
            }


            return psiMatrix;
        }


        hoNDArray<arma::Mat<std::complex<float>>>
        calculate_projection_matrices(const std::vector<float> &echo_times,
                                      const arma::Mat<std::complex<float>> &phiMatrix,
                                      const std::vector<float> &field_map_strengths,
                                      const std::vector<float> &r2stars) {

            auto num_fm = field_map_strengths.size();
            auto num_r2star = r2stars.size();
            hoNDArray<arma::Mat<std::complex<float>>> Ps(num_fm, num_r2star);
            size_t nte = echo_times.size();

#ifdef WIN32
#pragma omp parallel for
#else
#pragma omp parallel for
#endif
            for (int k3 = 0; k3 < num_fm; k3++) {

                float fm = field_map_strengths[k3];
                arma::Col<std::complex<float>> b_shifts(nte);
                for (int kt = 0; kt < nte; kt++)
                    b_shifts[kt] = std::exp(2if * PI * (echo_times[kt] - echo_times[0]) * fm);
                for (int k4 = 0; k4 < num_r2star; k4++) {
                    float r2star = r2stars[k4];

                    arma::Mat<std::complex<float>> psiMatrix = calculate_psi_matrix(echo_times, phiMatrix,
                                                                                    r2star);
                    Ps(k3, k4) = arma::diagmat(b_shifts) * (arma::eye<arma::Mat<std::complex<float>>>(nte, nte) -
                                                            psiMatrix * arma::pinv(psiMatrix)) *
                                 arma::diagmat(arma::conj(b_shifts));

                }
            }
            return Ps;
        }




        hoNDArray<float>
        calculate_r2star_map(const hoNDArray<std::complex<float> > &data, const hoNDArray<float> &field_map,
                             const std::vector<float> &r2star_values,
                             const arma::Mat<std::complex<float>> &phiMatrix, const std::vector<float> &echoTimes) {
            using cMat = arma::Mat<std::complex<float>>;
            uint16_t X = data.get_size(0);
            uint16_t Y = data.get_size(1);
            uint16_t Z = data.get_size(2);
            uint16_t CHA = data.get_size(3);
            uint16_t N = data.get_size(4);
            uint16_t S = data.get_size(5);
            uint16_t LOC = data.get_size(6);
            auto nte = phiMatrix.n_rows;


            auto data_corrected = data;

            for (int kS = 0; kS < S; kS++) {
                for (int kN = 0; kN < N; kN++) {
                    for (int kcha = 0; kcha < CHA; kcha++) {
                        for (int ky = 0; ky < Y; ky++) {
                            for (int kx = 0; kx < X; kx++) {
                                data_corrected(kx, ky, 0, kcha, kN, kS) *= std::exp(
                                        -2if * PI * field_map(kx, ky) * echoTimes[kS]);
                            }
                        }
                    }
                }
            }


            std::vector<arma::Mat<std::complex<float>>> Ps;
            for (auto r2star : r2star_values) {

                auto psiMatrix = calculate_psi_matrix(echoTimes, phiMatrix,
                                                      r2star);
                Ps.emplace_back(
                        arma::eye<arma::Mat<std::complex<float>>>(nte, nte) - psiMatrix * arma::pinv(psiMatrix)
                );

            }


            hoNDArray<float> r2star_map(field_map.dimensions());
#ifdef WIN32
#pragma omp parallel for
#else
#pragma omp parallel for collapse(2)
#endif
            for (int k2 = 0; k2 < Y; k2++) {
                for (int k1 = 0; k1 < X; k1++) {
                    // Get current signal
                    std::vector<cMat> signals(CHA, cMat(S, N));
                    for (int cha = 0; cha < CHA; cha++) {
                        auto &tempSignal = signals[cha];
                        for (int k4 = 0; k4 < N; k4++) {
                            for (int k5 = 0; k5 < S; k5++) {
                                tempSignal(k5, k4) = data_corrected(k1, k2, 0, cha, k4, k5, 0);

                            }
                        }
                    }


                    float minResidual = std::numeric_limits<float>::max();
                    std::vector<float> field_map_strengths = {field_map(k1, k2)};


                    for (int kr2 = 0; kr2 < r2star_values.size(); kr2++) {

                        float curResidual = 0;

                        for (int cha = 0; cha < CHA; cha++) {
                            // Apply projector
                            arma::Mat<std::complex<float>> projected = Ps[kr2] * signals[cha];
                            curResidual += std::accumulate(projected.begin(), projected.end(), 0.0f,
                                                           [](auto v1, auto v2) {
                                                               return v1 +
                                                                      std::norm(v2);
                                                           });
                        }
                        if (curResidual < minResidual) {
                            minResidual = curResidual;
                            r2star_map(k1, k2) = r2star_values[kr2];
                        }

                    }


                }
            }

            return r2star_map;
        }


        std::tuple<hoNDArray<float>, hoNDArray<uint16_t>>
        calculate_residual_and_r2star(const hoNDArray<std::complex<float>> &data, const Parameters &parameters,
                                      const arma::Mat<std::complex<float>> &phi,
                                      const std::vector<float> &field_strengths,
                                      const std::vector<float> &r2star_values) {
            using cMat = arma::Mat<std::complex<float>>;
            uint16_t X = data.get_size(0);
            uint16_t Y = data.get_size(1);
            uint16_t Z = data.get_size(2);
            uint16_t CHA = data.get_size(3);
            uint16_t N = data.get_size(4);
            uint16_t S = data.get_size(5);
            uint16_t LOC = data.get_size(6);


            auto projection_matrices = calculate_projection_matrices(parameters.echo_times_s, phi, field_strengths,
                                                                     r2star_values);

            auto result = std::make_tuple(hoNDArray<float>(field_strengths.size(), X, Y, Z),
                                          hoNDArray<uint16_t>(X, Y, Z, field_strengths.size()));

            auto &residual = std::get<0>(result);
            auto &r2starIndex = std::get<1>(result);

#ifdef WIN32
#pragma omp parallel for 
#else
#pragma omp parallel for collapse(3)
#endif
            for (int kz = 0; kz < Z; kz++) {
                for (int ky = 0; ky < Y; ky++) {
                    for (int kx = 0; kx < X; kx++) {

                        std::vector<cMat> signals(CHA, cMat(S, N));
                        for (int cha = 0; cha < CHA; cha++) {
                            auto &tempSignal = signals[cha];
                            for (int kn = 0; kn < N; kn++) {
                                for (int ks = 0; ks < S; ks++) {
                                    tempSignal(ks, kn) = data(kx, ky, kz, cha, kn, ks, 0);
                                }
                            }
                        }

                        for (int kf = 0; kf < field_strengths.size(); kf++) {

                            float minResidual = std::numeric_limits<float>::max();

                            for (int kr = 0; kr < r2star_values.size(); kr++) {
                                // Apply projector
                                float curResidual = 0;
                                for (int cha = 0; cha < CHA; cha++) {
                                    arma::Mat<std::complex<float>> projected =
                                            projection_matrices(kf, kr) * signals[cha];
                                    curResidual += std::accumulate(projected.begin(), projected.end(), 0.0f,
                                                                   [](auto v1, auto v2) {
                                                                       return v1 +
                                                                              std::norm(v2);
                                                                   });
                                }
                                if (curResidual < minResidual) {
                                    minResidual = curResidual;
                                    r2starIndex(kx, ky, kz, kf) = kr;
                                }
                            }
                            residual(kf, kx, ky, kz) = minResidual;

                        }
                    }
                }
            }

            return result;
        }

        hoNDArray<std::complex<float>>
        separate_species(const hoNDArray<std::complex<float>> &data,
                         const arma::Mat<std::complex<float>> &phiMatrix, hoNDArray<float> &r2star_map,
                         hoNDArray<float> &field_map, const Parameters &parameters) {

            using cMat = arma::Mat<std::complex<float>>;
            uint16_t X = data.get_size(0);
            uint16_t Y = data.get_size(1);
            uint16_t Z = data.get_size(2);
            uint16_t CHA = data.get_size(3);
            uint16_t N = data.get_size(4);
            uint16_t S = data.get_size(5);
            uint16_t LOC = data.get_size(6);
            hoNDArray<std::complex<float> > out(X, Y, Z, CHA, N, parameters.species.size(),
                                                LOC); // S dimension gets replaced by water/fat stuff

#ifdef WIN32
#pragma omp parallel for 
#else
#pragma omp parallel for collapse(3)
#endif
            for (int kz = 0; kz < Z; kz++) {
                for (int ky = 0; ky < Y; ky++) {
                    for (int kx = 0; kx < X; kx++) {
                        std::vector<cMat> signals(CHA, cMat(S, N));

                        // Get current signal
                        for (int cha = 0; cha < CHA; cha++) {
                            auto &tempSignal = signals[cha];
                            for (int kn = 0; kn < N; kn++) {
                                for (int ks = 0; ks < S; ks++) {
                                    tempSignal(ks, kn) = data(kx, ky, kz, cha, kn, ks, 0);
                                }
                            }
                        }

                        auto fm = field_map(kx, ky);
                        auto r2star = r2star_map(kx, ky);

                        arma::Mat<std::complex<float>> psiMatrix = calculate_psi_matrix(parameters.echo_times_s,
                                                                                        phiMatrix,r2star);

                        auto nte = parameters.echo_times_s.size();
                        arma::Col<std::complex<float>> b_shifts(nte);
                        for (int kt = 0; kt < nte; kt++)
                            b_shifts[kt] = std::exp(
                                    2if * PI * (parameters.echo_times_s[kt] - parameters.echo_times_s[0]) * fm);
                        psiMatrix = arma::diagmat(b_shifts) * psiMatrix;
                        // Solve for water and fat

                        for (int cha = 0; cha < CHA; cha++) {
                            arma::Mat<std::complex<float>> curWaterFat = arma::solve(psiMatrix, signals[cha]);
                            for (int kn = 0; kn < N; kn++) {
                                for (int kspecies = 0; kspecies <
                                                       parameters.species.size(); kspecies++) { // 2 elements for water and fat currently
                                    out(kx, ky, kz, cha, kn, kspecies, 0) = curWaterFat(kspecies, kn);
                                }
                            }
                        }

                    }
                }
            }
            return out;
        }


        arma::Mat<std::complex<float>>
        calculatePhiMatrix(const Parameters parameters) {

            auto echoTimes = parameters.echo_times_s;

            auto nte = echoTimes.size();
            auto nspecies = parameters.species.size();
            typedef arma::Mat<std::complex<float>> Cmat;
            Cmat phiMatrix = arma::zeros<Cmat>(nte, nspecies);
            for (int k1 = 0; k1 < nte; k1++) {
                for (int k2 = 0; k2 < nspecies; k2++) {
                    auto &species = parameters.species[k2];
                    auto npeaks = species.amplitude_frequency_pairs.size();
                    for (int k3 = 0; k3 < npeaks; k3++) {
                        auto relAmp = species.amplitude_frequency_pairs[k3].first;
                        auto freq_hz = parameters.field_strength_T * parameters.gyromagnetic_ratio_Mhz *
                                       species.amplitude_frequency_pairs[k3].second;
                        phiMatrix(k1, k2) += relAmp * exp(2if * PI * echoTimes[k1] * freq_hz);

                    }

                }
            }
            return phiMatrix;
        }


        std::vector<float> linspace(const std::pair<float, float> &range, unsigned int count) {
            std::vector<float> result(count);
            result[0] = range.first;
            for (int i = 1; i < count; i++) {
                result[i] = range.first + i * (range.second - range.first) / (count - 1);
            }
            return result;
        }

        std::complex<float> mean_frequency(const Parameters &parameters) {
            ChemicalSpecies fat = parameters.species[1];
            auto average_fat_freq =
                    accumulate(fat.amplitude_frequency_pairs.begin(), fat.amplitude_frequency_pairs.end(), std::complex<float>(0.0),
                               [](auto val, auto tup) {
                                   return val + std::get<0>(tup) * std::get<1>(tup);
                               }) /
                    accumulate(fat.amplitude_frequency_pairs.begin(), fat.amplitude_frequency_pairs.end(), std::complex<float>(0.0),
                               [](auto val, auto tup) { return val + std::get<0>(tup); });
            average_fat_freq *= parameters.field_strength_T * parameters.gyromagnetic_ratio_Mhz;
            return average_fat_freq;
        }

        std::tuple<hoNDArray<float>, hoNDArray<float>>
        calculate_field_map(const hoNDArray<std::complex<float>> &data, const Parameters &parameters,
                            const Config &config,
                            const std::complex<float> &average_fat_freq, const arma::Mat<std::complex<float>> &phi,
                            const std::vector<float> &field_map_strengths);

        FatWater::Output
        fatwater_separation(const hoNDArray<std::complex<float> > &data, Parameters parameters,
                            Config config) {

            uint16_t X = data.get_size(0);
            uint16_t Y = data.get_size(1);
            uint16_t Z = data.get_size(2);
            uint16_t CHA = data.get_size(3);
            uint16_t N = data.get_size(4);
            uint16_t S = data.get_size(5);
            uint16_t LOC = data.get_size(6);


            auto echo_times = parameters.echo_times_s;
            GDEBUG("In toolbox - Field Strength: %f T \n", parameters.field_strength_T);
            for (auto &te: echo_times) {
                GDEBUG("In toolbox - Echo time: %f seconds \n", te);
            }
            GDEBUG("In toolbox - PrecessionIsClockwise: %d \n", parameters.precession_is_clockwise);

            //Get or set some algorithm parameters

            std::complex<float> average_fat_freq = mean_frequency(parameters);


            arma::Mat<std::complex<float>> phi = calculatePhiMatrix(parameters);

            std::vector<float> field_map_strengths = linspace(config.frequency_range,
                                                              config.number_of_frequency_samples);


            hoNDArray<float> field_map, r2star_map;
            std::tie(field_map, r2star_map) = calculate_field_map(data, parameters, config, average_fat_freq, phi,
                                                                  field_map_strengths);


            auto species = separate_species(data, phi, r2star_map, field_map, parameters);

            if (config.do_gradient_descent) {
                bounded_field_map(field_map, data, parameters, (field_map_strengths[1] - field_map_strengths[0])*2);
                species = separate_species(data, phi, r2star_map, field_map, parameters);
            }

            return Output{std::move(species), std::move(field_map), std::move(r2star_map)};
        }

        std::tuple<hoNDArray<float>, hoNDArray<float>>
        calculate_field_map(const hoNDArray<std::complex<float>> &data, const Parameters &parameters,
                            const Config &config,
                            const std::complex<float> &average_fat_freq, const arma::Mat<std::complex<float>> &phi,
                            const std::vector<float> &field_map_strengths) {


            std::vector<float> r2stars_fine = linspace(config.r2_range, config.number_of_r2_fine_samples);
            std::vector<float> r2stars = linspace(config.r2_range, config.number_of_r2_samples);

            auto data_scaled = data;
            for (int downsamples = 0; downsamples < config.downsamples; downsamples++)
                data_scaled = downsample<std::complex<float>, 2>(data_scaled);

            hoNDArray<float> residual;
            hoNDArray<uint16_t> r2starIndex;
            std::tie(residual, r2starIndex) = calculate_residual_and_r2star(data_scaled, parameters, phi,
                                                                            field_map_strengths,
                                                                            r2stars);

            auto dF = field_map_strengths[1] - field_map_strengths[0];

            hoNDArray<std::vector<uint16_t>> local_min_indices = find_local_minima(residual);

            hoNDArray<float> lambda_map = approx_second_derivative(residual, local_min_indices, dF);
            lambda_map += mean(&lambda_map) * config.lambda_extra;
            lambda_map *= config.lambda * dF * dF;

            hoNDArray<uint16_t> fmIndex = solve_MRF(config, field_map_strengths, residual, local_min_indices,
                                                    lambda_map, abs(average_fat_freq), dF,
                                                    parameters.echo_times_s[1] - parameters.echo_times_s[0]);


            hoNDArray<float> field_map = create_field_map(fmIndex, field_map_strengths);

            if (config.downsamples) {
                field_map = upsample_spline<float, 2>(field_map,std::pow(2,config.downsamples));
            }
//            fmIndex = upsample<uint16_t,2>(&fmIndex);

            auto r2star_map = calculate_r2star_map(data, field_map, r2stars_fine, phi,
                                                   parameters.echo_times_s);


//            r2star_map = upsample<float,2>(&r2star_map);
            return std::make_tuple(std::move(field_map), std::move(r2star_map));
        }


    }
}
