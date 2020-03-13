#include "t1fit.h"
#include "HybridLM.h"
#include "hoArmadillo.h"
#include <vector>

#include "demons_registration.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::T1;

    template <class T> struct T1Residual_2param {
        const std::vector<T>& TI;
        const std::vector<T>& measurement;

        void operator()(const arma::Col<T>& params, arma::Col<T>& residual, arma::Mat<T>& jacobian) const {
            const auto& T1 = params[0];
            const auto& A  = params[1];

            for (int i = 0; i < residual.n_elem; i++) {
                residual(i) = T(2) * std::exp(-TI[i] / T1);
            }

            for (int i = 0; i < residual.n_elem; i++) {
                jacobian(i, 0) = A * TI[i] * residual[i] / (T1 * T1);
                jacobian(i, 1) = residual[i] - T(1);
            }

            for (int i = 0; i < residual.n_elem; i++) {
                residual(i) = measurement[i] - A * (T(1) - residual[i]);
            }
        }
    };

    template <class T> struct T1Residual_3param {
        const std::vector<T>& TI;
        const std::vector<T>& measurement;

        void operator()(const arma::Col<T>& params, arma::Col<T>& residual, arma::Mat<T>& jacobian) const {
            const auto& T1 = params[0];
            const auto& A  = params[1];
            const auto& B  = params[2];

            for (int i = 0; i < residual.n_elem; i++) {
                residual(i) = std::exp(-TI[i] / T1);
            }

            for (int i = 0; i < residual.n_elem; i++) {
                jacobian(i, 0) = B * TI[i] * residual[i] / (T1 * T1);
                jacobian(i, 1) = -T(1);
                jacobian(i, 2) = residual[i];
            }

            for (int i = 0; i < residual.n_elem; i++) {
                residual(i) = measurement[i] - (A - B * residual[i]);
            }
        }
    };

    template <class T> std::tuple<T, T> fit_T1_2param_single(const std::vector<T>& TI, const std::vector<T>& data) {

        T A  = *std::max_element(data.begin(), data.end()) - *std::min_element(data.begin(), data.end());
        T T1 = 800;

        T1Residual_2param<T> f{ TI, data };

        Solver::HybridLMSolver<T> solver(data.size(), 2);
        arma::Col<T> params{ T1, A };
        solver.solve(f, params);
        return { params[0], params[1] };
    }

    template <class T> std::tuple<T, T, T> fit_T1_3param_single(const std::vector<T>& TI, const std::vector<T>& data) {

        T A  = *std::max_element(data.begin(), data.end()) - *std::min_element(data.begin(), data.end());
        T B  = *std::max_element(data.begin(), data.end()) - *std::min_element(data.begin(), data.end());
        T T1 = 800;

        T1Residual_3param<T> f{ TI, data };

        Solver::HybridLMSolver<T> solver(data.size(), 3);
        arma::Col<T> params{ T1, A, B };
        solver.solve(f, params);
        return { params[0], params[1], params[2] };
    }

    hoNDArray<float> predict_signal(const T1_2param& params, const std::vector<float>& TI) {

        const auto& [A, T1] = params;

        auto dimensions = A.dimensions();
        dimensions.push_back(TI.size());
        auto result = hoNDArray<float>(dimensions);
        for (int cha = 0; cha < (int)TI.size(); cha++) {
            for (int y = 0; y < (int)A.get_size(1); y++) {
                for (int x = 0; x < (int)A.get_size(0); x++) {
                    result(x, y, cha) = A(x, y) * (2 - std::exp(-TI[cha] / T1(x, y)));
                }
            }
        }

        return result;
    }

    hoNDArray<float> predict_signal(const T1_3param& params, const std::vector<float>& TI) {

        const auto& [A, B, T1] = params;

        auto result = hoNDArray<float>(A.dimensions());
        for (int cha = 0; cha < (int)TI.size(); cha++) {
            for (int y = 0; y < (int)A.get_size(1); y++) {
                for (int x = 0; x < (int)A.get_size(0); x++) {
                    result(x, y, cha) = A(x, y) - B(x, y) * std::exp(-TI[cha] / T1(x, y));
                }
            }
        }

        return result;
    }

}

T1_2param Gadgetron::T1::fit_T1_2param(const hoNDArray<float>& data, const std::vector<float>& TI) {

    if (data.get_size(2) != TI.size()) {
        throw std::runtime_error("Data and TI do not match");
    }

    auto A  = hoNDArray<float>({ data.get_size(0), data.get_size(1) });
    auto T1 = A;

#pragma omp parallel
    {
        std::vector<float> data_view(TI.size());
#pragma omp for
        for (int y = 0; y < (int)data.get_size(1); y++) {
            for (int x = 0; x < (int)data.get_size(0); x++) {
                for (int t = 0; t < (int)TI.size(); t++) {
                    data_view[t] = data(x, y, t);
                }

                auto result = fit_T1_2param_single<float>(TI, data_view);

                A(x, y)  = std::get<1>(result);
                T1(x, y) = std::get<0>(result);
            }
        }
    }
    return { A, T1 };
};

T1_3param Gadgetron::T1::fit_T1_3param(const hoNDArray<float>& data, const std::vector<float>& TI) {

    if (data.get_size(2) != TI.size()) {
        throw std::runtime_error("Data and TI do not match");
    }

    auto A  = hoNDArray<float>({ data.get_size(0), data.get_size(1) });
    auto B  = hoNDArray<float>({ data.get_size(0), data.get_size(1) });
    auto T1 = hoNDArray<float>({ data.get_size(0), data.get_size(1) });

#pragma omp parallel
    {
        std::vector<float> data_view(TI.size());
#pragma omp for
        for (int y = 0; y < (int)data.get_size(1); y++) {
            for (int x = 0; x < (int)data.get_size(0); x++) {
                for (int t = 0; t < (int)TI.size(); t++) {
                    data_view[t] = data(x, y, t);
                }

                auto result = fit_T1_3param_single<float>(TI, data_view);

                A(x, y)  = std::get<1>(result);
                B(x, y)  = std::get<2>(result);
                T1(x, y) = std::get<0>(result);
            }
        }
    }
    return { A, B, T1 };
}
namespace {

    hoNDArray<float> phase_correct(const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI) {

        hoNDArray<float> result(data.dimensions());
        auto max_TI_index = std::distance(TI.begin(), std::max_element(TI.begin(), TI.end()));
        for (long long y = 0; y < (long long)data.get_size(1); y++) {
            for (long long x = 0; x < (long long)data.get_size(0); x++) {
                auto reference_phase = std::arg(data(x, y, max_TI_index));

                for (long long i = 0; i < (long long)data.get_size(2); i++) {
                    auto val        = data(x, y, i);
                    result(x, y, i) = std::polar(std::abs(val), std::arg(val) - reference_phase).real();
                }
            }
        }

        return result;
    }

    hoNDArray<std::complex<float>> register_and_deform_groups(const hoNDArray<float>& phase_corrected_data,
        const hoNDArray<float>& predicted, const hoNDArray<std::complex<float>>& data) {
        using namespace Indexing;
        hoNDArray<std::complex<float>> result(predicted.dimensions());

#pragma omp parallel for
        for (long long cha = 0; cha < (long long)data.get_size(2); cha++) {
            auto vector_field = Registration::diffeomorphic_demons<float, 2>(
                predicted(slice, slice, cha), phase_corrected_data(slice, slice, cha), 40, 2.0);
            result(slice, slice, cha)
                = Registration::deform_image<std::complex<float>, 2>(data(slice, slice, cha), vector_field);
        }

        return result;
    }

}

T1_3param Gadgetron::T1::motion_compensated_t1_fit(
    const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI, unsigned int iterations) {
    if (data.get_size(2) != TI.size()) {
        throw std::runtime_error("Data and TI do not match");
    }

    auto corrected = phase_correct(data, TI);

    auto corrected_orig = corrected;

    auto parameters = fit_T1_2param(corrected, TI);

    auto predicted = predict_signal(parameters, TI);

    for (int i = 0; i < iterations; i++) {
        std::cout << "Iteration " << i << std::endl;

        auto deformed_data = register_and_deform_groups(corrected_orig, predicted, data);
        corrected = phase_correct(deformed_data, TI);

        parameters = fit_T1_2param(corrected, TI);

        predicted = predict_signal(parameters, TI);

    }

    auto result = fit_T1_3param(corrected, TI);

    return result;
};
