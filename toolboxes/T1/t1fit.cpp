#include "t1fit.h"
#include "HybridLM.h"
#include "hoArmadillo.h"
#include <vector>

#include "demons_registration.h"

namespace {
    using namespace Gadgetron;

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
}

T1_2param Gadgetron::fit_T1_2param(const hoNDArray<float>& data, std::vector<float>& TI) {

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

T1_3param Gadgetron::fit_T1_3param(const hoNDArray<float>& data, std::vector<float>& TI) {

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

                auto result = fit_T1_3param_single<float>(TI, data_view);

                A(x, y)  = std::get<1>(result);
                T1(x, y) = std::get<0>(result);
            }
        }
    }
    return { A, T1 };
}
namespace {

    hoNDArray<float> phase_correct(const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI){

        hoNDArray<float> result(data.dimensions());
        auto max_TI_index = std::distance(TI.begin(),std::max_element(TI.begin(),TI.end()));
        for (long long y = 0; y < (long long)data.get_size(1); y++){
            for (long long x = 0; x < (long long)data.get_size(0); x++){
                auto reference_phase = std::arg(data(x,y,max_TI_index));

                for (long long i = 0; i < (long long)data.get_size(2); i++ ){
                    auto val = data(x,y,i);
                    result(x,y,i) = std::polar(std::abs(val),std::arg(val)-reference_phase).real();
                }

            }
        }
    }


}

hoNDArray<float> motion_compensated_t1_fit(const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI) {
    if (data.get_size(2) != TI.size()) {
        throw std::runtime_error("Data and TI do not match");
    }



};
