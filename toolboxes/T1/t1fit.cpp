#include "HybridLM.h"
#include <vector>
#include "hoArmadillo.h"

namespace Gadgetron {
    namespace {

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

        template <class T>
        std::tuple<T, T> fit_T1_2param_single(const std::vector<T>& TI, const std::vector<T>& data) {

            T A  = *std::max_element(data.begin(), data.end()) - *std::min_element(data.begin(), data.end());
            T T1 = 800;

            T1Residual_2param<T> f{ TI, data };

            Solver::HybridLMSolver<T> solver(data.size(), 2);
            arma::Col<T> params{ A, T1 };
            solver.solve(f, params);
            return { params[0], params[1] };
        }
    }

    hoNDArray<float> fit_T1_2param_arma(const hoNDArray<float>& data, std::vector<float>& TI) {

        if (data.get_size(2) != TI.size()) {
            throw std::runtime_error("Data and TI do not match");
        }

        if (raw_data.shape(0) != NRESIDUALS)
            throw std::runtime_error("Just no");
        std::array<float, NRESIDUALS> TI_view;
        for (int i = 0; i < (int)TI_view.size(); i++) {
            TI_view[i] = raw_TI[i];
        }

        auto Apy  = py::array_t<float>({ raw_data.shape(1), raw_data.shape(2) });
        auto T1py = py::array_t<float>({ raw_data.shape(1), raw_data.shape(2) });

        auto A  = Apy.mutable_unchecked<2>();
        auto T1 = T1py.mutable_unchecked<2>();

#pragma omp parallel
        {
            std::array<float, NRESIDUALS> data_view;
#pragma omp for
            for (int x = 0; x < (int)raw_data.shape(1); x++) {
                for (int y = 0; y < (int)raw_data.shape(2); y++) {
                    for (int t = 0; t < (int)raw_data.shape(0); t++) {
                        data_view[t] = raw_data(t, x, y);
                    }

                    auto result = fit_T1_2param_single_arma<float, 8>(TI_view, data_view);

                    A(x, y)  = std::get<1>(result);
                    T1(x, y) = std::get<0>(result);
                }
            }
        }
        return { T1py, Apy };
    };
}
