#include "t1fit.h"
#include "HybridLM.h"
#include "hoArmadillo.h"
#include <vector>
#include "hoNDArray_math.h"

#include "demons_registration.h"
#include "hoNDArray_fileio.h"

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
        auto status = solver.solve(f, params);
        switch (status) {
            case Solver::ReturnStatus::LINEAR_SOLVER_FAILED:
                return {std::numeric_limits<T>::quiet_NaN(),std::numeric_limits<T>::quiet_NaN()};
            case Solver::ReturnStatus::MAX_ITERATIONS_REACHED:
                return {0,0};
            case Solver::ReturnStatus::SUCCESS: break;
            }

        return { params[0], params[1] };
    }

    template <class T> std::tuple<T, T, T> fit_T1_3param_single(const std::vector<T>& TI, const std::vector<T>& data) {

        T A  = *std::max_element(data.begin(), data.end());
        T B  = A - *std::min_element(data.begin(), data.end());
        T T1 = 800;

        T1Residual_3param<T> f{ TI, data };

        Solver::HybridLMSolver<T> solver(data.size(), 3);
        arma::Col<T> params{ T1, A, B };
        auto status = solver.solve(f, params);
        switch (status) {
        case Solver::ReturnStatus::LINEAR_SOLVER_FAILED:
            return {std::numeric_limits<T>::quiet_NaN(),std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};
        case Solver::ReturnStatus::MAX_ITERATIONS_REACHED:
            return {0,0,0};
        case Solver::ReturnStatus::SUCCESS: break;
        }
        return { params[0], params[1], params[2] };
    }



}
    hoNDArray<float> Gadgetron::T1::predict_signal(const T1_2param& params, const std::vector<float>& TI) {

        const auto& [A, T1] = params;

        auto dimensions = A.dimensions();
        dimensions.push_back(TI.size());
        auto result = hoNDArray<float>(dimensions);
        for (int cha = 0; cha < (int)TI.size(); cha++) {
            for (int y = 0; y < (int)A.get_size(1); y++) {
                for (int x = 0; x < (int)A.get_size(0); x++) {
                    result(x, y, cha) = A(x, y) * (1 - 2*std::exp(-TI[cha] / T1(x, y)));
                }
            }
        }

        return result;
    }

    hoNDArray<float> Gadgetron::T1::predict_signal(const T1_3param& params, const std::vector<float>& TI) {

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
hoNDArray<float> Gadgetron::T1::phase_correct(const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI) {

        hoNDArray<float> result(data.dimensions());
        auto max_TI_index = std::distance(TI.begin(), std::max_element(TI.begin(), TI.end()));
        for (long long y = 0; y < (long long)data.get_size(1); y++) {
            for (long long x = 0; x < (long long)data.get_size(0); x++) {
                auto reference_phase = std::arg(data(x, y, max_TI_index));

                for (long long i = 0; i < (long long)data.get_size(2); i++) {
                    auto val        = data(x, y, i);
                    result(x, y, i) = std::abs(val)*std::polar(1.0f, std::arg(val) - reference_phase).real();
                }
            }
        }

        return result;
    }
namespace {

    




    hoNDArray<vector_td<float,2>> register_groups(const hoNDArray<float>& phase_corrected_data,
        const hoNDArray<float>& predicted, hoNDArray<vector_td<float,2>> vector_field, const registration_params& params) {
        using namespace Indexing;

        auto abs_corrected = abs(phase_corrected_data);
        auto abs_predicted = abs(predicted);

#pragma omp parallel for
        for (long long cha = 0; cha < (long long)predicted.get_size(2); cha++) {
            vector_field(slice,slice,cha) = Registration::diffeomorphic_demons<float, 2>(
                abs_corrected(slice, slice, cha), abs_predicted(slice, slice, cha),vector_field(slice,slice,cha), params.iterations,params.regularization_sigma, params.step_size,params.noise_sigma);
        }

        return vector_field;
    }
}


hoNDArray<std::complex<float>> Gadgetron::T1::deform_groups(const hoNDArray<std::complex<float>>& data,const hoNDArray<vector_td<float,2>>& vector_field){
    using namespace Indexing;
    assert(data.dimensions() == vector_field.dimensions());
    assert(data.dimensions().size() == 3);

    auto output = hoNDArray<std::complex<float>>(data.dimensions());
    for (long long cha = 0; cha < (long long)data.get_size(2); cha++){
        output(slice,slice,cha) = Registration::deform_image_bspline<std::complex<float>,2,float>(data(slice,slice,cha),vector_field(slice,slice,cha));
    }
    return output;

}



hoNDArray<vector_td<float,2>> Gadgetron::T1::t1_registration(
    const hoNDArray<std::complex<float>>& data, const std::vector<float>& TI, unsigned int iterations, registration_params params) {
    if (data.get_size(2) != TI.size()) {
        throw std::runtime_error("Data and TI do not match");
    }

    auto corrected = phase_correct(data, TI);

    auto corrected_orig = corrected;

    auto vector_field = hoNDArray<vector_td<float,2>>(data.dimensions());
    vector_field.fill(vector_td<float,2>(0));

    for (int i = 0; i < iterations; i++) {
        auto parameters = fit_T1_2param(corrected, TI);
        auto predicted = predict_signal(parameters, TI);
        auto updated_vector_field = register_groups(corrected_orig, predicted, vector_field, params);
        vector_field = updated_vector_field;
        auto deformed_data = deform_groups(data,vector_field);
        corrected = phase_correct(deformed_data, TI);
    }
    return vector_field;
};
