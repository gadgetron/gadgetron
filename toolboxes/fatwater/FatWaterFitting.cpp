//
// Created by dch on 18/04/18.
//

#include "FatWaterFitting.h"
#include "fatwater.h"

#include  <type_traits>
#include <numeric>
#include <boost/math/constants/constants.hpp>
#include <limits>
#include "complext.h"
#include "omp.h"

template<typename... Conds>
struct and_ : std::true_type {
};

template<typename Cond, typename... Conds>
struct and_<Cond, Conds...>
        : std::conditional<Cond::value, and_<Conds...>,
                std::false_type>::type {
};


using namespace Gadgetron;
using namespace std::complex_literals;
constexpr double PI = boost::math::constants::pi<double>();

#include <ceres/ceres.h>


using namespace Gadgetron;


namespace {

    struct complex_residual {
        template<class T>
        void operator()(T *e, const complext<T> &predicted, const complext<double> &data, int index) const {
            auto residual = predicted - data;
            e[2 * index] = residual.real();
            e[2 * index + 1] = residual.imag();
        }
    };


    struct abs_residual {
        template<class T>
        void operator()(T *e, const complext<T> &predicted, const complext<double> &data, int index) const {
            auto residual = abs(predicted) - abs(data);
            e[index] = residual;
        }
    };

    template<unsigned int NRESIDUALS, class RESIDUAL, class... SPECIES>
    class FatWaterModelCeres {
    public:

        static constexpr size_t NSPECIES = sizeof...(SPECIES);

        static constexpr size_t NPARAMS = NSPECIES * 2 + 2;

        FatWaterModelCeres(const FatWater::Parameters &parameters, const std::vector<float> &TEs,
                           const std::vector<complext<double>> &data, float fieldstrength, float r2star) : TEs_(TEs),
                                                                                                           data_(data),
                                                                                                           fieldstrength_(
                                                                                                                   fieldstrength),
                                                                                                           residual_function(),
                                                                                                           r2star_(r2star) {
            assert(parameters.species.size() == NSPECIES);
            assert(TEs.size() == data_.size());
            assert(data_.size() == NRESIDUALS);


            auto calc_amp = [&](auto &species, auto TE) {
                return std::accumulate(species.amplitude_frequency_pairs.begin(),
                                       species.amplitude_frequency_pairs.end(),
                                       complext<double>(0.0),
                                       [&](auto val, auto peak) {
                                           return val + complext<float>(peak.first) * exp(complext<double>(0.0, 2.0) *
                                                                                          (PI * peak.second *
                                                                                           fieldstrength_ *
                                                                                           parameters.gyromagnetic_ratio_Mhz *
                                                                                           TE));
                                       });
            };

            for (int k = 0; k < NSPECIES; k++) {
                for (int j = 0; j < NRESIDUALS; j++) {
                    const float TE = TEs[j];
                    species_amp[k][j] = calc_amp(parameters.species[k], TE);
                }
            }

        };


        template<class... Ts>
        bool operator()(const ceres::Jet<double, NPARAMS-1> *const fm_ptr,  Ts &&... args) const {
            return internal_implementation<ceres::Jet<double, NPARAMS-1>, float, ceres::Jet<SPECIES, NPARAMS-1>...>(fm_ptr,&this->r2star_,
                                                                                                         std::forward<Ts>(
                                                                                                                 args)...);
        }

        template<class... Ts>
        std::enable_if_t<sizeof...(Ts) == NSPECIES+1,bool>
        operator()(const double *const fm_ptr, Ts &&... args) const {
            return internal_implementation<double,float, SPECIES...>(fm_ptr,&r2star_, std::forward<Ts>(args)...);
        }


        template<class... Ts>
        bool operator()(const ceres::Jet<double, NPARAMS> *const fm_ptr, const ceres::Jet<double, NPARAMS> *const r2star, Ts &&... args) const {
            return internal_implementation<ceres::Jet<double, NPARAMS>, ceres::Jet<double, NPARAMS>, ceres::Jet<SPECIES, NPARAMS>...>(fm_ptr,r2star,
                                                                                                         std::forward<Ts>(
                                                                                                                 args)...);
        }

        template<class... Ts>
        std::enable_if_t<sizeof...(Ts) == NSPECIES+1,bool> operator()(const double *const fm_ptr,const double* const r2star, Ts &&... args) const {
            return internal_implementation<double,double, SPECIES...>(fm_ptr,r2star, std::forward<Ts>(args)...);
        }


    private:

        const std::vector<float> TEs_;
        const std::vector<complext<double>> data_;
        const float fieldstrength_;
        const float r2star_;
        RESIDUAL residual_function;
        std::array<std::array<complext<double>, NRESIDUALS>, NSPECIES> species_amp;


        template<class T, class... Ts>
        auto get_last(T, Ts &&... args) const {
            return get_last(args...);
        };

        template<class T>
        T get_last(T t) const { return t; }

        template<class T>
        struct base_type {
            using type=T;
        };
        template<class T>
        struct base_type<const T *const &> {
            using type=T;
        };
        template<class T>
        struct base_type<T *&> {
            using type=T;
        };

        template<class T> using base_type_t = typename base_type<T>::type;

        template<class T,class R, class... Ts>
        bool internal_implementation(const base_type_t<T> *const fm_ptr,const R* const r2star_ptr, const Ts *const ... args,
                                     base_type_t<T> *residual) const {


            const T &fm = *fm_ptr;
            const R &r2star = *r2star_ptr;
            auto species = std::array<complext<T>, NSPECIES>{extract_species(args)...};

            for (int j = 0; j < NRESIDUALS; j++) {
                const float TE = TEs_[j];
                complext<T> predicted = complext<T>(T(0));
                for (int i = 0; i < NSPECIES; i++)
                    predicted += species[i] * species_amp[i][j];

                predicted *= exp((T(-r2star) + complext<T>(T(0), T(2 * PI)) * fm) * T(TE));
                residual_function(residual, predicted, data_[j], j);

            }

            return true;


        }


        template<class T>
        static constexpr complext<T> extract_species(const T *const spec_ptr) {
            return complext<T>(spec_ptr[0], spec_ptr[1]);
        }


    };


    struct DiffLoss {
        DiffLoss(double scale1) : scale1_(scale1) {}

        template<class T>
        bool operator()(const T *const base, const T *const dx, T *residual) const {

            residual[0] = scale1_ * (base[0] - dx[0]);

            return true;
        }

    private:
        const double scale1_;

    };


    struct DiagLoss{
        DiagLoss(double scale1, double original) : scale1_(scale1), original_(original) {}

        template<class T>
        bool operator()(const T *const dx, T *residual) const {

            residual[0] = scale1_ * (original_ - dx[0]);

            return true;
        }

    private:
        const double scale1_;
        const double original_;

    };


    static void
    add_regularization(ceres::Problem &problem, hoNDArray<double> &field_map, float lambda,
                       ceres::LossFunction *loss = NULL) {


        auto add_term = [&](int x1, int y1,int z1,  int x2, int y2, int z2) {
            auto cost_function = new ceres::AutoDiffCostFunction<DiffLoss, 1, 1, 1>(new DiffLoss(lambda));
            std::vector<double *> ptrs = {&field_map(x1, y1,z1), &field_map(x2, y2,z2)};
            problem.AddResidualBlock(cost_function, loss, ptrs);
        };
        const size_t X = field_map.get_size(0);
        const size_t Y = field_map.get_size(1);
        const size_t Z = field_map.get_size(2);

        for (int kz = 0; kz < Z; kz++) {
            for (int ky = 0; ky < Y; ky++) {
                for (int kx = 0; kx < X; kx++) {

                    if (kx < X - 1) {
                        add_term(kx, ky, kz, kx + 1, ky, kz);
                    }
                    if (ky < Y - 1) {
                        add_term(kx, ky, kz, kx, ky + 1, kz);
                    }
                    if (kz < Z - 1) {
                        add_term(kx, ky, kz, kx, ky, kz + 1);
                    }

                }
            }
        }
    }


    static void
    add_regularization(ceres::Problem &problem, hoNDArray<double> &field_map, const hoNDArray<float> &lambda_map,
                       ceres::LossFunction *loss = NULL) {


        auto add_term = [&](int x1, int y1,int z1,  int x2, int y2, int z2) {
            auto weight = std::min(lambda_map(x1, y1,z1), lambda_map(x2, y2,z2));
            auto cost_function = new ceres::AutoDiffCostFunction<DiffLoss, 1, 1, 1>(new DiffLoss(weight));
            std::vector<double *> ptrs = {&field_map(x1, y1,z1), &field_map(x2, y2,z2)};
            problem.AddResidualBlock(cost_function, loss, ptrs);
        };
        const size_t X = field_map.get_size(0);
        const size_t Y = field_map.get_size(1);
        const size_t Z = field_map.get_size(2);

        for (int kz = 0; kz < Z; kz++) {
            for (int ky = 0; ky < Y; ky++) {
                for (int kx = 0; kx < X; kx++) {

                    if (kx < X - 1) {
                        add_term(kx, ky, kz, kx + 1, ky, kz);
                    }
                    if (ky < Y - 1) {
                        add_term(kx, ky, kz, kx, ky + 1, kz);
                    }
                    if (kz < Z - 1) {
                        add_term(kx, ky, kz, kx, ky, kz + 1);
                    }

                }
            }
        }
    }





    template<unsigned int ECHOES>
    void fat_water_fitting_echo(hoNDArray<float> &field_mapF, hoNDArray<float> &r2star_mapF,
                                hoNDArray<std::complex<float>> &fractionsF,
                                const hoNDArray<std::complex<float>> &input_data,
                                const hoNDArray<float> &lambda_map, const FatWater::Parameters &parameters)
                                 {

        hoNDArray<double> field_map;
        field_map.copyFrom(field_mapF);
        hoNDArray<double> r2star_map;
        r2star_map.copyFrom(r2star_mapF);
        hoNDArray<std::complex<double>> fractions;
        fractions.copyFrom(fractionsF);


        const size_t X = input_data.get_size(0);
        const size_t Y = input_data.get_size(1);
        const size_t Z = input_data.get_size(2);
        const size_t N = input_data.get_size(4);
        const size_t S = input_data.get_size(5);

        std::vector<float> TEs_repeated((S) * N);

        auto& TEs = parameters.echo_times_s;
        auto& field_strength = parameters.field_strength_T;

        for (int k3 = 0; k3 < S; k3++) {
            for (int k4 = 0; k4 < N; k4++) {
                TEs_repeated[k4 + (k3) * N] = double(TEs[k3]);
            }
        }

        ceres::Problem problem;
        ceres::Solver::Options options;
        options.minimizer_type = ceres::LINE_SEARCH;
        options.line_search_direction_type = ceres::LBFGS;
        options.linear_solver_type = ceres::ITERATIVE_SCHUR;
        options.num_threads = omp_get_max_threads();
        options.dense_linear_algebra_library_type = ceres::EIGEN;
        options.function_tolerance = 1e-6;
        options.gradient_tolerance = 1e-6;
        options.parameter_tolerance = 1e-6;


        for (int kz = 0; kz < Z; kz++) {
            for (int ky = 0; ky < Y; ky++) {
                for (int kx = 0; kx < X; kx++) {

                    std::vector<complext<double>> signal((S) * N);
                    auto &f = field_map(kx, ky, kz);
                    auto &r2 = r2star_map(kx, ky, kz);
                    std::complex<double> &water = fractions(kx, ky, kz, 0, 0, 0, 0);
                    std::complex<double> &fat = fractions(kx, ky, kz, 0, 0, 1, 0);

                    for (int k3 = 1; k3 < S; k3++) {
                        for (int k4 = 0; k4 < N; k4++) {
                            signal[k4 + (k3) * N] = input_data(kx, ky, kz, 0, k4, k3, 0);
                        }
                    }

                    auto cost_function = new ceres::AutoDiffCostFunction<FatWaterModelCeres<ECHOES, complex_residual, double, double>,
                            ECHOES * 2, 1,  2, 2>(
                            new FatWaterModelCeres<ECHOES, complex_residual, double, double>(parameters, TEs_repeated,
                                                                                             signal,
                                                                                             field_strength,
                                                                                             r2));
                    std::vector<double *> b = {&f,  (double *) &water, (double *) &fat};

                    problem.AddResidualBlock(cost_function, nullptr, b);
                }
            }
        }
        add_regularization(problem, field_map, lambda_map);
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << "Initial cost: " << summary.initial_cost << " Final cost:" << summary.final_cost << " Iterations "
                  << summary.iterations.size() << std::endl;

        field_mapF.copyFrom(field_map);
        r2star_mapF.copyFrom(r2star_map);
        fractionsF.copyFrom(fractions);


    }
}

void Gadgetron::FatWater::fat_water_fitting(hoNDArray<float> &field_mapF, hoNDArray<float> &r2star_mapF,
                                            hoNDArray<std::complex<float>> &fractionsF,
                                            const hoNDArray<std::complex<float>> &input_data,
                                            const hoNDArray<float> &lambda_map, const Parameters &parameters) {

    switch (parameters.echo_times_s.size()){
        case 3:
            fat_water_fitting_echo<3>(field_mapF,r2star_mapF,fractionsF,input_data,lambda_map,parameters);
            break;
        case 4:
            fat_water_fitting_echo<4>(field_mapF,r2star_mapF,fractionsF,input_data,lambda_map,parameters);
            break;
        case 5:
            fat_water_fitting_echo<5>(field_mapF,r2star_mapF,fractionsF,input_data,lambda_map,parameters);
            break;
        case 6:
            fat_water_fitting_echo<6>(field_mapF,r2star_mapF,fractionsF,input_data,lambda_map,parameters);
            break;
        case 7:
            fat_water_fitting_echo<7>(field_mapF,r2star_mapF,fractionsF,input_data,lambda_map,parameters);
            break;
        case 8:
            fat_water_fitting_echo<8>(field_mapF,r2star_mapF,fractionsF,input_data,lambda_map,parameters);
            break;
        default:
            throw std::invalid_argument("Fat water ftting only supported for 3 to 8 echoes");
    }


}




