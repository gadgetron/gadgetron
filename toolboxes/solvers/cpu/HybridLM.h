//
// Created by dchansen on 2/27/20.
//

#pragma once
#define ARMA_DONT_PRINT_ERRORS
#include "hoArmadillo.h"
namespace Gadgetron { namespace Solver {

    enum class ReturnStatus { SUCCESS, MAX_ITERATIONS_REACHED, LINEAR_SOLVER_FAILED };
    template <class Scalar> class HybridLMSolver {
    public:
        HybridLMSolver(size_t num_residuals, size_t num_params)
            : B{ arma::Mat<Scalar>(num_params, num_params) }
            , J{ arma::Mat<Scalar>(num_residuals, num_params) }
            , J_new{ arma::Mat<Scalar>(num_residuals, num_params) }
            , residuals{ arma::Col<Scalar>(num_residuals) }
            , JTJ{ arma::Mat<Scalar>(num_params, num_params) }
            , DTD{ arma::Col<Scalar>(num_params) } {}

        template <class F> ReturnStatus solve(F& func, arma::Col<Scalar>& params) {
            using namespace arma;
            constexpr int max_iterations = 1000;
            B.fill(fill::eye);

            DTD.fill(-std::numeric_limits<Scalar>::infinity());
            count = 0;

            func(params, residuals, J);

            JTJ = J.t() * J;

            mu               = 1e-4;
            v                = 2;
            bool use_lm_step = true;

            for (long long k = 0; k < max_iterations; k++) {
                auto new_params    = params;
                auto old_residuals = residuals;
                bool better        = false;
                StepStatus result;
                if (use_lm_step) {
                    result = lm_step(func, new_params, better, use_lm_step);
                } else {
                    result = qm_step(func, new_params, better, use_lm_step);
                }
                if (result == StepStatus::CONVERGED) return ReturnStatus::SUCCESS;
                if (result == StepStatus::LINEAR_SOLVER_FAILED) return ReturnStatus::LINEAR_SOLVER_FAILED;

                JTJ = J_new.t() * J_new;

                Col<Scalar> h = new_params - params;
                Col<Scalar> y = JTJ * h + (J_new - J).t() * residuals;
                Scalar hy     = dot(h, y);
                if (hy > 0) {
                    Col<Scalar> v    = B * h;
                    Mat<Scalar> Btmp = B + y * y.t() / dot(h, y) - v * v.t() / dot(h, v);
                    if (!Btmp.has_nan() && !Btmp.has_inf()) {
                        B = Btmp;
                    }
                }
                if (better) {
                    J      = J_new;
                    params = new_params;
                } else {
                    JTJ       = J.t() * J;
                    residuals = old_residuals;
                }
            }
            return ReturnStatus::MAX_ITERATIONS_REACHED;
        }

    private:

        enum class StepStatus { CONVERGED, CONTINUE, LINEAR_SOLVER_FAILED };
        template <class F> StepStatus lm_step(F& func, arma::Col<Scalar>& params, bool& better, bool& use_lm_step) {
            using namespace arma;
            const Col<Scalar> g = J.t() * residuals;

            for (long long i = 0; i < DTD.n_elem; i++) {
                DTD[i] = std::max(DTD[i], JTJ(i, i));
            }
            Mat<Scalar> h;
            if (!arma::solve(h, JTJ + mu * diagmat(DTD), -g, solve_opts::fast)) return StepStatus::LINEAR_SOLVER_FAILED;

            if (norm(h) < minimum_step_size * (norm(params) + minimum_step_size)) {
                return StepStatus::CONVERGED;
            }

            params += h;
            Scalar old_cost = dot(residuals, residuals) / 2;
            func(params, residuals, J_new);
            Scalar new_cost = dot(residuals, residuals) / 2;

            Scalar rho = (old_cost - new_cost) / (dot(h, mu * h - g) / 2);
            better     = rho > 0;
            if (better) {
                mu = mu * std::max(Scalar(1) / 3, Scalar(1) - std::pow(2 * rho - 1, Scalar(3)));
                v  = 2;

                Scalar gradient_norm = max(abs(J_new.t() * residuals));
                if (gradient_norm <= minimum_gradient) {
                    return StepStatus::CONVERGED;
                }
                if (gradient_norm < 0.02 * new_cost) {
                    count++;
                    if (count == 3) {
                        use_lm_step = false;
                        delta       = std::max<Scalar>(
                            1.5 * minimum_step_size * (norm(params - h) + minimum_step_size), norm(h) / 5);
                    }
                } else {
                    count = 0;
                }

            } else {
                mu    = mu * v;
                v     = 2 * v;
                count = 0;
            }
            return StepStatus::CONTINUE;
        }

        template <class F> StepStatus qm_step(F& func, arma::Col<Scalar>& params, bool& better, bool& use_lm_step) {
            using namespace arma;
            const Col<Scalar> g = J.t() * residuals;
            Col<Scalar> h;
            if (!arma::solve(h,B, -g, solve_opts::fast)) return StepStatus::LINEAR_SOLVER_FAILED;

            if (norm(h) < minimum_step_size * (norm(params) + minimum_step_size)) {
                return StepStatus::CONVERGED;
            }

            if (norm(h) > delta) {
                h *= delta / norm(h);
            }

            params += h;

            Scalar old_cost = dot(residuals, residuals) / 2;

            Scalar old_gradient_norm = max(abs(J.t() * residuals));
            func(params, residuals, J_new);
            Scalar new_cost          = dot(residuals, residuals) / 2;
            Scalar new_gradient_norm = max(abs(J_new.t() * residuals));

            if (new_gradient_norm <= minimum_gradient) {
                return StepStatus::CONVERGED;
            }

            Scalar rho = (old_cost - new_cost) / (dot(h, mu * h - g) / 2);
            if (rho < 0.25) {
                delta = delta / 2;
            } else if (rho > 0.75) {
                delta = std::max(delta, 3 * norm(h));
            }

            better = (new_cost < old_cost)
                     || ((new_cost <= (1 + std::sqrt(std::numeric_limits<Scalar>::epsilon())) * old_cost)
                         && (new_gradient_norm < old_gradient_norm));

            if (new_gradient_norm >= old_gradient_norm)
                use_lm_step = true;

            return StepStatus::CONTINUE;
        }
        arma::Mat<Scalar> B;
        arma::Mat<Scalar> J;
        arma::Mat<Scalar> J_new;
        arma::Mat<Scalar> JTJ;
        arma::Col<Scalar> DTD;
        arma::Col<Scalar> residuals;
        int count;
        Scalar mu;
        Scalar v;
        Scalar delta;

        Scalar minimum_step_size = Scalar(1e-6);
        Scalar minimum_gradient  = Scalar(1e-8);
    };

} // namespace Solver
}
