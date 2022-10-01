//
// Created by dch on 21/02/18.
//
#include "cmr_t1_mapping.h"
#include "curveFittingCostFunction.h"
#include "GadgetronTimer.h"
#include "hoNDArray_math.h"
#include "hoNDHarrWavelet.h"
#include "hoNDRedundantWavelet.h"
#include "ImageIOAnalyze.h"
#include "log.h"
#include "simplexLagariaSolver.h"
#include "twoParaExpDecayOperator.h"
#include "twoParaExpRecoveryOperator.h"

#include <chrono>
#include <gtest/gtest.h>
#include <boost/random.hpp>

#include <dlib/optimization.h>

#include <ceres/ceres.h>
#include <numeric>
#define ITERATIONS 10000


class twoParaExpRecovery  {
public:

    twoParaExpRecovery(std::vector<double> x,std::vector<double> y) : x_(x), y_(y) {}
    template <typename T>
            bool operator()(const T* const b, T* e) const {


//        T sign_b1 = (b[1])/ceres::abs(b[1]);
//        T rb = 1.0 / ( (ceres::abs(b[1])< T(FLT_EPSILON)) ? sign_b1*T(FLT_EPSILON) : b[1] );
        T rb = 1.0 /b[1];

        size_t ii;

//        T result = T(0);
        for (ii=0; ii< x_.size(); ii++)
        {
            T tmp = b[0] - b[0] * exp( -1 * x_[ii] * rb);
            tmp -= y_[ii];
//            result += tmp*tmp;
            e[ii] = tmp;
//            e[ii] = tmp*tmp;
        }
//        e[0] = result/T(x_.size());
        return true;
    }

private:
    std::vector<double> x_;
    std::vector<double> y_;
};


void time_ceres(){
    std::vector<double> b;

 auto start = std::chrono::high_resolution_clock::now();
    for (int  i = 0; i < ITERATIONS; i++) {
        ceres::Problem problem;

        std::vector<double> y = {178, 185, 182, 189, 178, 180, 187, 179, 177, 177, 471};
        auto x = std::vector<double>(11, 545);
        x[10] = 10000;
        b = {*std::max_element(y.begin(), y.end()), x[5]};
        auto cost_function = new ceres::AutoDiffCostFunction<twoParaExpRecovery, 11, 2>(
                new twoParaExpRecovery(x, y));
//        auto cost_function = new ceres::NumericDiffCostFunction<twoParaExpRecovery,ceres::RIDDERS,ceres::DYNAMIC,2>(
//                new twoParaExpRecovery(x,y),ceres::DO_NOT_TAKE_OWNERSHIP,1);
        problem.AddResidualBlock(cost_function, NULL, b.data());

        ceres::Solver::Options options;
//        options.max_num_iterations = 15000;
        options.linear_solver_type = ceres::DENSE_QR;
//        options.use_explicit_schur_complement = true;
    options.function_tolerance = 1e-4;
    options.gradient_tolerance = 1e-4;
        options.parameter_tolerance = 1e-4;
//        options.preconditioner_type = ceres::IDENTITY;
//    options.minimizer_type = ceres::LINE_SEARCH;
//        options.line_search_direction_type = ceres::BFGS;
//        options.trust_region_strategy_type = ceres::DOGLEG;
//    options.dogleg_type = ceres::SUBSPACE_DOGLEG;


        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
//    GINFO_STREAM(summary.FullReport() << std::endl);

//        GINFO_STREAM("B" << b[0] << " " << b[1] << std::endl);

    }

     typedef Gadgetron::twoParaExpRecoveryOperator<std::vector<double> > SignalType;
        typedef Gadgetron::leastSquareErrorCostFunction<std::vector<double> > CostType;

        // define solver

        // define signal model
        SignalType t1_sr;

        // define cost function
        CostType lse;
 std::vector<double> y = {178, 185, 182, 189, 178, 180, 187, 179, 177, 177, 471};
        auto x = std::vector<double>(11, 545);
        x[10] = 10000;
      auto cost_function2 = [&y, &x, &t1_sr,&lse](std::vector<double> b) {
//            auto b2 = std::vector<double>{b(0), b(1)};
            auto result = std::vector<double>(2);
            t1_sr.magnitude(x, b, result);

            double tmp = lse.eval(y,result);
            return tmp;


        };

    GINFO_STREAM("Cost " <<  cost_function2(b) << std::endl);

    auto end = std::chrono::high_resolution_clock::now();

    GINFO_STREAM("Fitting tookz " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl);
    GINFO_STREAM("B " << b[0] << " " << b[1] << std::endl);

}
void time_dlib(){

    double best_cost;

    auto start = std::chrono::high_resolution_clock::now();

    dlib::matrix<double, 2, 1> b;
    for (auto i  = 0; i < ITERATIONS; i++) {
        typedef Gadgetron::twoParaExpRecoveryOperator<std::vector<double> > SignalType;
        typedef Gadgetron::leastSquareErrorCostFunction<std::vector<double> > CostType;

        // define solver

        // define signal model
        SignalType t1_sr;

        // define cost function
        CostType lse;
//
//    solver.signal_model_ = &t1_sr;
//    solver.cf_ = &lse;
//
//    solver.max_iter_ = 150;
//    solver.max_fun_eval_ = 1000;
//    solver.thres_fun_ = 1e-4;
//
//    // set measured points
//    solver.x_.resize(11, 545); // echo time, in ms
//    solver.x_[10] = 10000;

        std::vector<double> y = {178, 185, 182, 189, 178, 180, 187, 179, 177, 177, 471};
        auto x = std::vector<double>(11, 545);
        x[10] = 10000;


        auto cost_function = [&y, &x, &t1_sr,&lse](const dlib::matrix<double, 2, 1> b) {
            auto b2 = std::vector<double>{b(0), b(1)};
            auto result = std::vector<double>(2);
            t1_sr.magnitude(x, b2, result);

            double tmp = lse.eval(y,result);
            return tmp;


        };

        b(0) = *std::max_element(y.begin(), y.end());
        b(1) = x[5];

//        dlib::find_min_bobyqa(cost_function, b, 4, dlib::uniform_matrix<double>(2, 1, 0),
//                              dlib::uniform_matrix<double>(2, 1, 1e5), 100, 1e-4, 150);
        dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),dlib::objective_delta_stop_strategy(1e-4),cost_function,b,-1);
//        dlib::find_min_
//        auto result = dlib::find_min_global(cost_function,{0,0},{1e5,1e5},dlib::max_function_calls(10000));
        best_cost = cost_function(b);

    }
    auto end = std::chrono::high_resolution_clock::now();

    GINFO_STREAM("Fitting tookz " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl);

    GINFO_STREAM("Best cost " << best_cost << " " << b << std::endl);
}
void time_gadgetron(){

    double best_cost;
    std::vector<float> b(2);
    auto start = std::chrono::system_clock::now();
    for (auto i = 0; i < ITERATIONS; i++) {
        typedef Gadgetron::twoParaExpRecoveryOperator<std::vector<float> > SignalType;
        typedef Gadgetron::leastSquareErrorCostFunction<std::vector<float> > CostType;

        // define solver
        Gadgetron::simplexLagariaSolver<std::vector<float>, SignalType, CostType> solver;

        // define signal model
        SignalType t1_sr;

        // define cost function
        CostType lse;

        solver.signal_model_ = &t1_sr;
        solver.cf_ = &lse;

        solver.max_iter_ = 1500;
        solver.max_fun_eval_ = 10000;
        solver.thres_fun_ = 1e-6;

        // set measured points
        solver.x_.resize(11, 545); // echo time, in ms
        solver.x_[10] = 10000;

        solver.y_.resize(11); // intensity
        solver.y_[0] = 178;
        solver.y_[1] = 185;
        solver.y_[2] = 182;
        solver.y_[3] = 189;
        solver.y_[4] = 178;
        solver.y_[5] = 180;
        solver.y_[6] = 187;
        solver.y_[7] = 179;
        solver.y_[8] = 177;
        solver.y_[9] = 177;
        solver.y_[10] = 471;

        std::vector<float> guess(2, 0);
        b[0] = 0;
        b[1] = 0;

        guess[0] = *std::max_element(solver.y_.begin(), solver.y_.end());
        guess[1] = solver.x_[solver.x_.size() / 2];

        solver.solve(b, guess);
        best_cost = solver.best_cost_;
    }
    auto end = std::chrono::system_clock::now();

    GINFO_STREAM("Fitting tookz " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl);
    GINFO_STREAM("Best cost " << best_cost << " " << b[0] << " " << b[1] <<  std::endl);
}
using namespace Gadgetron;
int main(){
    time_gadgetron();
    time_dlib();
    time_ceres();
}
