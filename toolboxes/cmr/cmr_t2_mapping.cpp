/** \file   cmr_t2_mapping.cpp
    \brief  Implement CMR T1 mapping for 2D acquisition
    \author Hui Xue
*/

#include "cmr_t2_mapping.h"
#include "log.h"

#include "hoNDArray_reductions.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_math.h"
#include "hoNDArray_linalg.h"

#include "simplexLagariaSolver.h"
#include "twoParaExpDecayOperator.h"
#include "curveFittingCostFunction.h"

#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron {

template <typename T>
CmrT2Mapping<T>::CmrT2Mapping() : BaseClass()
{
    max_iter_ = 150;
    max_fun_eval_ = 1000;
    thres_fun_ = 1e-4;

    max_map_value_ = 2500;
}

template <typename T>
CmrT2Mapping<T>::~CmrT2Mapping()
{
}

template <typename T>
void CmrT2Mapping<T>::get_initial_guess(const VectorType& ti, const VectorType& yi, VectorType& guess)
{
    if (guess.size() != this->get_num_of_paras())
    {
        guess.resize(this->get_num_of_paras(), 0);
    }

    // do a log linear fit
    size_t numpts = yi.size();
    GADGET_CHECK_THROW(numpts>0);

    T default_A = *std::max_element(yi.begin(), yi.end());
    T default_T2 = ti[ti.size() / 2];;

    try
    {
        size_t i;
        for (i = 0; i < numpts; i++)
        {
            if(yi[i]<=0)
            {
                guess[0] = default_A;
                guess[1] = default_T2;
                return;
            }
        }

        hoNDArray<T> log_yi(numpts), bi_est(2);
        hoNDArray<T> xi(numpts);
        for (size_t i = 0; i < numpts; i++)
        {
            xi[i] = ti[i];
            log_yi[i] = std::log(yi[i]);
        }

        T a, b;
        Gadgetron::linFit(xi, log_yi, a, b);

        guess[0] = std::exp(b);
        guess[1] = -1.0 / a;
    }
    catch(...)
    {
        guess[0] = default_A;
        guess[1] = default_T2;
    }

    if(guess[1]<0)
    {
        guess[0] = default_A;
        guess[1] = default_T2;
    }
}

template <typename T>
void CmrT2Mapping<T>::compute_map(const VectorType& ti, const VectorType& yi, const VectorType& guess, VectorType& bi, T& map_v)
{
    try
    {
        bi = guess;
        map_v = 0;

        typedef Gadgetron::twoParaExpDecayOperator< std::vector<T> > SignalType;
        typedef Gadgetron::leastSquareErrorCostFunction< std::vector<T> > CostType;

        // define solver
        Gadgetron::simplexLagariaSolver< VectorType, SignalType, CostType > solver;

        // define signal model
        SignalType t2;

        // define cost function
        CostType lse;

        solver.signal_model_ = &t2;
        solver.cf_ = &lse;

        solver.max_iter_ = max_iter_;
        solver.max_fun_eval_ = max_fun_eval_;
        solver.thres_fun_ = thres_fun_;

        solver.x_ = ti;
        solver.y_ = yi;

        solver.solve(bi, guess);

        if (bi[0] > 0 && bi[1] > 0)
        {
            map_v = bi[1];
            if (map_v >= max_map_value_) map_v = hole_marking_value_;
            if (map_v <= min_map_value_) map_v = hole_marking_value_;
        }
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in CmrT2Mapping<T>::compute_map(...) ... ");
    }
}

template <typename T>
void CmrT2Mapping<T>::compute_sd(const VectorType& ti, const VectorType& yi, const VectorType& bi, VectorType& sd, T& map_sd)
{
    try
    {
        sd.clear();
        sd.resize(bi.size(), 0);

        map_sd = 0;

        typedef Gadgetron::twoParaExpDecayOperator< std::vector<T> > SignalType;
        SignalType t2;

        // compute fitting values
        VectorType y;
        t2.magnitude(ti, bi, y);

        // compute residual
        VectorType res(y), abs_res(y);

        size_t num = ti.size();
        size_t N = this->get_num_of_paras();

        size_t n;
        for (n = 0; n < num; n++)
        {
            res[n] = y[n] - yi[n];
            abs_res[n] = std::abs(res[n]);
        }

        hoNDArray<T> grad;
        grad.create(N, num);
        Gadgetron::clear(grad);

        VectorType gradVec(N);
        for (n = 0; n < num; n++)
        {
            t2.gradient(ti[n], bi, gradVec);
            memcpy(grad.begin() + n*N, &gradVec[0], sizeof(T)*N);
        }

        GADGET_CATCH_THROW(this->compute_sd_impl(ti, yi, bi, abs_res, grad, sd));

        map_sd = sd[1];
        if (map_sd > max_map_value_) map_sd = this->hole_marking_value_;
    }
    catch (...)
    {
        GADGET_THROW("Exceptions happened in CmrT2Mapping<T>::compute_map(...) ... ");
    }
}

template <typename T>
size_t CmrT2Mapping<T>::get_num_of_paras() const
{
    return 2; // A and T2
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class CmrT2Mapping< float >;

}
