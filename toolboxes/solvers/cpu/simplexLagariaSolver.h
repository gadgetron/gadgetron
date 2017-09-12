/** \file       simplexLagariaSolver.h
    \brief      Implement the optimizer for the simplex method

                ref: Jeffrey C.Lagarias, James A.Reeds, Margaret H.Wright, Paul E.Wright, "Convergence Properties of the Nelder-Mead Simplex Method in Low Dimensions", SIAM Journal of Optimization, 9(1), 112-147, 1998.

    \author     Hui Xue
*/

#pragma once

#include "curveFittingSolver.h"

namespace Gadgetron { 

struct simplexLagariaSolverCompObj
{
    simplexLagariaSolverCompObj() {}
    ~simplexLagariaSolverCompObj() {}

    bool operator()(const std::pair<double, size_t>& m1, const std::pair<double, size_t>& m2) const
    {
        return !(m1.first >= m2.first);
    }
};

template <typename Array_Type, typename Singal_Type, typename Cost_Type>
class simplexLagariaSolver : public curveFittingSolver<Array_Type, Singal_Type, Cost_Type>
{
public:

    typedef simplexLagariaSolver<Array_Type, Singal_Type, Cost_Type> Self;
    typedef curveFittingSolver<Array_Type, Singal_Type, Cost_Type> BaseClass;

    typedef typename Array_Type::value_type ValueType;
    typedef ValueType T;
    typedef typename realType<ValueType>::Type value_type;

    simplexLagariaSolver(double thres_x=1e-4, double thres_fun=1e-4, size_t maxIter=600, size_t maxFunc=600);
    virtual ~simplexLagariaSolver();

    // Given the initial parametes bi, solve the problem and return best parameters back in b
    virtual void solve(Array_Type& b, const Array_Type& bi);

    /// threshold for minimal variable changes
    double thres_x_;
    /// threshold for minimal function value changes
    double thres_fun_;
    /// number of maximal iteration
    size_t max_iter_;
    /// number of maximal times of function evalution
    size_t max_fun_eval_;
    ///current number of iteration
    size_t iter_;
    /// current number of function evaluation
    size_t func_evals_;
    /// current best cost
    double best_cost_;

    using BaseClass::signal_model_;
    using BaseClass::cf_;
    using BaseClass::bi_;
    using BaseClass::x_;
    using BaseClass::y_;

protected:

    using BaseClass::output_mode_;
    using BaseClass::x0_;
    using BaseClass::y_est_;

    long long num_pt_;
    long long num_;

    Array_Type pt_try_;

};

template <typename Array_Type, typename Singal_Type, typename Cost_Type>
simplexLagariaSolver<Array_Type, Singal_Type, Cost_Type>::
simplexLagariaSolver(double thres_x, double thres_fun,size_t maxIter, size_t maxFunc) : BaseClass(), thres_x_(thres_x), thres_fun_(thres_fun), max_iter_(maxIter), max_fun_eval_(maxFunc), iter_(0), func_evals_(0)
{
}

template <typename Array_Type, typename Singal_Type, typename Cost_Type>
simplexLagariaSolver<Array_Type, Singal_Type, Cost_Type>::
~simplexLagariaSolver()
{
}

template <typename Array_Type, typename Singal_Type, typename Cost_Type>
void simplexLagariaSolver<Array_Type, Singal_Type, Cost_Type>::
solve(Array_Type& b, const Array_Type& bi)
{
    try
    {
        if (signal_model_ == NULL || cf_==NULL)
        {
            GADGET_THROW("simplexLagariaSolver solver has to have the signal model and cost function ... ");
        }

        T rho = 1; 
        T chi = 2; 
        T psi = (T)0.5; 
        T sigma = (T)0.5;

        // store the initial parameters
        bi_ = bi;

        size_t n = bi_.size();

        Array_Type xin(bi_);

        hoNDArray<T> v(n, n + 1);

        Array_Type fv(n+1, 0);

        memcpy(&v(0, 0), &bi_[0], sizeof(T)*n);

        fv[0] = this->func(bi_);

        func_evals_ = 1;
        size_t itercount = 0;
        bool is_shrink = false;

        T usual_delta = (T)0.05;
        T zero_term_delta = (T)0.00025;

        Array_Type x(n), y(n);

        size_t i, j;
        for (j = 0; j < n; j++)
        {
            y = xin;
            if (y[j] != 0)
                y[j] = (1 + usual_delta)*y[j];
            else
                y[j] = zero_term_delta;

            memcpy(&v(0, j + 1), &y[0], sizeof(T)*n);

            x = y; 
            T f = this->func(x);
            fv[j+1] = f;
        }

        std::vector< std::pair<double, size_t> > fvSorted(n + 1);
        for (j = 0; j <= n; j++)
        {
            fvSorted[j].first = fv[j];
            fvSorted[j].second = j;
        }

        std::sort(fvSorted.begin(), fvSorted.end(), simplexLagariaSolverCompObj());

        hoNDArray<T> vSorted(v);

        for (j = 0; j <= n; j++)
        {
            fv[j] = fvSorted[j].first;

            for (i = 0; i < n; i++)
            {
                vSorted(i, j) = v(i, fvSorted[j].second);
            }
        }
        v = vSorted;

        itercount = itercount + 1;
        func_evals_ = n + 1;

        int exitflag = 1;

        hoNDArray<T> fvDiff(n);
        hoNDArray<T> vDiff(n, n);

        Array_Type xbar(n, 0), xr(n, 0), xe(n, 0), xc(n, 0), xcc(n, 0);

        T fxe(0), fxc(0), fxcc(0);

        while ((func_evals_<max_fun_eval_) && (itercount<max_iter_))
        {
            for (j = 1; j <= n; j++)
            {
                fvDiff(j - 1) = std::abs(fv[0] - fv[j]);
            }

            T fDiff;
            size_t ind;
            Gadgetron::maxAbsolute(fvDiff, fDiff, ind);

            // check the thres_x_
            for (j = 1; j <= n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    vDiff(i, j - 1) = std::abs( v(i, j) - v(i, 0) );
                }
            }

            T xDiff;
            Gadgetron::maxAbsolute(vDiff, xDiff, ind);

            if (this->output_mode_ >= Self::OUTPUT_VERBOSE)
            {
                GDEBUG_STREAM("--> simplexLagariaSolver, itercount = " << itercount << " - fDiff : " << fDiff << " - xDiff : " << xDiff);
            }

            if ((fDiff < thres_fun_) && (xDiff < thres_x_))
            {
                break;
            }

            for (j = 0; j < n; j++)
            {
                xbar[j] = 0;
                for (i = 0; i < n; i++)
                {
                    xbar[j] += v(j, i);
                }
            }

            for (i = 0; i < n; i++)
            {
                xbar[i] /= (T)n;
            }

            for (i = 0; i < n; i++)
            {
                xr[i] = (1 + rho)*xbar[i] - rho*v(i, n);
            }

            bi_ = xr;
            T fxr = this->func(bi_);
            func_evals_ = func_evals_ + 1;

            if (fxr < fv[0])
            {
                for (i = 0; i < n; i++)
                {
                    xe[i] = (1 + rho*chi)*xbar[i] - rho*chi*v(i, n);
                }

                bi_ = xe;
                fxe = this->func(bi_);
                func_evals_ = func_evals_ + 1;

                if (fxe<fxr)
                {
                    memcpy(v.begin() + n*n, &xe[0], sizeof(T)*n);
                    fv[n] = fxe;
                }
                else
                {
                    memcpy(v.begin() + n*n, &xr[0], sizeof(T)*n);
                    fv[n] = fxr;
                }

            }
            else
            {
                if (fxr < fv[n - 1])
                {
                    memcpy(v.begin() + n*n, &xr[0], sizeof(T)*n);
                    fv[n] = fxr;
                }
                else
                {
                    if (fxr < fv[n])
                    {
                        for (i = 0; i < n; i++)
                        {
                            xc[i] = (1 + psi*rho)*xbar[i] - psi*rho*v(i, n);
                        }

                        bi_ = xc;
                        fxc = this->func(bi_);
                        func_evals_ = func_evals_ + 1;

                        if (fxc<=fxr)
                        {
                            memcpy(v.begin() + n*n, &xc[0], sizeof(T)*n);
                            fv[n] = fxc;
                        }
                        else
                        {
                            is_shrink = true;
                        }
                    }
                    else
                    {
                        for (i = 0; i < n; i++)
                        {
                            xcc[i] = (1 - psi)*xbar[i] + psi*v(i, n);
                        }

                        bi_ = xcc;
                        fxcc = this->func(bi_);
                        func_evals_ = func_evals_ + 1;

                        if (fxcc <= fv[n])
                        {
                            memcpy(v.begin() + n*n, &xcc[0], sizeof(T)*n);
                            fv[n] = fxcc;
                        }
                        else
                        {
                            is_shrink = true;
                        }
                    }

                    if (is_shrink)
                    {
                        for (j = 1; j<=n; j++)
                        {
                            for (i = 0; i < n; i++)
                            {
                                v(i, j) = v(i, 0) + sigma*(v(i, j) - v(i, 0));
                                bi_[i] = v(i, j); 
                            }
                            fv[j] = this->func(bi_);
                        }
                        func_evals_ = func_evals_ + n;
                    }
                }
            }

            for (j = 0; j <= n; j++)
            {
                fvSorted[j].first = fv[j];
                fvSorted[j].second = j;
            }

            std::sort(fvSorted.begin(), fvSorted.end(), simplexLagariaSolverCompObj());

            for (j = 0; j <= n; j++)
            {
                fv[j] = fvSorted[j].first;

                for (i = 0; i < n; i++)
                {
                    vSorted(i, j) = v(i, fvSorted[j].second);
                }
            }
            v = vSorted;

            itercount = itercount + 1;
        }

        this->iter_ = itercount;

        for (i = 0; i < n; i++)
        {
            bi_[i] = v(i, 0);
        }

        best_cost_ = fv[0];

        b = bi_;
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in simplexLagariaSolver<Array_Type, Func_Type>::solve(...) ... ");
    }
}

}
