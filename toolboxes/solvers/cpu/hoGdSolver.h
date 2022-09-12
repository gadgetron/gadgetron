/** \file       hoGdSolver.h
    \brief      Implement the optimizer for the gradient descent algorithm
    \author     Hui Xue
*/

#pragma once

#include "solver.h"
#include "linearOperator.h"

namespace Gadgetron { 

template <typename Array_Type, typename Proximal_Oper_Type> class hoGdSolver;

template <typename Array_Type, typename Proximal_Oper_Type>
class hoGdSolverCallBack
{
public:

    hoGdSolverCallBack() : solver_(NULL) {}
    virtual ~hoGdSolverCallBack() {}

    typedef hoGdSolver<Array_Type, Proximal_Oper_Type> SolverType;
    SolverType* solver_;

    virtual void execute(const Array_Type& b, Array_Type& x) = 0;
};

template <typename Array_Type, typename Proximal_Oper_Type>
class hoGdSolver : public solver<Array_Type, Array_Type>
{
public:

    typedef hoGdSolver<Array_Type, Proximal_Oper_Type> Self;
    typedef solver<Array_Type, Array_Type> BaseClass;

    typedef typename Array_Type::element_type ValueType;
    typedef typename realType<ValueType>::Type value_type;

    hoGdSolver();
    virtual ~hoGdSolver();

    virtual boost::shared_ptr<Array_Type> solve(Array_Type* x);
    virtual void solve(const Array_Type& b, Array_Type& x);

    /// number of max iterations
    size_t iterations_;

    /// threshold for detla change of solution gradient
    value_type grad_thres_;

    /// threshold for detla change of objective function
    value_type thres_;

    /// record the function values
    std::vector<value_type> func_value_;

    /// strength of proximity operation
    value_type proximal_strength_ratio_;

    /// the scale factor for regularization
    /// if < 0, then compute the scale factor in the solver
    value_type scale_factor_;

    /// whether to determine lamda from L1 term
    bool determine_proximal_strength_from_L1_term_;

    /// maximal number of inner iterations
    size_t iterations_inner_;

    /// maximal number of linear search steps
    size_t search_steps_;

    linearOperator<Array_Type>* oper_system_;
    Proximal_Oper_Type* oper_reg_;

    hoGdSolverCallBack<Array_Type, Proximal_Oper_Type>* call_back_;

protected:

};

template <typename Array_Type, typename Proximal_Oper_Type>
hoGdSolver<Array_Type, Proximal_Oper_Type>::
hoGdSolver() : BaseClass()
{
    iterations_ = 100;
    grad_thres_ = (value_type)1e-5;
    thres_ = (value_type)0.1;
    proximal_strength_ratio_ = 1e-3;

    scale_factor_ = -1;
    determine_proximal_strength_from_L1_term_ = false;

    iterations_inner_ = 50;
    search_steps_ = 10;

    oper_system_ = NULL;
    oper_reg_ = NULL;

    call_back_ = NULL;
}

template <typename Array_Type, typename Proximal_Oper_Type>
hoGdSolver<Array_Type, Proximal_Oper_Type>::
~hoGdSolver()
{
}

template <typename Array_Type, typename Proximal_Oper_Type>
boost::shared_ptr<Array_Type> hoGdSolver<Array_Type, Proximal_Oper_Type>::solve(Array_Type* x)
{
    boost::shared_ptr<Array_Type> b(new Array_Type);
    this->solve(*b, *x);
    return b;
}

template <typename Array_Type, typename Proximal_Oper_Type>
void hoGdSolver<Array_Type, Proximal_Oper_Type>::
solve(const Array_Type& b, Array_Type& x)
{
    try
    {
        if (oper_system_ == NULL || oper_reg_ == NULL)
        {
            GADGET_THROW("hoAdvancedGradientDescent solver can only handle two operators ... ");
        }

        GADGET_CHECK_THROW(this->x0_ != NULL);

        func_value_.reserve(iterations_);

        Array_Type ATb;
        Array_Type* pb = const_cast<Array_Type*>(&b);
        oper_system_->mult_MH(pb, &ATb);

        Array_Type WATb;
        if (determine_proximal_strength_from_L1_term_)
        {
            oper_reg_->mult_M(&ATb, &WATb);
        }

        value_type norm_length = 0.1;
        value_type norm_max;
        size_t indMax;
        hoNDArray<value_type> magWv, magATy;

        if (this->scale_factor_ < 0)
        {
            if (determine_proximal_strength_from_L1_term_)
            {
                Gadgetron::abs(WATb, magWv);
                Gadgetron::maxAbsolute(magWv, norm_max, indMax);
            }
            else
            {
                Gadgetron::abs(ATb, magATy);
                Gadgetron::maxAbsolute(magATy, norm_max, indMax);
            }
        }
        else
        {
            norm_max = scale_factor_;
        }

        value_type proximal_strength = proximal_strength_ratio_ * std::abs(norm_max);
        if (std::abs(proximal_strength) < FLT_EPSILON)
        {
            Gadgetron::abs(*(this->x0_), magATy);
            Gadgetron::maxAbsolute(magATy, norm_max, indMax);

            proximal_strength = proximal_strength_ratio_ * std::abs(norm_max);
        }

        if (this->output_mode_ >= Self::OUTPUT_VERBOSE)
        {
            GDEBUG_STREAM("---> hoAdvancedGradientDescent iteration : proximal_strength - " << proximal_strength);
        }

        if (proximal_strength<FLT_EPSILON) return;

        x = *(this->x0_);

        Array_Type Ax;
        oper_system_->mult_M(&x, &Ax);

        Array_Type bufX(x);
        Array_Type bufAx(Ax), bufAx2(Ax);
        Array_Type bufX2(x);
        Gadgetron::clear(bufX2);

        value_type stepA = 0;
        value_type stepB = 1;

        size_t nIter;

        Array_Type x2(x), diffx(x), xprev(x);
        Array_Type diffb(ATb), diffbNorm(ATb), ATAb(ATb);
        Array_Type proximal_WATb(WATb), proximal_WTATb(WATb), proximal_res(WATb), r(WATb);
        value_type diffA_norm, diffX_norm;

        for (nIter = 0; nIter<iterations_; nIter++)
        {
            value_type tt = (stepA - 1) / stepB;

            x2 = bufX2;
            Gadgetron::scal(tt, x2);
            Gadgetron::add(x, x2, x2);

            Gadgetron::subtract(Ax, bufAx, bufAx2);
            Gadgetron::scal(tt, bufAx2);
            Gadgetron::add(Ax, bufAx2, bufAx2);

            oper_system_->mult_MH(&bufAx2, &ATAb);
            Gadgetron::subtract(ATAb, ATb, diffb);

            bufX = x;

            size_t iterInner;
            for (iterInner = 0; iterInner<iterations_inner_; iterInner++)
            {
                diffbNorm = diffb;
                Gadgetron::scal(value_type(1.0) / norm_length, diffbNorm);
                Gadgetron::subtract(x2, diffbNorm, diffx);

                value_type proximal_strength_normalized = proximal_strength / norm_length;

                oper_reg_->mult_M(&diffx, &WATb);

                if (!proximal_WATb.dimensions_equal(&WATb))
                {
                    proximal_WATb.create(WATb.dimensions());
                    Gadgetron::clear(proximal_WATb);
                }

                if (!proximal_WTATb.dimensions_equal(&WATb))
                {
                    proximal_WTATb.create(WATb.dimensions());
                    Gadgetron::clear(proximal_WTATb);
                }

                size_t N = proximal_WATb.get_number_of_elements();
                ValueType* pProximal_WATb = proximal_WATb.begin();

                long long n;
#pragma omp parallel for default(none) private(n) shared(N, proximal_strength_normalized, pProximal_WATb)
                for (n = 0; n<N; n++)
                {
                    value_type mag = std::abs(pProximal_WATb[n]);
                    if (mag < FLT_EPSILON)
                        pProximal_WATb[n] = 0;
                    else
                    {
                        if (mag > proximal_strength_normalized)
                            pProximal_WATb[n] *= (proximal_strength_normalized / mag);
                    }
                }

                Gadgetron::subtract(WATb, proximal_WATb, WATb);
                Gadgetron::subtract(WATb, proximal_WTATb, WATb);

                size_t ii;
                for (ii = 0; ii<search_steps_; ii++)
                {
                    Gadgetron::add(WATb, proximal_WATb, proximal_res);

                    oper_reg_->proximity(proximal_res, proximal_strength_normalized);

                    Gadgetron::subtract(WATb, proximal_res, WATb);
                    Gadgetron::add(proximal_WATb, WATb, proximal_WATb);

                    value_type n1 = Gadgetron::nrm2(WATb);

                    if (oper_reg_->unitary())
                    {
                        Gadgetron::add(proximal_res, proximal_WTATb, WATb);
                    }
                    else
                    {
                        Gadgetron::add(proximal_res, proximal_WTATb, WATb);
                        oper_reg_->mult_MH(&WATb, &r);
                        oper_reg_->mult_M(&r, &WATb);
                    }

                    Gadgetron::subtract(proximal_res, WATb, r);
                    Gadgetron::add(proximal_WTATb, r, proximal_WTATb);

                    value_type nx = Gadgetron::nrm2(WATb);

                    if (n1 < nx*1e-4) break;
                }

                oper_reg_->mult_MH(&WATb, &x);

                Gadgetron::subtract(x, x2, diffx);

                oper_system_->mult_M(&x, &Ax);

                Gadgetron::subtract(Ax, bufAx2, bufAx);

                diffX_norm = Gadgetron::nrm2(diffx);
                diffX_norm = diffX_norm*diffX_norm;

                diffA_norm = Gadgetron::nrm2(bufAx);
                diffA_norm = diffA_norm*diffA_norm;

                if (diffX_norm <= FLT_EPSILON) break;
                if (diffA_norm <= diffX_norm * norm_length)
                {
                    break;
                }
                else
                {
                    norm_length = std::max((value_type)(1.5*norm_length), diffA_norm / diffX_norm);
                }
            }

            bufAx = Ax;

            stepA = stepB;
            stepB = (value_type)((1.0 + std::sqrt(4.0*stepA*stepA + 1.0)) / 2.0);

            Gadgetron::subtract(x, bufX, bufX2);
            Gadgetron::subtract(Ax, b, bufAx2);

            oper_reg_->mult_M(&x, &WATb);

            value_type error_data_fidelity = Gadgetron::nrm2(bufAx2);
            error_data_fidelity = error_data_fidelity*error_data_fidelity;
            error_data_fidelity *= 0.5;

            value_type error_image_reg = Gadgetron::asum(WATb);
            func_value_.push_back(error_data_fidelity + proximal_strength*error_image_reg);

            if (this->output_mode_ >= Self::OUTPUT_VERBOSE)
            {
                if (nIter > 0)
                {
                    GDEBUG_STREAM("---> iteration " << nIter << " - cost : " << error_data_fidelity << " - delta change : " << (func_value_[nIter] - func_value_[nIter - 1]) << " - " << (func_value_[nIter] - func_value_[nIter - 1]) / func_value_[nIter - 1]);
                }
                else
                {
                    GDEBUG_STREAM("---> iteration " << nIter << " - initial cost : " << error_data_fidelity);
                }
            }

            if (nIter >= 2)
            {
                if (func_value_[nIter] > func_value_[nIter - 1])
                {
                    x = xprev;
                    break;
                }

                if (std::abs(func_value_[nIter] - func_value_[nIter - 1]) <= thres_)
                {
                    break;
                }

                if (std::abs(func_value_[nIter] - func_value_[nIter - 1]) / func_value_[nIter - 1] <= grad_thres_)
                {
                    break;
                }
            }

            if(call_back_!=NULL)
            {
                call_back_->solver_ = this;
                call_back_->execute(b, x);
            }

            xprev = x;
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in hoGdSolver<Array_Type, Proximal_Oper_Type>::solve(...) ... ");
    }
}

}
