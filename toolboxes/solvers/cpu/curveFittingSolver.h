/** \file       curveFittingSolver.h
    \brief      Implement the base solver for curve fitting purpose.

    \author     Hui Xue
*/

#pragma once

#include "solver.h"

namespace Gadgetron { 

template <typename Array_Type, typename Singal_Type, typename Cost_Type>
class curveFittingSolver : public solver<Array_Type, Array_Type>
{
public:

    typedef curveFittingSolver<Array_Type, Singal_Type, Array_Type> Self;
    typedef solver<Array_Type, Array_Type> BaseClass;

    typedef typename Array_Type::value_type ValueType;
    typedef ValueType T;
    typedef typename realType<ValueType>::Type value_type;

    curveFittingSolver();
    virtual ~curveFittingSolver();

    // Given the initial parameters bi, solve the problem and return best parameters
    virtual boost::shared_ptr<Array_Type> solve(Array_Type* bi);
    // Given the initial parametes bi, solve the problem and return best parameters back in b
    virtual void solve(Array_Type& b, const Array_Type& bi) = 0;

    // the signal model object
    Singal_Type* signal_model_;
    // the cost function object
    Cost_Type* cf_;

    // initial parameters
    Array_Type bi_;

    // measured point (x, y)
    Array_Type x_;
    Array_Type y_;

protected:

    using BaseClass::output_mode_;
    using BaseClass::x0_;

    Array_Type y_est_;

    T func(Array_Type& b)
    {
        this->signal_model_->magnitude(x_, b, y_est_);
        return this->cf_->eval(y_, y_est_);
    }
};

template <typename Array_Type, typename Singal_Type, typename Cost_Type>
curveFittingSolver<Array_Type, Singal_Type, Cost_Type>::
curveFittingSolver() : BaseClass()
{
    signal_model_ = NULL;
    cf_ = NULL;
}

template <typename Array_Type, typename Singal_Type, typename Cost_Type>
curveFittingSolver<Array_Type, Singal_Type, Cost_Type>::
~curveFittingSolver()
{
}

template <typename Array_Type, typename Singal_Type, typename Cost_Type>
boost::shared_ptr<Array_Type> curveFittingSolver<Array_Type, Singal_Type, Cost_Type>::solve(Array_Type* bi)
{
    boost::shared_ptr<Array_Type> b(new Array_Type);
    this->bi_ = *bi;
    this->solve(*b, *bi);
    return b;
}

}
