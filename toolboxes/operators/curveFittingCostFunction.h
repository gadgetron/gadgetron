/** \file       curveFittingCostFunction.h
    \brief      Base class for cost function for the curve fitting purpose.
    \author     Hui Xue
*/

#pragma once

#include "log.h"
#include <cmath>

namespace Gadgetron {

    template <class ARRAY> class curveFittingCostFunction
    {
    public:

        typedef curveFittingCostFunction<ARRAY> Self;
        typedef typename ARRAY::value_type ELEMENT_TYPE;
        typedef typename realType<ELEMENT_TYPE>::Type REAL;

        curveFittingCostFunction() {}

        // evaluates cost function, given observed (y) and estimated (y_est) data
        virtual REAL eval(const ARRAY &y, const ARRAY& y_est) = 0;
        virtual ~curveFittingCostFunction() = default;
    };

    // simple LSE
    template <class ARRAY> class leastSquareErrorCostFunction : public curveFittingCostFunction< ARRAY >
    {
    public:
        typedef curveFittingCostFunction<ARRAY> BaseClass;
        typedef typename BaseClass::REAL REAL;

        leastSquareErrorCostFunction() {}
        virtual REAL eval(const ARRAY &y, const ARRAY& y_est)
        {
            GADGET_CHECK_THROW(y.size()==y_est.size());

            REAL cost(0);

            size_t num = y.size();
            for (size_t ii=0; ii<num; ii++)
            {
                REAL v = std::abs(y_est[ii] - y[ii]);
                cost += v*v;
            }

            if(num>1) cost /= num;

            return cost;
        }
    };

    // weighted LSE
    template <class ARRAY> class weightedLeastSquareErrorCostFunction : public curveFittingCostFunction< ARRAY >
    {
    public:
        typedef curveFittingCostFunction<ARRAY> BaseClass;
        typedef typename BaseClass::REAL REAL;

        weightedLeastSquareErrorCostFunction(const ARRAY& w) { w_ = w; }

        virtual REAL eval(const ARRAY &y, const ARRAY& y_est)
        {
            GADGET_CHECK_THROW(y.size()==y_est.size());

            REAL cost(0);

            size_t num = y.size();
            for (size_t ii=0; ii<num; ii++)
            {
                REAL v = std::abs(y_est[ii] - y[ii]);
                cost += std::abs(w_[ii]*w_[ii])*v*v;
            }

            if(num>1) cost /= num;

            return cost;
        }

    protected:
        ARRAY w_;
    };
}
