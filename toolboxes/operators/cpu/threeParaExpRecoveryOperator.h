/** \file threeParaExpRecoveryOperator.h
    \brief Implement three parameter exponential recovery operator. e.g. for T1 recovery with look-locker correction
    \author     Hui Xue
*/

#pragma once

#include "curveFittingOperator.h"
#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron {

    // y = bi[0] -  bi[1]* exp(-x/bi[2])

    template <class ARRAY> class threeParaExpRecoveryOperator : public curveFittingOperator<ARRAY>
    {
    public:

        typedef curveFittingOperator<ARRAY> BaseClass;
        typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
        typedef typename BaseClass::REAL REAL;

        threeParaExpRecoveryOperator();
        virtual ~threeParaExpRecoveryOperator();

        virtual void gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad);
        virtual void magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y);

    protected:
    };

    template <class ARRAY>
    threeParaExpRecoveryOperator<ARRAY>::threeParaExpRecoveryOperator() : BaseClass()
    {
    }

    template <class ARRAY>
    threeParaExpRecoveryOperator<ARRAY>::~threeParaExpRecoveryOperator()
    {
    }

    template <class ARRAY>
    void threeParaExpRecoveryOperator<ARRAY>::gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad)
    {
        try
        {
            size_t num = b.size();
            if(grad.size()!=num) grad.resize(num, 0);

            int sign_b2 = boost::math::sign(b[2]);
            ELEMENT_TYPE rb = 1.0 / ((std::abs(b[2])<FLT_EPSILON) ? sign_b2*FLT_EPSILON : b[2]);

            ELEMENT_TYPE val = std::exp(-1 * xi * rb);
            grad[0] = 1;
            grad[1] = -val;
            grad[2] = -1 * b[1] * val * xi * rb * rb;
        }
        catch(...)
        {
            GADGET_THROW("Errors happened in threeParaExpRecoveryOperator<ARRAY>::gradient(...) ... ");
        }
    }

    template <class ARRAY>
    void threeParaExpRecoveryOperator<ARRAY>::magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y)
    {
        size_t num = x.size();
        if(y.size()!=x.size()) y.resize(num, 0);

        int sign_b2 = boost::math::sign(b[2]);
        ELEMENT_TYPE rb = 1.0 / ( (std::abs(b[2])<FLT_EPSILON) ? sign_b2*FLT_EPSILON : b[2] );

        size_t ii;
        for (ii=0; ii<num; ii++)
        {
            y[ii] = b[0] - b[1] * exp( -1 * x[ii] * rb);
        }
    }
}
