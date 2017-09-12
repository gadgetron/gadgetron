/** \file twoParaExpDecayOperator.h
    \brief Implement two parameter exponential decay operator.
    \author     Hui Xue
*/

#pragma once

#include "curveFittingOperator.h"
#include <boost/math/special_functions/sign.hpp>

namespace Gadgetron {

    // y = bi[0] * exp(-x/bi[1])

    template <class ARRAY> class twoParaExpDecayOperator : public curveFittingOperator<ARRAY>
    {
    public:

        typedef curveFittingOperator<ARRAY> BaseClass;
        typedef typename BaseClass::ELEMENT_TYPE ELEMENT_TYPE;
        typedef typename BaseClass::REAL REAL;

        twoParaExpDecayOperator();
        virtual ~twoParaExpDecayOperator();

        virtual void gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad);
        virtual void magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y);

    protected:
    };

    template <class ARRAY>
    twoParaExpDecayOperator<ARRAY>::twoParaExpDecayOperator() : BaseClass()
    {
    }

    template <class ARRAY>
    twoParaExpDecayOperator<ARRAY>::~twoParaExpDecayOperator()
    {
    }

    template <class ARRAY>
    void twoParaExpDecayOperator<ARRAY>::gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad)
    {
        try
        {
            size_t num = b.size();
            if(grad.size()!=num) grad.resize(num, 0);

            int sign_b1 = boost::math::sign(b[1]);
            ELEMENT_TYPE rb = 1.0 / ( (std::abs(b[1])<FLT_EPSILON) ? sign_b1*FLT_EPSILON : b[1] );

            ELEMENT_TYPE val = std::exp(-1 * xi * rb);
            grad[0] = val;
            grad[1] = b[0] * val * xi * rb * rb;
        }
        catch(...)
        {
            GADGET_THROW("Errors happened in twoParaExpDecayOperator<ARRAY>::gradient(...) ... ");
        }
    }

    template <class ARRAY>
    void twoParaExpDecayOperator<ARRAY>::magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y)
    {
        size_t num = x.size();
        if(y.size()!=x.size()) y.resize(num, 0);

        int sign_b1 = boost::math::sign(b[1]);
        ELEMENT_TYPE rb = 1.0 / ( (std::abs(b[1])<FLT_EPSILON) ? sign_b1*FLT_EPSILON : b[1] );

        size_t ii;
        for (ii=0; ii<num; ii++)
        {
            y[ii] = b[0] * exp( -1 * x[ii] * rb);
        }
    }
}
