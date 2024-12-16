/** \file curveFittingOperator.h
    \brief Base class for all operators for curve fitting purpose.
    \author     Hui Xue
*/

#pragma once

#include "log.h"

namespace Gadgetron {

    template <class ARRAY> class curveFittingOperator
    {
    public:

        typedef curveFittingOperator<ARRAY> Self;
        typedef typename ARRAY::value_type ELEMENT_TYPE;
        typedef typename realType<ELEMENT_TYPE>::Type REAL;

        curveFittingOperator();
        virtual ~curveFittingOperator();

        /// Calculates the gradient of the operator at point xi with parameter bi
        virtual void gradient(const ELEMENT_TYPE& xi, const ARRAY& b, ARRAY& grad) = 0;

        /// evaluate the operator at parameter bi for every xi
        virtual void magnitude(const ARRAY& x, const ARRAY& b, ARRAY& y) = 0;

    protected:
    };

    template <class ARRAY>
    curveFittingOperator<ARRAY>::curveFittingOperator()
    {
    }

    template <class ARRAY>
    curveFittingOperator<ARRAY>::~curveFittingOperator()
    {
    }
}
