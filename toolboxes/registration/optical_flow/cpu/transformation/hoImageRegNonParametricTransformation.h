/** \file   hoImageRegNonParametricTransformation.h
    \brief  Define the base class for the non-parametric geometry transformation in gadgetron registration
    \author Hui Xue
*/

#ifndef hoImageRegNonParametricTransformation_H_
#define hoImageRegNonParametricTransformation_H_

#pragma once

#include "hoImageRegTransformation.h"

namespace Gadgetron {

    /// non-parametric transformation, e.g. deformation field
    template<typename ValueType, unsigned int DIn, unsigned int DOut> 
    class hoImageRegNonParametricTransformation : public hoImageRegTransformation<ValueType, DIn, DOut>
    {
    public:

        typedef hoImageRegTransformation<ValueType, DIn, DOut> BaseClass;
        typedef hoImageRegNonParametricTransformation<ValueType, DIn, DOut> Self;

        typedef ValueType T;

        typedef typename BaseClass::input_point_type input_point_type;
        typedef typename BaseClass::output_point_type output_point_type;

        typedef typename BaseClass::jacobian_position_type jacobian_position_type;

        hoImageRegNonParametricTransformation() : BaseClass() {}
        virtual ~hoImageRegNonParametricTransformation() {}

        virtual bool invertTransformation() = 0;

        virtual bool setIdentity() = 0;

        virtual bool transform(const T* pt_in, T* pt_out) const = 0;

        virtual bool transform(const T& xi, const T& yi, T& xo, T& yo) const = 0;

        virtual bool transform(const T& xi, const T& yi, const T& zi, T& xo, T& yo, T& zo) const = 0;

        virtual bool transform(const size_t* pt_in, T* pt_out) const = 0;
        virtual bool transform(const size_t* pt_in, size_t N, T* pt_out) const = 0;
        virtual bool transform(const size_t& xi, const size_t& yi, T& xo, T& yo) const = 0;
        virtual bool transform(const size_t* xi, const size_t* yi, size_t N, T* xo, T* yo) const = 0;
        virtual bool transform(const size_t& xi, const size_t& yi, const size_t& zi, T& xo, T& yo, T& zo) const = 0;
        virtual bool transform(const size_t* xi, const size_t* yi, const size_t* zi, size_t N, T* xo, T* yo, T* zo) const = 0;

        /// compute jacobian matrix to spatial position
        /// DOut*DIn matrix
        virtual bool jacobianPosition(const input_point_type& /*pos*/, jacobian_position_type& jac)
        {
            jac.createMatrix(DOut, DIn);
            jac.setIdentity();
            return true;
        }

        virtual void print(std::ostream& os) const
        {
            using namespace std;
            os << "--------------Gagdgetron non-parametric geometry transformation -------------" << endl;
            os << "Input dimension is : " << DIn << endl;
            os << "Output dimension is : " << DOut << endl;

            std::string elemTypeName = std::string(typeid(T).name());
            os << "Transformation data type is : " << elemTypeName << std::endl;
        }

        virtual std::string transformationName() const
        {
            return std::string("hoImageRegNonParametricTransformation"); 
        }

        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

    protected:
    };
}
#endif // hoImageRegNonParametricTransformation_H_
