/** \file       gtPlusLinearSolver.h
    \brief      Define the base class for linear solver
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSolver.h"

namespace Gadgetron { namespace gtPlus {

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
class gtPlusLinearSolver : public gtPlusSolver<Array_Type_I, Array_Type_O>
{
public:

    typedef gtPlusSolver<Array_Type_I, Array_Type_O> BaseClass;

    typedef typename BaseClass::ValueType ValueType;
    typedef typename realType<ValueType>::Type value_type;

    gtPlusLinearSolver();
    virtual ~gtPlusLinearSolver();

    Oper_Type* get();
    void set(Oper_Type& op);

    virtual bool solve(const Array_Type_I& b, Array_Type_O& x) = 0;

    virtual void printInfo(std::ostream& os) const;

    // number of max iterations
    size_t iterMax_;

    // threshold for detla change of residual
    value_type thres_;

    // initial guess for the solver
    Array_Type_O* x0_;

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::printIter_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;

protected:

    using BaseClass::callback_;
    Oper_Type* oper_;
};

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
gtPlusLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::
gtPlusLinearSolver() : oper_(NULL)
{

}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
gtPlusLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::
~gtPlusLinearSolver() 
{

}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
Oper_Type* gtPlusLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::get()
{
    return oper_;
}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
void gtPlusLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::set(Oper_Type& op)
{
    oper_ = &op;
}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
void gtPlusLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::
printInfo(std::ostream& os) const
{
    using namespace std;
    os << "-------------- GTPlus ISMRMRD linear solver --------------------" << endl;
    os << "Linear solver for GtPlus ISMRMRD package ... " << endl;
    os << "----------------------------------------------------------------" << endl;
}

}}
