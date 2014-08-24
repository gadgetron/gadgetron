/** \file       gtPlusNonLinearSolver.h
    \brief      Define the base class for GtPlus non-linear solvers
    \author     Hui Xue
*/

#pragma once

#include "gtPlusSolver.h"

namespace Gadgetron { namespace gtPlus {

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
class gtPlusNonLinearSolver : public gtPlusSolver<Array_Type_I, Array_Type_O>
{
public:

    typedef gtPlusSolver<Array_Type_I, Array_Type_O> BaseClass;

    typedef typename BaseClass::ValueType ValueType;
    typedef typename realType<ValueType>::Type value_type;

    // one operator is related to a weight
    typedef std::pair<Oper_Type*, ValueType> Oper_Elem_Type;
    // multiple operator can be added to a solver
    typedef std::vector< Oper_Elem_Type > Oper_List_Type;

    gtPlusNonLinearSolver();
    virtual ~gtPlusNonLinearSolver();

    Oper_List_Type getOperList();
    void setOperList(Oper_List_Type& opero);

    void add(Oper_Type& op, ValueType a);
    void remove(Oper_Type*& op, ValueType& a);

    bool set(size_t ind, Oper_Type& op, ValueType a);

    // main function to perform the solver
    virtual bool solve(const Array_Type_I& b, Array_Type_O& x) = 0;

    virtual void printInfo(std::ostream& os) const;

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::printIter_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_complex_;
    using BaseClass::gtPlus_mem_manager_;

protected:

    using BaseClass::callback_;
    Oper_List_Type operList_;
};

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
gtPlusNonLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::
gtPlusNonLinearSolver() : BaseClass()
{

}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
gtPlusNonLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::
~gtPlusNonLinearSolver() 
{
    operList_.clear();
}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
void gtPlusNonLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::
add(Oper_Type& op, ValueType a)
{
    operList_.push_back(Oper_Elem_Type(&op, a));
}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
void gtPlusNonLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::
remove(Oper_Type*& op, ValueType& a)
{
    if ( operList_.empty() )
    {
        op = NULL;
        a = 0;
    }
    else
    {
        op = operList_[operList_.size()-1].first;
        a = operList_[operList_.size()-1].second;
        operList_.pop_back();
    }
}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
bool gtPlusNonLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::
set(size_t ind, Oper_Type& op, ValueType a)
{
    if ( ind >= operList_.size() )
    {
        GADGET_WARN_MSG("ind >= operList_.size()");
    }

    operList_[ind].first = &op;
    operList_[ind].second = a;

    return true;
}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
void gtPlusNonLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>::
printInfo(std::ostream& os) const
{
    using namespace std;
    os << "-------------- GTPlus ISMRMRD nonlinear solver -----------------" << endl;
    os << "ISMRMRD nonlinear solver ... " << endl;
    os << "----------------------------------------------------------------" << endl;
}

}}
