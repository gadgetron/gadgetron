/** \file       gtPlusLSQRSolver.h
    \brief      Implement the LSQR linear solver for Ax=b
    \author     Hui Xue

    Ref to:
    http://www.stanford.edu/group/SOL/software/lsqr.html
    C. C. Paige and M. A. Saunders, LSQR: An algorithm for sparse linear equations and sparse least squares, TOMS 8(1), 43-71 (1982). 
    C. C. Paige and M. A. Saunders, Algorithm 583; LSQR: Sparse linear equations and least-squares problems, TOMS 8(2), 195-209 (1982).
*/

#pragma once

#include "gtPlusLinearSolver.h"

namespace Gadgetron { namespace gtPlus {

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
class gtPlusLSQRSolver : public gtPlusLinearSolver<Array_Type_I, Array_Type_O, Oper_Type>
{
public:

    typedef gtPlusLinearSolver<Array_Type_I, Array_Type_O, Oper_Type> BaseClass;

    typedef typename BaseClass::ValueType ValueType;

    typedef typename realType<ValueType>::Type value_type;

    gtPlusLSQRSolver();
    virtual ~gtPlusLSQRSolver();

    virtual bool solve(const Array_Type_I& b, Array_Type_O& x);

    virtual void printInfo(std::ostream& os) const;

    using BaseClass::iterMax_;
    using BaseClass::thres_;
    using BaseClass::x0_;
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
    using BaseClass::oper_;
};

// ===================================================================================== //
//                           Implementation of template function                         //
// ===================================================================================== //

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
gtPlusLSQRSolver<Array_Type_I, Array_Type_O, Oper_Type>::
gtPlusLSQRSolver() : BaseClass()
{
    iterMax_ = 70;
    thres_ = (value_type)1e-4;
}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
gtPlusLSQRSolver<Array_Type_I, Array_Type_O, Oper_Type>::
~gtPlusLSQRSolver() 
{

}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
bool gtPlusLSQRSolver<Array_Type_I, Array_Type_O, Oper_Type>::
solve(const Array_Type_I& b, Array_Type_O& x)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(oper_!=NULL);

        x = *x0_;

        // Set up for the method
        value_type n2b;
        Gadgetron::norm2(b, n2b);

        int flag = 1;

        value_type tolb = thres_ * n2b;
        Array_Type_I u(b);

        // u = u - A(x, varargin{:}, 'notransp');
        // u = b - A*x0
        GADGET_CHECK_RETURN_FALSE(oper_->forwardOperator(x, u));
        GADGET_CHECK_RETURN_FALSE(Gadgetron::subtract(b, u, u));

        value_type beta;
        Gadgetron::norm2(u, beta);

        value_type normr(beta);
        if (std::abs(beta)>0)
        {
            GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( value_type(1.0)/beta, u));
        }

        value_type c = 1;
        value_type s = 0;
        value_type phibar = beta;

        // v = A(u, varargin{:},'transp');
        Array_Type_I v(x);
        GADGET_CHECK_RETURN_FALSE(oper_->adjointOperator(u, v));

        value_type alpha;
        Gadgetron::norm2(v, alpha);
        if (std::abs(alpha)>0)
        {
            GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( value_type(1.0)/alpha, v));
        }

        Array_Type_I d(x);
        Gadgetron::clear(d);

        value_type normar;
        normar = alpha * beta;

        // Check for all zero solution
        if ( std::abs(normar) < DBL_EPSILON )
        {
            Gadgetron::clear(x);
            return true;
        }

        value_type norma(0);
        value_type sumnormd2 = 0;
        size_t stag = 0;
        size_t iter = iterMax_;
        size_t  maxstagsteps = 3;

        // loop over maxit iterations (unless convergence or failure)

        Array_Type_I z(v), dtmp(d), ztmp(v), vt(v), utmp(u);
        Array_Type_I normaVec(3);

        value_type thet, rhot, rho, phi, tmp, tmp2;

        size_t ii;
        for ( ii=0; ii<iterMax_; ii++ )
        {
            // z = v;
            memcpy(z.begin(), v.begin(), v.get_number_of_bytes());

            // u = A(z, varargin{:},'notransp') - alpha*u;
            GADGET_CHECK_RETURN_FALSE(oper_->forwardOperator(z, utmp));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( alpha, u));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::subtract( utmp, u, u));

            Gadgetron::norm2(u, beta);
            GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( value_type(1.0)/beta, u));

            normaVec(0) = norma;
            normaVec(1) = alpha;
            normaVec(2) = beta;
            Gadgetron::norm2(normaVec, norma);

            thet = - s * alpha;
            rhot = c * alpha;
            rho = (value_type)( std::sqrt( (double)(rhot*rhot + beta*beta) ));
            c = rhot / rho;
            s = - beta / rho;
            phi = c * phibar;
            if ( std::abs(phi)< DBL_EPSILON )
            {
                stag = 1;
            }

            phibar = s * phibar;

            // d = (z - thet * d) / rho;
            //dtmp = d;
            memcpy(dtmp.begin(), d.begin(), d.get_number_of_bytes());
            GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( thet, dtmp));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::subtract( z, dtmp, ztmp));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( value_type(1.0)/rho, ztmp));
            //d = ztmp;
            memcpy(d.begin(), ztmp.begin(), d.get_number_of_bytes());

            // sumnormd2 = sumnormd2 + (norm(d(:)))^2;
            Gadgetron::norm2(d, tmp);
            sumnormd2 += (tmp*tmp);

            // Check for stagnation of the method
            Gadgetron::norm2(x, tmp2);

            if ( std::abs(phi)*std::abs(tmp) < DBL_EPSILON*std::abs(tmp2) )
            {
                stag++;
            }
            else
            {
                stag = 0;
            }

            // check for convergence in min{|b-A*x|}
            if ( std::abs(normar/(norma*normr)) <= thres_ )
            {
                flag = 0;
                break;
            }

            // check for convergence in A*x=b
            if (std::abs(normr) <= std::abs(tolb) )
            {
                flag = 0;
                break;
            }

            if (stag >= maxstagsteps)
            {
                flag = 3;
                break;
            }

            //if (printIter_)
            //{
            //    GADGET_MSG("Iteration " << ii << " - normar/(norma*normr) = " << std::abs(normar/(norma*normr)) << " - normr = " << std::abs(normr) );
            //}

            // x = x + phi * d;
            //dtmp = d;
            memcpy(dtmp.begin(), d.begin(), d.get_number_of_bytes());
            GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( phi, dtmp));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::add( x, dtmp, x));

            normr = (value_type)(std::abs( (double)s) * normr);

            // vt = A(u, varargin{:},'transp');
            GADGET_CHECK_RETURN_FALSE(oper_->adjointOperator(u, vt));

            // v = vt - beta * v;
            GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( beta, v));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::subtract( vt, v, v));

            Gadgetron::norm2(v, alpha);

            GADGET_CHECK_RETURN_FALSE(Gadgetron::scal( value_type(1.0)/alpha, v));

            normar = alpha * std::abs( (value_type)s * phi);
        }

        if (printIter_)
        {
            GADGET_MSG("Total iteration number is  " << ii << " - relative norm is " << std::abs(normar/(norma*normr)) << " ... ");
        }

        if (flag == 1)
        {
            if ( normar/(norma*normr) <= thres_ )
            {
                flag = 0;
            }

            if (std::abs(normr) <= std::abs(tolb) )
            {
                flag = 0;
            }
        }

        //if (printIter_)
        //{
        //    value_type relres = normr/n2b;
        //    GADGET_MSG("Flag = " << flag << " - relres = " << std::abs(relres) );
        //}
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in gtPlusLSQRSolver<Array_Type_I, Array_Type_O, Oper_Type>::solve(...) ... ");
        return false;
    }

    return true;
}

template <typename Array_Type_I, typename Array_Type_O, typename Oper_Type>
void gtPlusLSQRSolver<Array_Type_I, Array_Type_O, Oper_Type>::
printInfo(std::ostream& os) const
{
    os << "-------------- GTPlus ISMRMRD linear LSQR solver -------------" << std::endl;
    os << "The linear solver solves Ax=b problem" << std::endl;
    os << "------------------------------------------------------------" << std::endl;
}

}}
