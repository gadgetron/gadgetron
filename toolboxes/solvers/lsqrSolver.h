#pragma once

#include "linearOperator.h"
#include "linearOperatorSolver.h"
#include "cgPreconditioner.h"
#include "real_utilities.h"

#include <vector>
#include <iostream>
#include "encodingOperatorContainer.h"

namespace Gadgetron {

template <class ARRAY_TYPE> class lsqrSolver: public linearOperatorSolver<ARRAY_TYPE>
{
protected:

    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;

public:

    lsqrSolver()
    {
        iterations_ = 10;
        tc_tolerance_ = (REAL)1e-3;
    }

    virtual ~lsqrSolver() {}

    virtual void set_tc_tolerance( REAL tolerance ) { tc_tolerance_ = tolerance; }
    virtual REAL get_tc_tolerance() { return tc_tolerance_; }

    virtual void set_max_iterations( unsigned int iterations ) { iterations_ = iterations; }
    virtual unsigned int get_max_iterations() { return iterations_; }

    virtual void solve(ARRAY_TYPE* x, ARRAY_TYPE* b)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > image_dims = this->encoding_operator_->get_domain_dimensions();

            GADGET_CHECK_THROW(x != NULL);
            GADGET_CHECK_THROW(b != NULL);

            GADGET_CHECK_THROW(b->dimensions_equal(image_dims.get()));

            if (this->x0_ != NULL)
            {
                GADGET_CHECK_THROW(this->x0_->dimensions_equal(image_dims.get()));
                *x = *(this->x0_);
            }
            else
            {
                x->create(*image_dims);
                Gadgetron::clear(*x);
            }

            REAL n2b = Gadgetron::nrm2(b);

            int flag = 1;

            REAL tolb = tc_tolerance_ * n2b;
            ARRAY_TYPE u(*b);

            this->encoding_operator_->mult_M(x, &u);
            Gadgetron::subtract(*b, u, u);

            REAL beta = Gadgetron::nrm2(&u);

            REAL normr(beta);
            if (std::abs(beta)>0)
            {
                Gadgetron::scal(REAL(1.0) / beta, u);
            }

            REAL c = 1;
            REAL s = 0;
            REAL phibar = beta;

            ARRAY_TYPE v(*x);
            this->encoding_operator_->mult_MH(&u, &v);

            REAL alpha = Gadgetron::nrm2(&v);
            if (std::abs(alpha)>0)
            {
                Gadgetron::scal(REAL(1.0) / alpha, v);
            }

            ARRAY_TYPE d(*x);
            Gadgetron::clear(d);

            REAL normar;
            normar = alpha * beta;

            // Check for all zero solution
            if (std::abs(normar) < DBL_EPSILON)
            {
                Gadgetron::clear(x);
                return;
            }

            REAL norma(0);
            REAL sumnormd2 = 0;
            size_t stag = 0;
            size_t iter = iterations_;
            size_t  maxstagsteps = 3;

            ARRAY_TYPE z(v), dtmp(d), ztmp(v), vt(v), utmp(u);
            ARRAY_TYPE normaVec(3);

            REAL thet, rhot, rho, phi, tmp, tmp2;

            size_t ii;
            for (ii = 0; ii<iterations_; ii++)
            {
                z = v;

                this->encoding_operator_->mult_M(&z, &utmp);
                Gadgetron::scal(alpha, u);
                Gadgetron::subtract(utmp, u, u);

                beta = Gadgetron::nrm2(&u);
                Gadgetron::scal(REAL(1.0) / beta, u);

                normaVec(0) = norma;
                normaVec(1) = alpha;
                normaVec(2) = beta;
                norma = Gadgetron::nrm2(&normaVec);

                thet = -s * alpha;
                rhot = c * alpha;
                rho = (REAL)(std::sqrt((double)(rhot*rhot + beta*beta)));
                c = rhot / rho;
                s = -beta / rho;
                phi = c * phibar;
                if (std::abs(phi)< DBL_EPSILON)
                {
                    stag = 1;
                }

                phibar = s * phibar;

                dtmp = d;
                Gadgetron::scal(thet, dtmp);
                Gadgetron::subtract(z, dtmp, ztmp);
                Gadgetron::scal(REAL(1.0) / rho, ztmp);

                d = ztmp;
                tmp = Gadgetron::nrm2(&d);
                sumnormd2 += (tmp*tmp);

                // Check for stagnation of the method
                tmp2 = Gadgetron::nrm2(x);

                if (std::abs(phi)*std::abs(tmp) < DBL_EPSILON*std::abs(tmp2))
                {
                    stag++;
                }
                else
                {
                    stag = 0;
                }

                // check for convergence in min{|b-A*x|}
                if (std::abs(normar / (norma*normr)) <= tc_tolerance_)
                {
                    flag = 0;
                    break;
                }

                // check for convergence in A*x=b
                if (std::abs(normr) <= std::abs(tolb))
                {
                    flag = 0;
                    break;
                }

                if (stag >= maxstagsteps)
                {
                    flag = 3;
                    break;
                }

                // memcpy(dtmp.begin(), d.begin(), d.get_number_of_bytes());
                dtmp = d;

                Gadgetron::scal(phi, dtmp);
                Gadgetron::add(*x, dtmp, *x);

                normr = (REAL)(std::abs((double)s) * normr);

                this->encoding_operator_->mult_MH(&u, &vt);

                Gadgetron::scal(beta, v);
                Gadgetron::subtract(vt, v, v);

                alpha = Gadgetron::nrm2(&v);

                Gadgetron::scal(REAL(1.0) / alpha, v);

                normar = alpha * std::abs((REAL)s * phi);
            }

            if (this->output_mode_ >= solver<ARRAY_TYPE, ARRAY_TYPE>::OUTPUT_VERBOSE)
            {
                GDEBUG_STREAM("Total iteration number is  " << ii << " - relative norm is " << std::abs(normar / (norma*normr)) << " ... ");
            }

            if (flag == 1)
            {
                if (normar / (norma*normr) <= tc_tolerance_)
                {
                    flag = 0;
                }

                if (std::abs(normr) <= std::abs(tolb))
                {
                    flag = 0;
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in lsqrSolver<ARRAY_TYPE>::solve(x, b) ... ");
        }
    }

    virtual boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE *b)
    {
        boost::shared_ptr<ARRAY_TYPE> x = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE());
        this->solve(x.get(), b);
        return x;
    }

protected:

    unsigned int iterations_;
    REAL tc_tolerance_;
};

}
