
#include "hoLSQRSolver.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"

namespace Gadgetron{

    template <typename ARRAY_TYPE>
    hoLSQRSolver<ARRAY_TYPE>::hoLSQRSolver() : BaseClass(), x0_(NULL), verbose_(false), iterations_(70), tc_tolerance_( (REAL)1e-4 )
    {
    }

    template <typename ARRAY_TYPE>
    hoLSQRSolver<ARRAY_TYPE>::~hoLSQRSolver()
    {
    }

    template <typename ARRAY_TYPE>
    void hoLSQRSolver<ARRAY_TYPE>::solve(ARRAY_TYPE* x, ARRAY_TYPE* b)
    {
        try
        {
            boost::shared_ptr< std::vector<size_t> > image_dims = this->encoding_operator_->get_domain_dimensions();

            GADGET_CHECK_THROW(x!=NULL);
            GADGET_CHECK_THROW(b != NULL);

            GADGET_CHECK_THROW(b->dimensions_equal(image_dims.get()));

            if (x0_ != NULL)
            {
                GADGET_CHECK_THROW(x0_->dimensions_equal(image_dims.get()));
                *x = *x0_;
            }
            else
            {
                x->create(*image_dims);
                Gadgetron::clear(*x);
            }

            REAL n2b;
            Gadgetron::norm2(*b, n2b);

            int flag = 1;

            REAL tolb = tc_tolerance_ * n2b;
            ARRAY_TYPE u(*b);

            encoding_operator_->mult_M(x, &u);
            Gadgetron::subtract(*b, u, u);

            REAL beta;
            Gadgetron::norm2(u, beta);

            REAL normr(beta);
            if (std::abs(beta)>0)
            {
                Gadgetron::scal(REAL(1.0) / beta, u);
            }

            REAL c = 1;
            REAL s = 0;
            REAL phibar = beta;

            ARRAY_TYPE v(x);
            encoding_operator_->mult_MH(&u, &v);

            REAL alpha;
            Gadgetron::norm2(v, alpha);
            if (std::abs(alpha)>0)
            {
                Gadgetron::scal(REAL(1.0) / alpha, v);
            }

            ARRAY_TYPE d(x);
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
                memcpy(z.begin(), v.begin(), v.get_number_of_bytes());

                encoding_operator_->mult_M(&z, &utmp);
                Gadgetron::scal(alpha, u);
                Gadgetron::subtract(utmp, u, u);

                Gadgetron::norm2(u, beta);
                Gadgetron::scal(REAL(1.0) / beta, u);

                normaVec(0) = norma;
                normaVec(1) = alpha;
                normaVec(2) = beta;
                Gadgetron::norm2(normaVec, norma);

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

                memcpy(dtmp.begin(), d.begin(), d.get_number_of_bytes());
                Gadgetron::scal(thet, dtmp);
                Gadgetron::subtract(z, dtmp, ztmp);
                Gadgetron::scal(REAL(1.0) / rho, ztmp);

                memcpy(d.begin(), ztmp.begin(), d.get_number_of_bytes());

                Gadgetron::norm2(d, tmp);
                sumnormd2 += (tmp*tmp);

                // Check for stagnation of the method
                Gadgetron::norm2(*x, tmp2);

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

                memcpy(dtmp.begin(), d.begin(), d.get_number_of_bytes());
                Gadgetron::scal(phi, dtmp);
                Gadgetron::add(*x, dtmp, *x);

                normr = (REAL)(std::abs((double)s) * normr);

                encoding_operator_->mult_MH(&u, &vt);

                Gadgetron::scal(beta, v);
                Gadgetron::subtract(vt, v, v);

                Gadgetron::norm2(v, alpha);

                Gadgetron::scal(REAL(1.0) / alpha, v);

                normar = alpha * std::abs((REAL)s * phi);
            }

            if (verbose_)
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
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoLSQRSolver<ARRAY_TYPE>::solve(x, b) ... ");
        }
    }

    template <typename ARRAY_TYPE>
    boost::shared_ptr<ARRAY_TYPE> hoLSQRSolver<ARRAY_TYPE>::solve(ARRAY_TYPE *b)
    {
        boost::shared_ptr<ARRAY_TYPE> x = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE());
        this->solve(x.get(), b);
        return x;
    }

    template class EXPORTCPUSOLVER hoLSQRSolver< hoNDArray< std::complex<float> > >;
    template class EXPORTCPUSOLVER hoLSQRSolver< hoNDArray< std::complex<double> > >;
}
