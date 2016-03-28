/** \file   hoLSQRSolver.h
    \brief  Implementation a cpu version LSQR solver
    \author Hui Xue

    Ref to:
    http://www.stanford.edu/group/SOL/software/lsqr.html
    C. C. Paige and M. A. Saunders, LSQR: An algorithm for sparse linear equations and sparse least squares, TOMS 8(1), 43-71 (1982).
    C. C. Paige and M. A. Saunders, Algorithm 583; LSQR: Sparse linear equations and least-squares problems, TOMS 8(2), 195-209 (1982).
*/

#pragma once

#include "linearOperatorSolver.h"
#include "cpusolver_export.h"

namespace Gadgetron{

    template <class ARRAY_TYPE> class EXPORTCPUSOLVER hoLSQRSolver : public linearOperatorSolver<ARRAY_TYPE>
    {
    public:

        typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
        typedef typename realType<ELEMENT_TYPE>::Type REAL;

        typedef linearOperatorSolver<ARRAY_TYPE> BaseClass;

        hoLSQRSolver();
        virtual ~hoLSQRSolver();

        virtual void set_tc_tolerance(REAL tolerance) { tc_tolerance_ = tolerance; }
        virtual REAL get_tc_tolerance() { return tc_tolerance_; }

        virtual void set_max_iterations(unsigned int iterations) { iterations_ = iterations; }
        virtual unsigned int get_max_iterations() { return iterations_; }

        void set_x0(ARRAY_TYPE* x0) { x0_ = x0; }

        void set_verbose(bool v) { verbose_ = v; }
        bool get_verbose() { return verbose_; }

        // solve the linear equation and results are stored in x
        virtual void solve(ARRAY_TYPE* x, ARRAY_TYPE* b);

        // return the results as a shared array
        virtual boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE *b);

    protected:

        using BaseClass::encoding_operator_;

        // initial guess for x
        // if not set, the all zeros are used
        ARRAY_TYPE* x0_;

        bool verbose_;

        unsigned int iterations_;
        REAL tc_tolerance_;
    };
}
