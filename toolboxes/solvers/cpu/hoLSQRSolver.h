/** \file   hoLSQRSolver.h
    \brief  Implementation a cpu version LSQR solver
    \author Hui Xue

    Ref to:
    http://www.stanford.edu/group/SOL/software/lsqr.html
    C. C. Paige and M. A. Saunders, LSQR: An algorithm for sparse linear equations and sparse least squares, TOMS 8(1), 43-71 (1982).
    C. C. Paige and M. A. Saunders, Algorithm 583; LSQR: Sparse linear equations and least-squares problems, TOMS 8(2), 195-209 (1982).
*/

#pragma once

#include "cpusolver_export.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "lsqrSolver.h"

namespace Gadgetron{

    template <class ARRAY_TYPE> class EXPORTCPUSOLVER hoLSQRSolver : public lsqrSolver<ARRAY_TYPE>
    {
    public:

        typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
        typedef typename realType<ELEMENT_TYPE>::Type REAL;

        typedef lsqrSolver<ARRAY_TYPE> BaseClass;

        hoLSQRSolver();
        virtual ~hoLSQRSolver();

        void set_verbose(bool v) { verbose_ = v; }
        bool get_verbose() { return verbose_; }

        // solve the linear equation and results are stored in x
        virtual void solve(ARRAY_TYPE* x, ARRAY_TYPE* b);

        // return the results as a shared array
        virtual boost::shared_ptr<ARRAY_TYPE> solve(ARRAY_TYPE *b);

    protected:

        using BaseClass::encoding_operator_;
        using BaseClass::iterations_;
        using BaseClass::tc_tolerance_;

        bool verbose_;
    };
}
