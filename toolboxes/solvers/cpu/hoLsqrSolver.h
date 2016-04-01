/** \file   hoLsqrSolver.h
    \brief  Implementation a cpu version LSQR solver
    \author Hui Xue

    Ref to:
    http://www.stanford.edu/group/SOL/software/lsqr.html
    C. C. Paige and M. A. Saunders, LSQR: An algorithm for sparse linear equations and sparse least squares, TOMS 8(1), 43-71 (1982).
    C. C. Paige and M. A. Saunders, Algorithm 583; LSQR: Sparse linear equations and least-squares problems, TOMS 8(2), 195-209 (1982).
*/

#pragma once

#include "cpusolver_export.h"
#include "hoNDArray_math.h"
#include "lsqrSolver.h"

namespace Gadgetron{

    template <class T> class EXPORTCPUSOLVER hoLsqrSolver : public lsqrSolver< hoNDArray<T> >
    {
    public:

        typedef lsqrSolver< hoNDArray<T> > BaseClass;

        hoLsqrSolver();
        virtual ~hoLsqrSolver();

    protected:

        using BaseClass::encoding_operator_;
        using BaseClass::iterations_;
        using BaseClass::tc_tolerance_;
    };
}
