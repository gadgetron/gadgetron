
#include "hoLsqrSolver.h"

namespace Gadgetron{

    template <typename ARRAY_TYPE>
    hoLsqrSolver<ARRAY_TYPE>::hoLsqrSolver() : BaseClass()
    {
    }

    template <typename ARRAY_TYPE>
    hoLsqrSolver<ARRAY_TYPE>::~hoLsqrSolver()
    {
    }

    template class EXPORTCPUSOLVER hoLsqrSolver< std::complex<float> >;
    template class EXPORTCPUSOLVER hoLsqrSolver< std::complex<double> >;
}
