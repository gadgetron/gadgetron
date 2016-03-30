
#include "hoLSQRSolver.h"

namespace Gadgetron{

    template <typename ARRAY_TYPE>
    hoLSQRSolver<ARRAY_TYPE>::hoLSQRSolver() : BaseClass()
    {
    }

    template <typename ARRAY_TYPE>
    hoLSQRSolver<ARRAY_TYPE>::~hoLSQRSolver()
    {
    }

    template class EXPORTCPUSOLVER hoLSQRSolver< std::complex<float> >;
    template class EXPORTCPUSOLVER hoLSQRSolver< std::complex<double> >;
}
