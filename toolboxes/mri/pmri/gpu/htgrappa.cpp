#include "htgrappa.h"
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_linalg.h"

#ifdef USE_OMP
#include <omp.h>
#endif


/*
  This file is used to hide certain Armadillo calls from the nvcc compiler. If Armadillo functions need to
  be called in a *.cu file, it is preferably to wrap the calls in a function and place that function in
  a *.cpp file so that Armadillo code will not be compiled by nvcc.

  Some error handling may be needed in these functions, but eventually SymmetricHermitianPositiveDefiniteLinearSystem_posv
  will be renamed and made to throw exceptions and then it should be handled.

 */



namespace Gadgetron
{
  template <class T> void ht_grappa_solve_spd_system(hoNDArray<T> *A, hoNDArray<T> *B) {
    /*
      We are swithcing off OpenMP threading before this call to posv. There seems to be a bad interaction between openmp, cuda, and BLAS.
      So far this problem has only been observed from *.cu files (or in functions called from *.cu files) but the problem may be more general.

      This is a temporary fix that we should keep an eye on.
     */

    hoNDArray<T> A_ori;
    A_ori = *A;

#ifdef USE_OMP
    int num_threads = omp_get_num_threads();
    omp_set_num_threads(1);
#endif //USE_OMP

    try
    {
        posv(*A, *B);
    }
    catch(...)
    {
        // it is found that if signal is very very high, the posv can throw exceptions due to ill-conditioned matrix of A
        // hesv does not require A to be a positive-definite matrix, but an n-by-n symmetric matrix

        GERROR_STREAM("ht_grappa_solve_spd_system : posv(*A, *B) throws exceptions ... ");
        *A = A_ori;
        hesv(*A, *B);
        GERROR_STREAM("ht_grappa_solve_spd_system : hesv(*A, *B) is called ");
    }

#ifdef USE_OMP
    omp_set_num_threads(num_threads);
#endif //USE_OMP

  }

  template void ht_grappa_solve_spd_system< float_complext >(hoNDArray< float_complext > *A, hoNDArray< float_complext > *B);

}
