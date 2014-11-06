#include "htgrappa.h"
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_linalg.h"

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
    hoMatrix< std::complex<float> > Am(A->get_size(0),A->get_size(1),(std::complex<float>*)A->get_data_ptr(),false);
    hoMatrix< std::complex<float> > Bm(B->get_size(0),B->get_size(1),(std::complex<float>*)B->get_data_ptr(),false);
    posv(Am, Bm);
  }

  template void ht_grappa_solve_spd_system< float_complext >(hoNDArray< float_complext > *A, hoNDArray< float_complext > *B);
  
}
