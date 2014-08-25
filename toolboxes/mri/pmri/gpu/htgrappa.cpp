#include "htgrappa.h"
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"
#include "hoMatrix_util.h"

/*
  This file is used to hide certain Armadillo calls from the nvcc compiler. If Armadillo functions need to 
  be called in a *.cu file, it is preferably to wrap the calls in a function and place that function in 
  a *.cpp file so that Armadillo code will not be compiled by nvcc. 
 */

namespace Gadgetron
{
  template <class T> void ht_grappa_solve_spd_system(hoNDArray<T> *A, hoNDArray<T> *B) {
    hoMatrix< std::complex<float> > Am(A->get_size(0),A->get_size(1),(std::complex<float>*)A->get_data_ptr(),false);
    hoMatrix< std::complex<float> > Bm(B->get_size(0),B->get_size(1),(std::complex<float>*)B->get_data_ptr(),false);
    SymmetricHermitianPositiveDefiniteLinearSystem_posv(Am, Bm);
  }

  template void ht_grappa_solve_spd_system< float_complext >(hoNDArray< float_complext > *A, hoNDArray< float_complext > *B);
  
}
