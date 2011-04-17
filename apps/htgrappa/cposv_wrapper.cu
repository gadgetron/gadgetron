#include "cposv_wrapper.h"

#include "magma.h"

int cposv_wrapper(cuNDArray<float2>* A, cuNDArray<float2>* rhs)
{
  
  magma_int_t mret, m, nrhs, info;
  
  if (A->get_number_of_dimensions() > 2) {
    std::cerr << "cgels_wrapper: matrix A must be a 2 dimensional matrix" << std::endl;
    return -1;
  }

  if (rhs->get_number_of_dimensions() > 2) {
    std::cerr << "cgels_wrapper: matrix rhs must be a 1 or 2 dimensional matrix" << std::endl;
    return -1;
  }

  m = A->get_size(0);

  magma_int_t rows_rhs;
  if (rhs->get_number_of_dimensions() == 1) {
    nrhs = 1;
    rows_rhs = rhs->get_number_of_elements();
  } else {
    nrhs = rhs->get_size(1);
    rows_rhs = rhs->get_size(0);
  }

  if (rows_rhs != m) {
    std::cerr << "cgels_wrapper: Matrix dimensions mismatch" << std::endl;
    return -1;
  }

  info = 0;
  
  mret = magma_cposv_gpu( 'U', m, nrhs, 
			  A->get_data_ptr(), m, 
			  rhs->get_data_ptr(), m, &info );

  if (mret) {
    std::cout << "htgrappa_calculate_grappa_unmixing: Cholesky failed" << std::endl;
    return -1;
  }

  return 0;
}
