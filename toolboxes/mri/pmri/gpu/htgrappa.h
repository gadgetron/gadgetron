/** \file htgrappa.h
    \brief Utilities to calibrate grappa weights and corresponding unmixing coefficients - GPU-based.
*/

#pragma once
#ifndef HTGRAPPA_H
#define HTGRAPPA_H

#include "cuNDArray.h"

#include <list>

namespace Gadgetron
{

  template <class T>
  int htgrappa_calculate_grappa_unmixing(cuNDArray<T>* ref_data,
                                         cuNDArray<T>* b1,
                                         unsigned int acceleration_factor,
                                         std::vector<unsigned int>* kernel_size,
                                         cuNDArray<T>* out_mixing_coeff,
                                         std::vector< std::pair<unsigned int, unsigned int> >* sampled_region = 0,
                                         std::list< unsigned int >* uncombined_channels = 0);

  template <class T>
  int inverse_clib_matrix(cuNDArray<T>* A,
                          cuNDArray<T>* b,
                          cuNDArray<T>* out_mixing_coeff,
                          double lamda);

  template <class T> void ht_grappa_solve_spd_system(hoNDArray<T> *A, hoNDArray<T> *B);

    template <class T> void ht_grappa_solve_spd_system(cuNDArray<T>& A, cuNDArray<T>& B);
}

#endif //HTGRAPPA_H
