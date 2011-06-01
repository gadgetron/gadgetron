#ifndef HTGRAPPA_H
#define HTGRAPPA_H

#include "cuNDArray.h"

template <class T> EXPORTGPUPMRI int htgrappa_calculate_grappa_unmixing(cuNDArray<T>* ref_data, 
							  cuNDArray<T>* b1,
							  unsigned int acceleration_factor,
							  std::vector<unsigned int> kernel_size,
							  cuNDArray<T>* out_mixing_coeff);

#endif //HTGRAPPA_H
