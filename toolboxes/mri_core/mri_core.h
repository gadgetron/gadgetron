#ifndef MRI_CORE_H_
#define MRI_CORE_H_

#include "hoNDArray.h"
#include "mri_core_export.h"

#include <complex>

#ifdef USE_OMP
    #include "omp.h"
#endif //USE_OMP

namespace Gadgetron{

// combine (point-wise inner product)
template <typename T> EXPORTMRICORE
void coilmap_combine(const hoNDArray<std::complex<T> > & coilmap, const hoNDArray<std::complex<T> > & data, hoNDArray<std::complex<T> > & result);

// scale (point-wise outer product)
template <typename T> EXPORTMRICORE
void coilmap_scale(const hoNDArray<std::complex<T> > & coilmap, const hoNDArray<std::complex<T> > & rho, hoNDArray<std::complex<T> > & result);

// norm (point-wise 2norm)
template <typename T> EXPORTMRICORE
void coilmap_norm(const hoNDArray<std::complex<T> > & coilmap, hoNDArray<T> & result);


}

#endif //MRI_CORE_H_
