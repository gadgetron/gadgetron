/*
 * matrix_decomposition.h
 *
 *  Created on: Dec 10, 2011
 *      Author: Michael S. Hansen
 */

#ifndef MATRIX_DECOMPOSITION_H_
#define MATRIX_DECOMPOSITION_H_

#include "hoNDArray.h"
#include <complex>
#include "complext.h"

#include "linalg_export.h"

namespace Gadgetron 
{

/**
 *   Perform Cholesky decomposition of matrix.
 *   hoNDArray should be symmetric/hermitian positive definite.
 *   Matrix will be replaced with lower triangular matrix.
 *   Calls LAPACK subroutine _POTRF
 *
 */
template <typename T> EXPORTLINALG void hoNDArray_choldc(hoNDArray<T>* A);

/**
 * Invert matrix assuming it is lower trinagular
 */
template <typename T> EXPORTLINALG void hoNDArray_inv_lower_triangular(hoNDArray<T>* A);

/**
 *  SVD
 */
template <typename T> EXPORTLINALG void hoNDArray_svd(hoNDArray< T >* A, hoNDArray< T >* U, hoNDArray<typename realType<T>::Type>* S, hoNDArray< T >* VT);

} //Namespace Gadgetron


#endif /* MATRIX_DECOMPOSITION_H_ */
