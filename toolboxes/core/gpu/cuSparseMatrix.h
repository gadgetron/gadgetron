/*
 * CUSPARSE.h
 *
 *  Created on: Jan 28, 2015
 *      Author: u051747
 */

#pragma once
#include "cusparse_v2.h"
#include "gpucore_export.h"
#include <thrust/device_vector.h>
#include "cuNDArray.h"
#include "cudaDeviceManager.h"

namespace Gadgetron{

EXPORTGPUCORE std::string gadgetron_getCusparseErrorString(cusparseStatus_t err);

template<class T> struct cuCsrMatrix {

	cuCsrMatrix(){
	CUSPARSE_CALL(cusparseCreateMatDescr(&descr));
	}

	~cuCsrMatrix(){
		CUSPARSE_CALL(cusparseDestroyMatDescr(descr));
	}


	int m,n, nnz;
	thrust::device_vector<int> csrRow, csrColdnd;
	thrust::device_vector<T> data;
	cusparseMatDescr_t descr;
};

/**
 * Performs a sparse matrix vector multiplication: vec_out = alpha*Mat * beta*vec_in
 * @param alpha
 * @param beta
 * @param mat
 * @param vec_in
 * @param vec_out
 * @param adjoint
 */
template<class T> EXPORTGPUCORE void sparseMV(T alpha,T beta, const cuCsrMatrix<T> & mat, const cuNDArray<T> & vec_in, cuNDArray<T>& vec_out, bool adjoint=false);



}
