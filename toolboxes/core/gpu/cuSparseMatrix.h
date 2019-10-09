/*
 * CUSPARSE.h
 *
 *  Created on: Jan 28, 2015
 *      Author: u051747
 */

#pragma once
#include "cusparse.h"
#include "gpucore_export.h"
#include <thrust/device_vector.h>
#include "cuNDArray.h"
#include "cudaDeviceManager.h"

namespace Gadgetron{



EXPORTGPUCORE std::string gadgetron_getCusparseErrorString(cusparseStatus_t err);

inline void copysparseMatDescr(cusparseMatDescr_t& destination, const cusparseMatDescr_t& source){
    CUSPARSE_CALL(cusparseSetMatType(destination,cusparseGetMatType(source)));
    CUSPARSE_CALL(cusparseSetMatFillMode(destination,cusparseGetMatFillMode(source)));
    CUSPARSE_CALL(cusparseSetMatDiagType(destination,cusparseGetMatDiagType(source)));
    CUSPARSE_CALL(cusparseSetMatIndexBase(destination,cusparseGetMatIndexBase(source)));
}

template<class T> struct cuCsrMatrix {

	cuCsrMatrix(){
	CUSPARSE_CALL(cusparseCreateMatDescr(&descr));
	}

	~cuCsrMatrix(){
		cusparseDestroyMatDescr(descr);
	}


	cuCsrMatrix(cuCsrMatrix&& other){
		CUSPARSE_CALL(cusparseCreateMatDescr(&descr));
		*this = std::move(other);
	}

	cuCsrMatrix(const cuCsrMatrix& other){
		CUSPARSE_CALL(cusparseCreateMatDescr(&descr));
		*this = other;
	}

	cuCsrMatrix& operator=(cuCsrMatrix&& other){
		copysparseMatDescr(descr,other.descr);
		csrRow = std::move(other.csrRow);
		csrColdnd = std::move(other.csrColdnd);
		data = std::move(other.data);
		m = other.m;
		n = other.n;
		nnz = other.nnz;
		return *this;
	}

	cuCsrMatrix& operator=(const cuCsrMatrix& other){
		copysparseMatDescr(descr,other.descr);
		csrRow = other.csrRow;
		csrColdnd = other.csrColdnd;
		data = other.data;
		m = other.m;
		n = other.n;
		nnz = other.nnz;
		return *this;
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

template<class T> EXPORTGPUCORE void sparseMM(T alpha,T beta, const cuCsrMatrix<T> & mat, const cuNDArray<T> & mat_in, cuNDArray<T>& mat_out, bool adjoint=false);

template<class T> EXPORTGPUCORE cuCsrMatrix<T> transpose(const cuCsrMatrix<T>& matrix);

}
