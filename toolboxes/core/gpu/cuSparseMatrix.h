/*
 * CUSPARSE.h
 *
 *  Created on: Jan 28, 2015
 *      Author: u051747
 */

#pragma once
#include "cusparse.h"
#include <thrust/device_vector.h>
#include "cuNDArray.h"
#include "cudaDeviceManager.h"

namespace Gadgetron
{


	template <class T>
	struct cuCsrMatrix
	{

		cuCsrMatrix(size_t rows, size_t cols, thrust::device_vector<int> csrRow, thrust::device_vector<int> csrColdnd, thrust::device_vector<T> data) : csrRow{std::move(csrRow)}, csrColdnd{std::move(csrColdnd)}, data{std::move(data)}, rows{rows}, cols{cols}
		{
			cusparseCreateCsr(&descr, rows, cols, this->data.size(), 
				thrust::raw_pointer_cast(this->csrRow.data()), thrust::raw_pointer_cast(this->csrColdnd.data()), thrust::raw_pointer_cast(this->data.data()), 
				CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, Gadgetron::cuda_datatype<T>());
		}

		~cuCsrMatrix()
		{
			if (this->descr)
				cusparseDestroySpMat(this->descr);
		}

		cuCsrMatrix(cuCsrMatrix &&other)
		{
			*this = std::move(other);
		}

		cuCsrMatrix &operator=(cuCsrMatrix &&other)
		{
			this->descr = other.descr;
			other.descr = nullptr;
			this->csrColdnd = std::move(other.csrColdnd);
			this->csrRow = std::move(other.csrRow);
			this->data = std::move(this->data);
			return *this;
		}

		size_t rows, cols;
		thrust::device_vector<int> csrRow, csrColdnd;
		thrust::device_vector<T> data;
		cusparseSpMatDescr_t descr;
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
	template <class T>
	void sparseMV(T alpha, T beta, const cuCsrMatrix<T> &mat, const cuNDArray<T> &vec_in, cuNDArray<T> &vec_out, bool adjoint = false);

	template <class T>
	void sparseMM(T alpha, T beta, const cuCsrMatrix<T> &mat, const cuNDArray<T> &mat_in, cuNDArray<T> &mat_out, bool adjoint = false);

} // namespace Gadgetron
