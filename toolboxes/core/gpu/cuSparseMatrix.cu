#include "cuSparseMatrix.h"
#include <thrust/extrema.h>
#include <thrust/device_ptr.h>
#include "cuNDArray_math.h"
using namespace Gadgetron;

template<class T>
static auto create_DnVec( cuNDArray<T>& vec){

	cusparseDnVecDescr_t dnvec;	
	cusparseCreateDnVec(&dnvec,vec.size(),vec.data(),cuda_datatype<T>());
	auto deleter = [](cusparseDnVecDescr_t val){cusparseDestroyDnVec(val);};

	return std::unique_ptr<std::decay_t<decltype(*dnvec)>,decltype(deleter)>(dnvec,deleter);
}



template<class T> EXPORTGPUCORE void Gadgetron::sparseMV(T alpha,T beta, const cuCsrMatrix<T> & mat, const cuNDArray<T> & vec_in, cuNDArray<T>& vec_out, bool adjoint){

	if (vec_in.get_number_of_elements() != (adjoint ? mat.rows : mat.cols))
		throw std::runtime_error("Matrix and input vector have mismatching dimensions");
	if (vec_out.get_number_of_elements() != (adjoint ? mat.rows : mat.cols))
		throw std::runtime_error("Matrix and output vector have mismatching dimensions");

	cusparseOperation_t trans = adjoint ?  CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;
	//cusparseStatus_t status = sparseCSRMV(cudaDeviceManager::Instance()->lockSparseHandle(),trans,mat.m,mat.n,mat.nnz,&alpha, mat.descr,
	//		thrust::raw_pointer_cast(&mat.data[0]),thrust::raw_pointer_cast(&mat.csrRow[0]),thrust::raw_pointer_cast(&mat.csrColdnd[0]),vec_in.get_data_ptr(),&beta,vec_out.get_data_ptr());


	auto dnvec_in = create_DnVec(const_cast<cuNDArray<T>&>(vec_in));
	auto dnvec_out = create_DnVec(vec_out);

	size_t bufferSize;
	auto handle =  cudaDeviceManager::Instance()->lockSparseHandle();
	cusparseSpMV(handle, trans, &alpha,mat.descr,dnvec_in.get(),&beta,dnvec_out.get(),cuda_datatype<T>(),CUSPARSE_CSRMV_ALG2,&bufferSize);
	cuNDArray<char> buffer(bufferSize);

	cusparseStatus_t status = cusparseSpMV(handle, trans, &alpha,mat.descr,dnvec_in.get(),&beta,dnvec_out.get(),cuda_datatype<T>(),CUSPARSE_CSRMV_ALG2, buffer.data());

	cudaDeviceManager::Instance()->unlockSparseHandle();
	if (status != CUSPARSE_STATUS_SUCCESS){
		std::stringstream ss;
		ss << "Sparse Matrix Vector multiplication failed. Error: ";
		ss << gadgetron_getCusparseErrorString(status);
		throw cuda_error(ss.str());
	}



}
template<class T>
static auto create_DnMat( cuNDArray<T>& mat){

	cusparseDnMatDescr_t dnmat;	
	cusparseCreateDnMat(&dnmat,mat.get_size(0),mat.get_size(1),mat.get_size(0),mat.data(),cuda_datatype<T>(), CUSPARSE_ORDER_COL);
	auto deleter = [](cusparseDnMatDescr_t val){cusparseDestroyDnMat(val);};
	return std::unique_ptr<std::decay_t<decltype(*dnmat)>,decltype(deleter)>(dnmat,deleter);
	//return std::unique_ptr<std::decay_t<decltype(*dnmat)>,decltype(&cusparseDestroyDnMat)>(dnmat);
}



template<class T> EXPORTGPUCORE void Gadgetron::sparseMM(T alpha,T beta, const cuCsrMatrix<T> & mat, const cuNDArray<T> & mat_in, cuNDArray<T>& mat_out, bool adjoint) {

	if (mat_in.get_size(1) != mat_out.get_size(1)) throw std::runtime_error("In and out dense matrix must have same second dimension");
	if (mat_in.get_size(0) != mat.rows) throw std::runtime_error("Input matrix and sparse matrix have mismatched dimensions");
	if (mat_out.get_size(0) != mat.cols) throw std::runtime_error("Output matrix and sparse matrix have mismatched dimensions");

	cusparseOperation_t trans = adjoint ?  CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;
	auto handle = cudaDeviceManager::Instance()->lockSparseHandle();

	auto dnmat_in = create_DnMat(const_cast<cuNDArray<T>&>(mat_in));
	auto dnmat_out = create_DnMat(mat_out);
	size_t bufferSize;
	CUSPARSE_CALL(cusparseSpMM_bufferSize(handle, trans, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat.descr, dnmat_in.get(), &beta, dnmat_out.get(), cuda_datatype<T>(),CUSPARSE_CSRMM_ALG1, &bufferSize));
	cuNDArray<char> buffer(bufferSize);

	CUSPARSE_CALL(cusparseSpMM(handle, trans, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat.descr, dnmat_in.get(), &beta, dnmat_out.get(), cuda_datatype<T>(),CUSPARSE_CSRMM_ALG1, buffer.data()));
	cudaDeviceManager::Instance()->unlockSparseHandle();

}


template EXPORTGPUCORE void Gadgetron::sparseMV<float>(float alpha,float beta, const cuCsrMatrix<float> & mat, const cuNDArray<float> & vec_in, cuNDArray<float>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMV<double>(double alpha,double beta, const cuCsrMatrix<double> & mat, const cuNDArray<double> & vec_in, cuNDArray<double>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMV<complext<float> >(complext<float> alpha,complext<float> beta, const cuCsrMatrix<complext<float> > & mat, const cuNDArray<complext<float> > & vec_in, cuNDArray<complext<float> >& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMV<complext<double> >(complext<double> alpha,complext<double> beta, const cuCsrMatrix<complext<double> > & mat, const cuNDArray<complext<double> > & vec_in, cuNDArray<complext<double> >& vec_out, bool adjoint);

template EXPORTGPUCORE void Gadgetron::sparseMM<float>(float alpha,float beta, const cuCsrMatrix<float> & mat, const cuNDArray<float> & vec_in, cuNDArray<float>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMM<double>(double alpha,double beta, const cuCsrMatrix<double> & mat, const cuNDArray<double> & vec_in, cuNDArray<double>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMM<complext<float> >(complext<float> alpha,complext<float> beta, const cuCsrMatrix<complext<float> > & mat, const cuNDArray<complext<float> > & vec_in, cuNDArray<complext<float> >& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMM<complext<double> >(complext<double> alpha,complext<double> beta, const cuCsrMatrix<complext<double> > & mat, const cuNDArray<complext<double> > & vec_in, cuNDArray<complext<double> >& vec_out, bool adjoint);
