#include "cuSparseMatrix.h"
#include <thrust/extrema.h>
#include <thrust/device_ptr.h>
#include "cuNDArray_math.h"
using namespace Gadgetron;


static cusparseStatus_t  sparseCSRMV(cusparseHandle_t handle, cusparseOperation_t transA,int m, int n, int nnz,
		const float * alpha, const cusparseMatDescr_t descrA, const float * csrValA,
		const int * csrRowPtrA, const int * csrColndA, const float *x, const float* beta, float* y){

	return cusparseScsrmv( handle, transA,m, n, nnz,  alpha, descrA,  csrValA,  csrRowPtrA,  csrColndA, x,  beta, y);
}


static cusparseStatus_t  sparseCSRMV(cusparseHandle_t handle, cusparseOperation_t transA,int m, int n, int nnz,
		const double * alpha, const cusparseMatDescr_t descrA, const double * csrValA,
		const int * csrRowPtrA, const int * csrColndA, const double *x, const double* beta, double* y){

	return cusparseDcsrmv( handle, transA,m, n, nnz,  alpha, descrA,  csrValA,  csrRowPtrA,  csrColndA, x,  beta, y);
}


static cusparseStatus_t  sparseCSRMV(cusparseHandle_t handle, cusparseOperation_t transA,int m, int n, int nnz,
		const complext<float> * alpha, const cusparseMatDescr_t descrA, const complext<float> * csrValA,
		const int * csrRowPtrA, const int * csrColndA, const complext<float> *x, const complext<float>* beta, complext<float>* y){

	thrust::device_ptr<const int> csrRow(csrRowPtrA);

	thrust::device_ptr<const int> csrCol(csrColndA);

	thrust::device_ptr<const complext<float> >  csrVal(csrValA);

	thrust::device_ptr<const complext<float> > xptr(x);
	thrust::device_ptr<complext<float> > yptr(y);


	if (transA == CUSPARSE_OPERATION_NON_TRANSPOSE){
		std::cout << "In sum " << thrust::reduce(xptr,xptr+n) << " out sum " << thrust::reduce(yptr,yptr+m) << std::endl;
	} else {
		std::cout << "T In sum " << thrust::reduce(xptr,xptr+m) << " out sum " << thrust::reduce(yptr,yptr+n) << std::endl;
	}



	return cusparseCcsrmv( handle, transA,m, n, nnz,  (cuComplex*) alpha, descrA,  (cuComplex*) csrValA,  csrRowPtrA,  csrColndA, (cuComplex*) x,  (cuComplex*) beta,  (cuComplex*)y);
}


static cusparseStatus_t  sparseCSRMV(cusparseHandle_t handle, cusparseOperation_t transA,int m, int n, int nnz,
		const complext<double> * alpha, const cusparseMatDescr_t descrA, const complext<double> * csrValA,
		const int * csrRowPtrA, const int * csrColndA, const complext<double> *x, const complext<double>* beta, complext<double>* y){
	return cusparseZcsrmv( handle, transA,m, n, nnz,  (cuDoubleComplex*) alpha, descrA,  (cuDoubleComplex*) csrValA,  csrRowPtrA,  csrColndA, (cuDoubleComplex*) x,  (cuDoubleComplex*) beta,  (cuDoubleComplex*)y);
}

template<class T> EXPORTGPUCORE void Gadgetron::sparseMV(T alpha,T beta, const cuCsrMatrix<T> & mat, const cuNDArray<T> & vec_in, cuNDArray<T>& vec_out, bool adjoint=false){

	if (vec_in.get_number_of_elements() != (adjoint ? mat.m : mat.n))
		throw std::runtime_error("Matrix and input vector have mismatching dimensions");
	if (vec_out.get_number_of_elements() != (adjoint ? mat.n : mat.m))
		throw std::runtime_error("Matrix and output vector have mismatching dimensions");

	cusparseOperation_t trans = adjoint ?  CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;
	cusparseStatus_t status = sparseCSRMV(cudaDeviceManager::Instance()->lockSparseHandle(),trans,mat.m,mat.n,mat.nnz,&alpha, mat.descr,
			thrust::raw_pointer_cast(&mat.data[0]),thrust::raw_pointer_cast(&mat.csrRow[0]),thrust::raw_pointer_cast(&mat.csrColdnd[0]),vec_in.get_data_ptr(),&beta,vec_out.get_data_ptr());

	cudaDeviceManager::Instance()->unlockSparseHandle();
	if (status != CUSPARSE_STATUS_SUCCESS){
		std::stringstream ss;
		ss << "Sparse Matrix Vector multiplication failed. Error: ";
		ss << gadgetron_getCusparseErrorString(status);
		throw cuda_error(ss.str());
	}



}




EXPORTGPUCORE std::string Gadgetron::gadgetron_getCusparseErrorString(cusparseStatus_t err)
{
  switch (err){
  case CUSPARSE_STATUS_NOT_INITIALIZED:
    return "NOT INITIALIZED";
  case CUSPARSE_STATUS_ALLOC_FAILED:
    return "ALLOC FAILED";
  case CUSPARSE_STATUS_INVALID_VALUE:
    return "INVALID VALUE";
  case CUSPARSE_STATUS_ARCH_MISMATCH:
    return "ARCH MISMATCH";
  case CUSPARSE_STATUS_MAPPING_ERROR:
    return "MAPPING ERROR";
  case CUSPARSE_STATUS_EXECUTION_FAILED:
    return "EXECUTION FAILED";
  case CUSPARSE_STATUS_INTERNAL_ERROR:
    return "INTERNAL ERROR";
  case CUSPARSE_STATUS_SUCCESS:
    return "SUCCES";
  case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
  	return "MATRIX TYPE NOT SUPPORTED";
  default:
    return "UNKNOWN CUSPARSE ERROR";
  }
}



template EXPORTGPUCORE void Gadgetron::sparseMV<float>(float alpha,float beta, const cuCsrMatrix<float> & mat, const cuNDArray<float> & vec_in, cuNDArray<float>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMV<double>(double alpha,double beta, const cuCsrMatrix<double> & mat, const cuNDArray<double> & vec_in, cuNDArray<double>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMV<complext<float> >(complext<float> alpha,complext<float> beta, const cuCsrMatrix<complext<float> > & mat, const cuNDArray<complext<float> > & vec_in, cuNDArray<complext<float> >& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMV<complext<double> >(complext<double> alpha,complext<double> beta, const cuCsrMatrix<complext<double> > & mat, const cuNDArray<complext<double> > & vec_in, cuNDArray<complext<double> >& vec_out, bool adjoint);
