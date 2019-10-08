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



	return cusparseCcsrmv_mp( handle, transA,m, n, nnz,  (cuComplex*) alpha, descrA,  (cuComplex*) csrValA,  csrRowPtrA,  csrColndA, (cuComplex*) x,  (cuComplex*) beta,  (cuComplex*)y);
}


static cusparseStatus_t  sparseCSRMV(cusparseHandle_t handle, cusparseOperation_t transA,int m, int n, int nnz,
		const complext<double> * alpha, const cusparseMatDescr_t descrA, const complext<double> * csrValA,
		const int * csrRowPtrA, const int * csrColndA, const complext<double> *x, const complext<double>* beta, complext<double>* y){
	return cusparseZcsrmv_mp( handle, transA,m, n, nnz,  (cuDoubleComplex*) alpha, descrA,  (cuDoubleComplex*) csrValA,  csrRowPtrA,  csrColndA, (cuDoubleComplex*) x,  (cuDoubleComplex*) beta,  (cuDoubleComplex*)y);
}

template<class T> EXPORTGPUCORE void Gadgetron::sparseMV(T alpha,T beta, const cuCsrMatrix<T> & mat, const cuNDArray<T> & vec_in, cuNDArray<T>& vec_out, bool adjoint){

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

static cusparseStatus_t sparseCSRMM(cusparseHandle_t handle, cusparseOperation_t transA, int m, int n, int k, int nnz, const float* alpha,
		cusparseMatDescr_t descrA, const float* csrValA, const int* csrRowPtrA, const int* csrColIndA, const float* B,
		int ldb, const float* beta, float* C, int ldc){
	return cusparseScsrmm(handle,transA, m, n, k, nnz, alpha, descrA, csrValA, csrRowPtrA, csrColIndA, B, ldb, beta, C, ldc);
}
static cusparseStatus_t sparseCSRMM(cusparseHandle_t handle, cusparseOperation_t transA, int m, int n, int k, int nnz, const double* alpha,
		cusparseMatDescr_t descrA, const double* csrValA, const int* csrRowPtrA, const int* csrColIndA, const double* B,
		int ldb, const double* beta, double* C, int ldc){
	return cusparseDcsrmm(handle,transA, m, n, k, nnz, alpha, descrA, csrValA, csrRowPtrA, csrColIndA, B, ldb, beta, C, ldc);
}
static cusparseStatus_t sparseCSRMM(cusparseHandle_t handle, cusparseOperation_t transA, int m, int n, int k, int nnz, const complext<float>* alpha,
		cusparseMatDescr_t descrA, const complext<float>* csrValA, const int* csrRowPtrA, const int* csrColIndA, const complext<float>* B,
		int ldb, const complext<float>* beta, complext<float>* C, int ldc){
	return cusparseCcsrmm(handle, transA, m, n, k, nnz, (cuFloatComplex*) alpha, descrA, (cuFloatComplex*) csrValA,
			csrRowPtrA, csrColIndA, (cuFloatComplex*) B, ldb, (cuFloatComplex*) beta, (cuFloatComplex*)C, ldc);
}

static cusparseStatus_t sparseCSRMM(cusparseHandle_t handle, cusparseOperation_t transA, int m, int n, int k, int nnz, const complext<double>* alpha,
		cusparseMatDescr_t descrA, const complext<double>* csrValA, const int* csrRowPtrA, const int* csrColIndA, const complext<double>* B,
		int ldb, const complext<double>* beta, complext<double>* C, int ldc){
	return cusparseZcsrmm(handle, transA, m, n, k, nnz, (cuDoubleComplex*) alpha, descrA, (cuDoubleComplex*) csrValA,
			csrRowPtrA, csrColIndA, (cuDoubleComplex*) B, ldb, (cuDoubleComplex*) beta, (cuDoubleComplex*)C, ldc);
}
template<class T> EXPORTGPUCORE void Gadgetron::sparseMM(T alpha,T beta, const cuCsrMatrix<T> & mat, const cuNDArray<T> & mat_in, cuNDArray<T>& mat_out, bool adjoint) {

	if (mat_in.get_size(1) != mat_out.get_size(1)) throw std::runtime_error("In and out dense matrix must have same second dimension");
	if (mat_in.get_size(0) != mat.n) throw std::runtime_error("Input matrix and sparse matrix have mismatched dimensions");
	if (mat_out.get_size(0) != mat.m) throw std::runtime_error("Output matrix and sparse matrix have mismatched dimensions");

	cusparseOperation_t trans = adjoint ?  CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;

	CUSPARSE_CALL(sparseCSRMM(cudaDeviceManager::Instance()->lockSparseHandle(),trans,
			mat.m,int(mat_in.get_size(1)),mat.n,mat.nnz,&alpha,mat.descr,
			thrust::raw_pointer_cast(&mat.data[0]),thrust::raw_pointer_cast(&mat.csrRow[0]),
			thrust::raw_pointer_cast(&mat.csrColdnd[0]),thrust::raw_pointer_cast(&mat_in.get_device_ptr()[0]),
			int(mat_in.get_size(0)),&beta,thrust::raw_pointer_cast(&mat_out.get_device_ptr()[0]),mat_out.get_size(0)));

	cudaDeviceManager::Instance()->unlockSparseHandle();

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


static cusparseStatus_t sparseCsr2Csc(cusparseHandle_t handle, int m, int n, int nnz,const  float* csrVal,
		const int* csrRowPtr,const int* csrColIndx, float* cscVal, int* cscRowInd, int* cscColPtr,
		cusparseAction_t copyValues, cusparseIndexBase_t idxBase) {
	return cusparseScsr2csc(handle, m, n, nnz, csrVal,csrRowPtr,csrColIndx,cscVal,cscRowInd,cscColPtr,copyValues,idxBase);
}


static cusparseStatus_t sparseCsr2Csc(cusparseHandle_t handle, int m, int n, int nnz,const  double* csrVal,
		const int* csrRowPtr,const int* csrColIndx, double* cscVal, int* cscRowInd, int* cscColPtr,
		cusparseAction_t copyValues, cusparseIndexBase_t idxBase) {
	return cusparseDcsr2csc(handle, m, n, nnz, csrVal,csrRowPtr,csrColIndx,cscVal,cscRowInd,cscColPtr,copyValues,idxBase);
}

static cusparseStatus_t sparseCsr2Csc(cusparseHandle_t handle, int m, int n, int nnz,const complext<float>* csrVal,
		const int* csrRowPtr,const int* csrColIndx, complext<float>* cscVal, int* cscRowInd, int* cscColPtr,
		cusparseAction_t copyValues, cusparseIndexBase_t idxBase) {
	return cusparseCcsr2csc(handle, m, n, nnz, (cuFloatComplex*) csrVal,csrRowPtr,csrColIndx,(cuFloatComplex*)cscVal,cscRowInd,
		 cscColPtr,copyValues,idxBase);
}


static cusparseStatus_t sparseCsr2Csc(cusparseHandle_t handle, int m, int n, int nnz,const complext<double>* csrVal,
		const int* csrRowPtr,const int* csrColIndx, complext<double>* cscVal, int* cscRowInd, int* cscColPtr,
		cusparseAction_t copyValues, cusparseIndexBase_t idxBase) {
	return cusparseZcsr2csc(handle, m, n, nnz, (cuDoubleComplex*) csrVal,csrRowPtr,csrColIndx,(cuDoubleComplex*) cscVal,cscRowInd,
			cscColPtr,copyValues,idxBase);
}

template<class T> EXPORTGPUCORE cuCsrMatrix<T> Gadgetron::transpose(const Gadgetron::cuCsrMatrix<T> &matrix) {
	cuCsrMatrix<T> transposed;
	copysparseMatDescr(transposed.descr,matrix.descr);

	transposed.m = matrix.n;
	transposed.n = matrix.m;
	transposed.nnz = matrix.nnz;
	transposed.data = thrust::device_vector<T>(matrix.nnz);
	transposed.csrRow = thrust::device_vector<int>(transposed.m+1);
	transposed.csrColdnd = thrust::device_vector<int>(transposed.nnz);

	CUSPARSE_CALL(sparseCsr2Csc(cudaDeviceManager::Instance()->lockSparseHandle(),transposed.n,transposed.m,transposed.nnz,
			thrust::raw_pointer_cast(&matrix.data[0]),thrust::raw_pointer_cast(&matrix.csrRow[0]),
			thrust::raw_pointer_cast(&matrix.csrColdnd[0]),thrust::raw_pointer_cast(&transposed.data[0]),
			thrust::raw_pointer_cast(&transposed.csrColdnd[0]),thrust::raw_pointer_cast(&transposed.csrRow[0]),
			CUSPARSE_ACTION_NUMERIC,CUSPARSE_INDEX_BASE_ZERO));

	cudaDeviceManager::Instance()->unlockSparseHandle();
	cudaDeviceSynchronize();

	return transposed;



}


template EXPORTGPUCORE void Gadgetron::sparseMV<float>(float alpha,float beta, const cuCsrMatrix<float> & mat, const cuNDArray<float> & vec_in, cuNDArray<float>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMV<double>(double alpha,double beta, const cuCsrMatrix<double> & mat, const cuNDArray<double> & vec_in, cuNDArray<double>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMV<complext<float> >(complext<float> alpha,complext<float> beta, const cuCsrMatrix<complext<float> > & mat, const cuNDArray<complext<float> > & vec_in, cuNDArray<complext<float> >& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMV<complext<double> >(complext<double> alpha,complext<double> beta, const cuCsrMatrix<complext<double> > & mat, const cuNDArray<complext<double> > & vec_in, cuNDArray<complext<double> >& vec_out, bool adjoint);


template EXPORTGPUCORE Gadgetron::cuCsrMatrix<float> Gadgetron::transpose(const Gadgetron::cuCsrMatrix<float> &matrix);
template EXPORTGPUCORE Gadgetron::cuCsrMatrix<double> Gadgetron::transpose(const Gadgetron::cuCsrMatrix<double> &matrix);
template EXPORTGPUCORE Gadgetron::cuCsrMatrix<complext<float>> Gadgetron::transpose(const Gadgetron::cuCsrMatrix<complext<float>> &matrix);
template EXPORTGPUCORE Gadgetron::cuCsrMatrix<complext<double>> Gadgetron::transpose(const Gadgetron::cuCsrMatrix<complext<double>> &matrix);

template EXPORTGPUCORE void Gadgetron::sparseMM<float>(float alpha,float beta, const cuCsrMatrix<float> & mat, const cuNDArray<float> & vec_in, cuNDArray<float>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMM<double>(double alpha,double beta, const cuCsrMatrix<double> & mat, const cuNDArray<double> & vec_in, cuNDArray<double>& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMM<complext<float> >(complext<float> alpha,complext<float> beta, const cuCsrMatrix<complext<float> > & mat, const cuNDArray<complext<float> > & vec_in, cuNDArray<complext<float> >& vec_out, bool adjoint);
template EXPORTGPUCORE void Gadgetron::sparseMM<complext<double> >(complext<double> alpha,complext<double> beta, const cuCsrMatrix<complext<double> > & mat, const cuNDArray<complext<double> > & vec_in, cuNDArray<complext<double> >& vec_out, bool adjoint);
