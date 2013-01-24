/*
 * matrix_decomposition.cpp
 *
 *  Created on: Dec 10, 2011
 *      Author: Michael S. Hansen
 */

#include "matrix_decomposition.h"
#include <complex>

void potrf_wrapper(char* UPLO, int* N, float* A, int* LDA, int* info)
{
	spotrf_(UPLO, N, A, LDA, info);
}

void potrf_wrapper(char* UPLO, int* N, double* A, int* LDA, int* info)
{
	dpotrf_(UPLO, N, A, LDA, info);
}

void potrf_wrapper(char* UPLO, int* N, std::complex<float>* A, int* LDA, int* info)
{
	cpotrf_(UPLO, N, A, LDA, info);
}

void potrf_wrapper(char* UPLO, int* N, std::complex<double>* A, int* LDA, int* info)
{
	zpotrf_(UPLO, N, A, LDA, info);
}

template <typename T> void hoNDArray_choldc(hoNDArray<T>* A)
{
	/*
	 *  We are specifying Upper Triangular,
	 *  but matrix comes in transposed (row-major) compared to
	 *  Fortran column-major order. As a result, we will get the lower
	 *  triangular matrix.
	 */
	char UPLO = 'U';
	if (A->get_number_of_dimensions() != 2) {
		BOOST_THROW_EXCEPTION(gt_runtime_error("This is not a matrix, only two dimensions allowed"));
	}

	int N = A->get_size(0);
	if (N != A->get_size(1)) {
		BOOST_THROW_EXCEPTION(gt_runtime_error("Matrix is not symmetric."));
	}

	int LDA = N;
	int info = 0;

	potrf_wrapper(&UPLO, &N, A->get_data_ptr(), &LDA, &info);

	if (info != 0) {
		BOOST_THROW_EXCEPTION(gt_runtime_error("Error calling _potrf wrapper routine."));
	}

	/* Temp code to zero upper triangular */
	/*
	T* d = A->get_data_ptr();
	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = i+1; j < N; j++) {
			d[i*N+j] = 0;
		}
	}
	*/

}


//Template instanciations
template EXPORTLINALG void hoNDArray_choldc(hoNDArray< std::complex<float> >* A);
template EXPORTLINALG void hoNDArray_choldc(hoNDArray< std::complex<double> >* A);
template EXPORTLINALG void hoNDArray_choldc(hoNDArray< float >* A);
template EXPORTLINALG void hoNDArray_choldc(hoNDArray< double >* A);


void trtri_wrapper(char* UPLO, char* DIAG, int* N, float* A, int* LDA, int* info)
{
	strtri_(UPLO, DIAG, N, A, LDA, info);
}

void trtri_wrapper(char* UPLO, char* DIAG, int* N, double* A, int* LDA, int* info)
{
	dtrtri_(UPLO, DIAG, N, A, LDA, info);
}

void trtri_wrapper(char* UPLO, char* DIAG, int* N, std::complex<float>* A, int* LDA, int* info)
{
	ctrtri_(UPLO, DIAG, N, A, LDA, info);
}

void trtri_wrapper(char* UPLO, char* DIAG, int* N, std::complex<double>* A, int* LDA, int* info)
{
	ztrtri_(UPLO, DIAG, N, A, LDA, info);
}

template <typename T> void hoNDArray_inv_lower_triangular(hoNDArray<T>* A)
{
	const char* fname = "hoNDArray_inv_lower_triangular(hoNDArray<T>* A)";

	if (A->get_number_of_dimensions() != 2) {
		BOOST_THROW_EXCEPTION(gt_runtime_error("Error array is not 2 dimensional."));
	}

	int N = A->get_size(0);

	if (N != A->get_size(1)) {
		BOOST_THROW_EXCEPTION(gt_runtime_error("Error array is not 2 dimensional."));
	}

	int LDA = N;
	char UPLO = 'U'; //We are passing upper, but matrix is really lower. This is do deal with row and column major order differences
	char DIAG = 'N';
	int info;

	trtri_wrapper(&UPLO, &DIAG, &N, A->get_data_ptr(), &LDA, &info);

	if (info != 0) {
		BOOST_THROW_EXCEPTION(gt_runtime_error("Error inverting triangular matrix."));
	}

}

template EXPORTLINALG void hoNDArray_inv_lower_triangular(hoNDArray<float>* A);
template EXPORTLINALG void hoNDArray_inv_lower_triangular(hoNDArray<double>* A);
template EXPORTLINALG void hoNDArray_inv_lower_triangular(hoNDArray< std::complex<float> >* A);
template EXPORTLINALG void hoNDArray_inv_lower_triangular(hoNDArray< std::complex<double> >* A);


template<typename T>
boost::shared_ptr<hoNDArray<T> > hoNDArray_transpose(hoNDArray<T> *A, bool copy_data = true)
{
	const char* fname = "hoNDArray_transpose(hoNDArray<T> *A, bool copy_data = true)";

	boost::shared_ptr<hoNDArray<T> > ret_val;

	if (A->get_number_of_dimensions() != 2) {
		BOOST_THROW_EXCEPTION(gt_runtime_error("Error array is not 2 dimensional."));
		return ret_val;
	}

	std::vector<unsigned int> permute_order(2);
	permute_order[0] = 1;permute_order[1] = 0;

	std::vector<unsigned int> perm_dims(2);
	perm_dims[0] = A->get_size(1);
	perm_dims[1] = A->get_size(0);

	ret_val.reset(new hoNDArray<T>);
	ret_val.get()->create(&perm_dims);

	if (copy_data) {
		A->permute(&permute_order, ret_val.get());
	}
	return ret_val;
}

void gesvd_wrapper(char* JOBU, char* JOBVT, int* M, int* N, float* A,
		int* LDA, float* S, float* U, int* LDU, float* VT,
		int* LDVT, float* WORK, int* LWORK, float* RWORK, int* INFO)
{
	sgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
}

void gesvd_wrapper(char* JOBU, char* JOBVT, int* M, int* N, double* A,
		int* LDA, double* S, double* U, int* LDU, double* VT,
		int* LDVT, double* WORK, int* LWORK, double* RWORK, int* INFO)
{
	dgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
}

void gesvd_wrapper(char* JOBU, char* JOBVT, int* M, int* N, std::complex<float>* A,
		int* LDA, float* S, std::complex<float>* U, int* LDU, std::complex<float>* VT,
		int* LDVT, std::complex<float>* WORK, int* LWORK, float* RWORK, int* INFO)
{
	cgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
}

void gesvd_wrapper(char* JOBU, char* JOBVT, int* M, int* N, std::complex<double>* A,
		int* LDA, double* S, std::complex<double>* U, int* LDU, std::complex<double>* VT,
		int* LDVT, std::complex<double>* WORK, int* LWORK, double* RWORK, int* INFO)
{
	zgesvd_(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO);
}

template<typename T> void hoNDArray_svd(hoNDArray<T > *A,
		hoNDArray< T > *U, hoNDArray< typename realType<T>::type > *S, hoNDArray< T > *VT)
{


	if (A->get_number_of_dimensions() != 2) {
		BOOST_THROW_EXCEPTION(gt_runtime_error(" Error array A is not 2 dimensional."));
	}

	boost::shared_ptr< hoNDArray< T > > A_t = hoNDArray_transpose(A);
	if (!A_t.get()) {
		BOOST_THROW_EXCEPTION(gt_runtime_error("Transpose of input failed."));
	}

	int M = A_t->get_size(0);
	int N = A_t->get_size(1);
    int min_M_N = M > N ? N : M;
    int max_M_N = M < N ? N : M;

	T* A_ptr = A_t.get()->get_data_ptr();
	if (!A_ptr) {
		BOOST_THROW_EXCEPTION(gt_runtime_error( "Data array pointer is undefined."));

	}


	boost::shared_ptr< hoNDArray< T > > U_t;
	boost::shared_ptr< hoNDArray< T > > VT_t;

	char JOBU, JOBVT;
	T* U_ptr = 0;
	T* VT_ptr = 0;
	typename realType<T>::type * S_ptr = 0;

	if (S) {
		if (S->get_number_of_elements() < min_M_N) {
			BOOST_THROW_EXCEPTION(gt_runtime_error("S is too small."));

		}
		S_ptr = S->get_data_ptr();
	} else {
		BOOST_THROW_EXCEPTION(gt_runtime_error("Null pointer detected for S."));

	}

	int LDU = 1;
	if (U) {
		if (U->get_number_of_dimensions() != 2) {
			BOOST_THROW_EXCEPTION(gt_runtime_error("Error array U is not 2 dimensional."));

		}

		U_t = hoNDArray_transpose(U, false);

		if (U_t.get()->get_size(0) != M) {
			BOOST_THROW_EXCEPTION(gt_runtime_error("Number of rows in U is not equal to number of rows in A"));
		}

		if (U_t.get()->get_size(1) == M) {
			JOBU = 'A';
		} else if (U_t.get()->get_size(1) == min_M_N) {
			JOBU = 'S';
		} else {
			BOOST_THROW_EXCEPTION(gt_runtime_error("Invalid number of columns of U"));
		}

		U_ptr = U_t.get()->get_data_ptr();
		LDU = U_t.get()->get_size(0);
	} else {
		JOBU = 'N';
	}

	int LDVT = 1;
	if (VT) {
		if (VT->get_number_of_dimensions() != 2) {
			BOOST_THROW_EXCEPTION(gt_runtime_error("Error array VT is not 2 dimensional."));
		}

		VT_t = hoNDArray_transpose(VT, false);

		if (VT_t.get()->get_size(0) == N) {
			JOBVT = 'A';
		} else if (VT_t.get()->get_size(0) == min_M_N) {
			JOBVT = 'S';
		} else {
			BOOST_THROW_EXCEPTION(gt_runtime_error("Invalid number of rows of VT"));
		}

		VT_ptr = VT_t.get()->get_data_ptr();
		LDVT = VT_t.get()->get_size(0);

	} else {
		JOBVT = 'N';
	}

	//Lets allocate some work storage
	std::vector<unsigned int> work_dim(1);

	int LWORK = 5*2*min_M_N+max_M_N;
	work_dim[0] = LWORK;

	hoNDArray< T > WORK(&work_dim);
	work_dim[0] = 5*min_M_N;
	hoNDArray< typename realType<T>::type > RWORK(&work_dim);

	//Now we are finally ready to call the SVD
	int INFO;

	gesvd_wrapper(&JOBU, &JOBVT, &M, &N, A_ptr,
			&M, S_ptr, U_ptr, &LDU, VT_ptr,
			&LDVT, WORK.get_data_ptr(),
			&LWORK, RWORK.get_data_ptr(), &INFO);

	if (INFO != 0) {
		std::stringstream ss;
		ss << "Call to gesvd failed, INFO = " << INFO << "";
		BOOST_THROW_EXCEPTION(gt_runtime_error(ss.str()));
	}

	std::vector<unsigned int> permute_order(2);
	permute_order[0] = 1;permute_order[1] = 0;

	if (U) {
		U_t.get()->permute(&permute_order,U);
	}

	if (VT) {
		VT_t.get()->permute(&permute_order,VT);
	}
}

//Template instanciations
template EXPORTLINALG void hoNDArray_svd(hoNDArray< float > *A, hoNDArray< float > *U, hoNDArray< float > *S, hoNDArray< float > *VT);
template EXPORTLINALG void hoNDArray_svd(hoNDArray< double > *A, hoNDArray< double > *U, hoNDArray< double > *S, hoNDArray< double > *VT);
template EXPORTLINALG void hoNDArray_svd(hoNDArray< std::complex<float> > *A, hoNDArray< std::complex<float> > *U, hoNDArray< float > *S, hoNDArray< std::complex<float> > *VT);
template EXPORTLINALG void hoNDArray_svd(hoNDArray< std::complex<double> > *A, hoNDArray< std::complex<double> > *U, hoNDArray< double > *S, hoNDArray< std::complex<double> > *VT);

