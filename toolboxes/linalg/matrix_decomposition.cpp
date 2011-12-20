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

template <typename T> int hoNDArray_choldc(hoNDArray<T>* A)
{
	const char* fname = "hoNDArray_choldc(hoNDArray<T>* A)";
	/*
	 *  We are specifying Upper Triangular,
	 *  but matrix comes in transposed (row-major) compared to
	 *  Fortran column-major order. As a result, we will get the lower
	 *  triangular matrix.
	 */
	char UPLO = 'U';
	if (A->get_number_of_dimensions() != 2) {
		std::cout << fname << ": This is not a matrix, only two dimensions allowed\n" << std::endl;
		return -1;
	}

	int N = A->get_size(0);
	if (N != A->get_size(1)) {
		std::cout << fname << ": Matrix is not symmetric.\n" << std::endl;
		return -1;
	}

	int LDA = N;
	int info = 0;

	potrf_wrapper(&UPLO, &N, A->get_data_ptr(), &LDA, &info);

	if (info != 0) {
		std::cout << fname << ": Error calling _potrf wrapper routine.\n" << std::endl;
		return -1;
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

	return info;
}


//Template instanciations
template EXPORTLINALG int hoNDArray_choldc(hoNDArray< std::complex<float> >* A);
template EXPORTLINALG int hoNDArray_choldc(hoNDArray< std::complex<double> >* A);
template EXPORTLINALG int hoNDArray_choldc(hoNDArray< float >* A);
template EXPORTLINALG int hoNDArray_choldc(hoNDArray< double >* A);


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

template <typename T> int hoNDArray_inv_lower_triangular(hoNDArray<T>* A)
{
	const char* fname = "hoNDArray_inv_lower_triangular(hoNDArray<T>* A)";

	if (A->get_number_of_dimensions() != 2) {
		std::cout << fname << ": Error array is not 2 dimensional.\n" << std::endl;
		return -1;
	}

	int N = A->get_size(0);

	if (N != A->get_size(1)) {
		std::cout << fname << ": Error array is not 2 dimensional.\n" << std::endl;
		return -1;
	}

	int LDA = N;
	char UPLO = 'U'; //We are passing upper, but matrix is really lower. This is do deal with row and column major order differences
	char DIAG = 'N';
	int info;

	trtri_wrapper(&UPLO, &DIAG, &N, A->get_data_ptr(), &LDA, &info);

	if (info != 0) {
		std::cout << fname << ": Error inverting triangular matrix.\n" << std::endl;
		return -1;
	}

	return 0;
}

template EXPORTLINALG int hoNDArray_inv_lower_triangular(hoNDArray<float>* A);
template EXPORTLINALG int hoNDArray_inv_lower_triangular(hoNDArray<double>* A);
template EXPORTLINALG int hoNDArray_inv_lower_triangular(hoNDArray< std::complex<float> >* A);
template EXPORTLINALG int hoNDArray_inv_lower_triangular(hoNDArray< std::complex<double> >* A);


template<typename T>
boost::shared_ptr<hoNDArray<T> > hoNDArray_transpose(hoNDArray<T> *A, bool copy_data = true)
{
	const char* fname = "hoNDArray_transpose(hoNDArray<T> *A, bool copy_data = true)";

	boost::shared_ptr<hoNDArray<T> > ret_val;

	if (A->get_number_of_dimensions() != 2) {
		std::cout << fname << ": Error array is not 2 dimensional.\n" << std::endl;
		return ret_val;
	}

	std::vector<unsigned int> permute_order(2);
	permute_order[0] = 1;permute_order[1] = 0;

	std::vector<unsigned int> perm_dims(2);
	perm_dims[0] = A->get_size(1);
	perm_dims[1] = A->get_size(0);

	ret_val.reset(new hoNDArray<T>);
	if (!ret_val.get()->create(&perm_dims)) {
		std::cout << fname << ": Unable to allocate transposed array.\n" << std::endl;
		ret_val.reset();
		return ret_val;
	}

	if (copy_data) {
		if (A->permute(&permute_order, ret_val.get()) != 0) {
			std::cout << fname << ": Unable to transpose array.\n" << std::endl;
			ret_val.reset();
			return ret_val;
		}
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

template<typename T, typename Y> int hoNDArray_svd(hoNDArray<T > *A,
		hoNDArray< T > *U, hoNDArray< Y > *S, hoNDArray< T > *VT)
{
	const char* fname = "hoNDArray_svd(hoNDArray *A, hoNDArray *U, hoNDArray *S, hoNDArray *VT)";

	if (A->get_number_of_dimensions() != 2) {
		std::cout << fname << ": Error array A is not 2 dimensional.\n" << std::endl;
		return -1;
	}

	boost::shared_ptr< hoNDArray< T > > A_t = hoNDArray_transpose(A);
	if (!A_t.get()) {
		std::cout << fname << ": Transpose of input failed.\n" << std::endl;
		return -1;
	}

	int M = A_t->get_size(0);
	int N = A_t->get_size(1);
    int min_M_N = M > N ? N : M;
    int max_M_N = M < N ? N : M;

	T* A_ptr = A_t.get()->get_data_ptr();
	if (!A_ptr) {
		std::cout << fname << ": Data array pointer is undefined.\n" << std::endl;
		return -1;
	}


	boost::shared_ptr< hoNDArray< T > > U_t;
	boost::shared_ptr< hoNDArray< T > > VT_t;

	char JOBU, JOBVT;
	T* U_ptr = 0;
	T* VT_ptr = 0;
	Y* S_ptr = 0;

	if (S) {
		if (S->get_number_of_elements() < min_M_N) {
			std::cout << fname << ": S is too small.\n" << std::endl;
			return -1;
		}
		S_ptr = S->get_data_ptr();
	} else {
		std::cout << fname << ": Null pointer detected for S.\n" << std::endl;
		return -1;
	}

	int LDU = 1;
	if (U) {
		if (U->get_number_of_dimensions() != 2) {
			std::cout << fname << ": Error array U is not 2 dimensional.\n" << std::endl;
			return -1;
		}

		U_t = hoNDArray_transpose(U, false);

		if (U_t.get()->get_size(0) != M) {
			std::cout << fname << ": Number of rows in U is not equal to number of rows in A\n" << std::endl;
			return -1;
		}

		if (U_t.get()->get_size(1) == M) {
			JOBU = 'A';
		} else if (U_t.get()->get_size(1) == min_M_N) {
			JOBU = 'S';
		} else {
			std::cout << fname << ": Invalid number of columns of U\n" << std::endl;
			return -1;
		}

		U_ptr = U_t.get()->get_data_ptr();
		LDU = U_t.get()->get_size(0);
	} else {
		JOBU = 'N';
	}

	int LDVT = 1;
	if (VT) {
		if (VT->get_number_of_dimensions() != 2) {
			std::cout << fname << ": Error array VT is not 2 dimensional.\n" << std::endl;
			return -1;
		}

		VT_t = hoNDArray_transpose(VT, false);

		if (VT_t.get()->get_size(0) == N) {
			JOBVT = 'A';
		} else if (VT_t.get()->get_size(0) == min_M_N) {
			JOBVT = 'S';
		} else {
			std::cout << fname << ": Invalid number of rows of VT\n" << std::endl;
			return -1;
		}

		VT_ptr = VT_t.get()->get_data_ptr();
		LDVT = VT_t.get()->get_size(0);

	} else {
		JOBVT = 'N';
	}

	//Lets allocate some work storage
	std::vector<unsigned int> work_dim(1);

	int LWORK = 5*2*min_M_N+max_M_N;

	hoNDArray< T > WORK;
	work_dim[0] = LWORK;

    if (!WORK.create(&work_dim)) {
		std::cout << fname << ": Unable to create temporary WORK storage\n" << std::endl;
		return -1;
    }

	hoNDArray< Y > RWORK;
	work_dim[0] = 5*min_M_N;
	if (!RWORK.create(&work_dim)) {
		std::cout << fname << ": Unable to create temporary RWORK storage\n" << std::endl;
		return -1;
	}

	//Now we are finally ready to call the SVD
	int INFO;

	gesvd_wrapper(&JOBU, &JOBVT, &M, &N, A_ptr,
			&M, S_ptr, U_ptr, &LDU, VT_ptr,
			&LDVT, WORK.get_data_ptr(),
			&LWORK, RWORK.get_data_ptr(), &INFO);

	if (INFO != 0) {
		std::cout << fname << ": Call to gesvd failed, INFO = " << INFO << "\n" << std::endl;
		return -1;
	}

	std::vector<unsigned int> permute_order(2);
	permute_order[0] = 1;permute_order[1] = 0;

	if (U) {
		if (U_t.get()->permute(&permute_order,U) != 0) {
			std::cout << fname << ": Failed to permute result (U)\n" << std::endl;
			return -1;
		}
	}

	if (VT) {
		if (VT_t.get()->permute(&permute_order,VT) != 0) {
			std::cout << fname << ": Failed to permute result (VT)\n" << std::endl;
			return -1;
		}
	}


	return 0;
}

//Template instanciations
template EXPORTLINALG int hoNDArray_svd(hoNDArray< float > *A, hoNDArray< float > *U, hoNDArray< float > *S, hoNDArray< float > *VT);
template EXPORTLINALG int hoNDArray_svd(hoNDArray< double > *A, hoNDArray< double > *U, hoNDArray< double > *S, hoNDArray< double > *VT);
template EXPORTLINALG int hoNDArray_svd(hoNDArray< std::complex<float> > *A, hoNDArray< std::complex<float> > *U, hoNDArray< float > *S, hoNDArray< std::complex<float> > *VT);
template EXPORTLINALG int hoNDArray_svd(hoNDArray< std::complex<double> > *A, hoNDArray< std::complex<double> > *U, hoNDArray< double > *S, hoNDArray< std::complex<double> > *VT);

