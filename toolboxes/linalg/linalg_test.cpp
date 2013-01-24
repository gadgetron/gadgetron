/*
 * linalg_test.cpp
 *
 *  Created on: Dec 9, 2011
 *      Author: Michael S. Hansen
 */

#include <iostream>
#include <hoNDArray.h>
#include <hoNDArray_fileio.h>
#include <matrix_vector_op.h>
#include <matrix_decomposition.h>

#define DIFF_LIMIT 1e-6

double mcompare(hoNDArray< std::complex<float> >* A, hoNDArray< std::complex<float> >* B)
{
	double comp = 0.0;
	double root_sum = 0.0;
	if (A->get_number_of_elements() != B->get_number_of_elements()) {
		std::cout << "Wrong number of elements in comparison" << std::endl;
		return 9e30;
	}

	std::complex<float>* a = A->get_data_ptr();
	std::complex<float>* b = B->get_data_ptr();
	for (unsigned long int i = 0; i < A->get_number_of_elements(); i++) {
		comp += abs(a[i]-b[i]);
		root_sum += abs(a[i]*conj(b[i]));
	}
	comp /= root_sum;
	return comp;
}

/*
 *   Simple test program for linear algebra routines.
 */
int main(int argc, char** argv)
{
	std::cout << "Simple test of linear algebra routines" << std::endl;
	if (argc != 2) {
		std::cout << "Usage: linalg_test <folder_with_test_data>" << std::endl;
		return -1;
	}

	std::string filenameA = std::string(argv[1]) + std::string("/A.cplx");
	std::string filenameB = std::string(argv[1]) + std::string("/B.cplx");
	std::string filenameC1 = std::string(argv[1]) + std::string("/C1.cplx");
	std::string filenameC2 = std::string(argv[1]) + std::string("/C2.cplx");
	std::string filenameS = std::string(argv[1]) + std::string("/S.cplx");
	std::string filenameS_chol = std::string(argv[1]) + std::string("/S_chol.cplx");
	std::string filenameS_chol_inv = std::string(argv[1]) + std::string("/S_chol_inv.cplx");

	boost::shared_ptr< hoNDArray<std::complex<float> > > A = read_nd_array< std::complex<float> >(filenameA.c_str());
	boost::shared_ptr< hoNDArray<std::complex<float> > > B = read_nd_array< std::complex<float> >(filenameB.c_str());
	boost::shared_ptr< hoNDArray<std::complex<float> > > C1 = read_nd_array< std::complex<float> >(filenameC1.c_str());
	boost::shared_ptr< hoNDArray<std::complex<float> > > C2 = read_nd_array< std::complex<float> >(filenameC2.c_str());

	boost::shared_ptr< hoNDArray<std::complex<float> > > S = read_nd_array< std::complex<float> >(filenameS.c_str());
	boost::shared_ptr< hoNDArray<std::complex<float> > > S_chol = read_nd_array< std::complex<float> >(filenameS_chol.c_str());
	boost::shared_ptr< hoNDArray<std::complex<float> > > S_chol_inv = read_nd_array< std::complex<float> >(filenameS_chol_inv.c_str());

	std::complex<float> alpha(1.0,0);
	std::complex<float> beta(1.0,0);
	hoNDArray_gemm( A.get(), B.get(), alpha,  C1.get(), beta);

	write_nd_array< std::complex<float> >(C1.get(), "C2_calc.cplx");

	double diff = mcompare(C1.get(),C2.get());
	if (diff > DIFF_LIMIT) {
		std::cout << "Complex GEMM FAILED with diff: " << diff << std::endl;
		return -1;
	} else {
		std::cout << "Complex GEMM SUCCESS with diff: " << diff << std::endl;
	}


	hoNDArray_choldc(S.get());

	write_nd_array< std::complex<float> >(S.get(), "S_chol_calc.cplx");

	diff = mcompare(S.get(),S_chol.get());
	if (diff > DIFF_LIMIT) {
		std::cout << "Complex Cholesky decomposition FAILED with diff: " << diff << std::endl;
		return -1;
	} else {
		std::cout << "Complex Cholesky decomposition SUCCESS with diff: " << diff << std::endl;
	}

	hoNDArray_inv_lower_triangular(S.get());

	write_nd_array< std::complex<float> >(S.get(), "S_chol_inv_calc.cplx");

	diff = mcompare(S.get(),S_chol_inv.get());
	if (diff > DIFF_LIMIT) {
		std::cout << "Complex Triangular inversion FAILED with diff: " << diff << std::endl;
		return -1;
	} else {
		std::cout << "Complex Triangular inversion SUCCESS with diff: " << diff << std::endl;
	}

	return 0;
}




