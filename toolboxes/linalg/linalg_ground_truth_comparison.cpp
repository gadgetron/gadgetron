/*
* linalg_test.cpp
*
*  Created on: Dec 9, 2011
*      Author: Michael S. Hansen
*/

#include <iostream>
#include <hoNDArray.h>
#include <hoNDArray_fileio.h>
#include <hoNDArray_utils.h>
#include <matrix_vector_op.h>
#include <matrix_decomposition.h>
#include "GadgetronTimer.h"
#include "hoNDArray_math_util.h"
#include "hoNDArray_elemwise.h"
#include "hoMatrix_util.h"
#include "hoNDMath_util.h"
#include "hoNDFFT.h"
#include <fftw3.h>
#include <valarray>
#include <omp.h>

#define DIFF_LIMIT 1e-5

using namespace Gadgetron;

double mcompare(hoNDArray< std::complex<float> >* A, hoNDArray< std::complex<float> >* B)
{
    float comp = 0.0;
    float root_sum = 0.0;
    if (A->get_number_of_elements() != B->get_number_of_elements()) {
        std::cout << "Wrong number of elements in comparison" << std::endl;
        return 9e30;
    }

    hoNDArray< std::complex<float> > diff;
    Gadgetron::subtract(*A, *B, diff);

    /*std::complex<float>* a = A->get_data_ptr();
    std::complex<float>* b = B->get_data_ptr();
    for (unsigned long int i = 0; i < A->get_number_of_elements(); i++) {
    comp += abs(a[i]-b[i]);
    root_sum += abs(a[i]*conj(b[i]));
    }
    comp /= root_sum;*/

    Gadgetron::norm1(diff, comp);

    std::complex<float> r;
    Gadgetron::dotc(*A, *B, r);
    if ( std::abs(r) > 0 ) comp /= std::abs(r);
    return comp;
}

double mcompare(hoNDArray< float >* A, hoNDArray< float >* B)
{
    float comp = 0.0;
    float root_sum = 0.0;
    if (A->get_number_of_elements() != B->get_number_of_elements()) {
        std::cout << "Wrong number of elements in comparison" << std::endl;
        return 9e30;
    }

    hoNDArray< float > diff;
    Gadgetron::subtract(*A, *B, diff);

    /*std::complex<float>* a = A->get_data_ptr();
    std::complex<float>* b = B->get_data_ptr();
    for (unsigned long int i = 0; i < A->get_number_of_elements(); i++) {
    comp += abs(a[i]-b[i]);
    root_sum += abs(a[i]*conj(b[i]));
    }
    comp /= root_sum;*/

    Gadgetron::norm1(diff, comp);

    float r;
    Gadgetron::math::dotu(A->get_number_of_elements(), A->begin(), B->begin(), r);
    if ( std::abs(r) > 0 )  comp /= std::abs(r);
    return comp;
}

void compare_result(hoNDArray< std::complex<float> >& res, hoNDArray< std::complex<float> >& res_math, const std::string& msg)
{
    double diff = mcompare(&res, &res_math);
    if (diff > DIFF_LIMIT)
    {
        std::cout << msg << " - FAILED with diff: " << diff << std::endl;
    }
    else
    {
        std::cout << msg << " - SUCCESS with diff: " << diff << std::endl;
    }
}

void compare_result(hoNDArray< float >& res, hoNDArray< float >& res_math, const std::string& msg)
{
    double diff = mcompare(&res, &res_math);
    if (diff > DIFF_LIMIT)
    {
        std::cout << msg << " - FAILED with diff: " << diff << std::endl;
    }
    else
    {
        std::cout << msg << " - SUCCESS with diff: " << diff << std::endl;
    }
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

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("matrix multiplication");
    GADGET_MSG("------------------------------------------------------------------");

    std::complex<float> alpha(1.0,0);
    std::complex<float> beta(1.0,0);
    double diff;

    /*{
    GadgetronTimer t("GEMM Time (system)", true);
    hoNDArray_gemm( A.get(), B.get(), alpha,  C1.get(), beta);
    std::cout << C1->get_size(0) << ", " << C1->get_size(1) << ", " << C1->get_number_of_elements() << std::endl;
    }

    {
    GadgetronTimer t("GEMM Time (MKL)", true);
    GeneralMatrixProduct_gemm_CXFL( *C1.get(), *B.get(), *A.get());
    std::cout << C1->get_size(0) << ", " << C1->get_size(1) << ", " << C1->get_number_of_elements() << std::endl;
    }

    {
    GadgetronTimer t("Write time", true);
    write_nd_array< std::complex<float> >(C1.get(), "C2_calc.cplx");
    }

    double diff;
    {
    GadgetronTimer compare("CompareTime", true);
    diff = mcompare(C1.get(),C2.get());
    }

    if (diff > DIFF_LIMIT) {
    std::cout << "Complex GEMM FAILED with diff: " << diff << std::endl;
    return -1;
    } else {
    std::cout << "Complex GEMM SUCCESS with diff: " << diff << std::endl;
    }

    hoNDArray_choldc(S.get());
    zero_tril(S.get());

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
    }*/

    hoNDArray<std::complex<float> > a(*A);
    hoNDArray<std::complex<float> > b(*A);

    Gadgetron::scal( std::complex<float>(2), b);

    float r;

    hoNDArray<std::complex<float> > res, res_math;
    hoNDArray<float > res_f, res_f_math;

    {
        GadgetronTimer t("allocate res", true);
        res = a;
        res_math = a;

        res_f.create(a.get_dimensions());
        res_f_math.create(a.get_dimensions());
    }

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("axpy");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("axpy Time (MKL)", true);
        Gadgetron::axpy( alpha, a, b, res);
    }

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector add");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("vzAdd Time (MKL)", true);
        Gadgetron::add( *A.get(), *A.get(), res);
    }

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector subtract");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("vzSub Time (MKL)", true);
        Gadgetron::subtract( a, b, res);
    }


    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector multiplication");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("vzMul Time (MKL)", true);
        Gadgetron::multiply( a, b, res);
    }

    compare_result(res, res_math, "multiply");

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector addEpsilon");
    GADGET_MSG("------------------------------------------------------------------");

    res = a;
    res_math = a;

    {
        GadgetronTimer t("addEpsilon Time (MKL)", true);
        Gadgetron::addEpsilon( res );
    }

    compare_result(res, res_math, "addEpsilon");

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector divide");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("divide Time (MKL)", true);
        Gadgetron::divide( a, res, res);
    }

    compare_result(res, res_math, "divide");

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector sqrt");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("sqrt Time (MKL)", true);
        Gadgetron::sqrt( a, res);
    }

    compare_result(res, res_math, "sqrt");

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector conjugate");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("conjugate Time (MKL)", true);
        Gadgetron::conjugate( a, res);
    }

    compare_result(res, res_math, "conjugate");

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector conjugate multiplication");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("vcMulByConj Time (MKL)", true);
        Gadgetron::multiplyConj( a, b, res);
    }

    compare_result(res, res_math, "multiplyConj");

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector scal");
    GADGET_MSG("------------------------------------------------------------------");

    res = a;
    res_math = a;

    {
        GadgetronTimer t("scal Time (MKL)", true);
        Gadgetron::scal( alpha, a);
    }

    compare_result(res, res_math, "scal");

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector dotc");
    GADGET_MSG("------------------------------------------------------------------");

    std::complex<float> rdotc(0);

    {
        GadgetronTimer t("dotc Time (MKL)", true);
        rdotc = Gadgetron::dotc( a, b);
    }
    std::cout << "dotc = " << rdotc << std::endl;

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector dotu");
    GADGET_MSG("------------------------------------------------------------------");

    std::complex<float> rdotu;

    {
        GadgetronTimer t("dotu Time (MKL)", true);
        rdotu = Gadgetron::dotu( a, b );
    }
    std::cout << "dotu = " << rdotu << std::endl;

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector absolute");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("absolute Time (MKL)", true);
        Gadgetron::absolute( a, res);
    }

    compare_result(res, res_math, "absolute");

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector argument");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("argument Time (MKL)", true);
        Gadgetron::argument( a, res_f);
    }

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("vector inv");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("inv Time (MKL)", true);
        Gadgetron::inv( a, res);
    }

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("norm2");
    GADGET_MSG("------------------------------------------------------------------");

    float rn;

    {
        GadgetronTimer t("Time (MKL)", true);
        Gadgetron::norm2( a, rn);
    }
    std::cout << "nrm2 = " << rn << std::endl;

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("norm1");
    GADGET_MSG("------------------------------------------------------------------");

    {
        GadgetronTimer t("Time (MKL)", true);
        Gadgetron::norm1( a, rn);
    }
    std::cout << "nrm1 = " << rn << std::endl;

    GADGET_MSG("------------------------------------------------------------------");
    GADGET_MSG("conv2");
    GADGET_MSG("------------------------------------------------------------------");

    hoNDArray<std::complex<float> > ker;
    ker.create(3, 3);
    Gadgetron::fill(ker, std::complex<float>(1) );

    {
        GadgetronTimer t("conv2 Time (MKL)", true);
        Gadgetron::conv2( a, ker, res);
    }
}
