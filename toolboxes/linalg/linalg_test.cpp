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
#include "hoNDArray_blas.h"
#include "hoNDArray_operators.h"
#include "hoNDArray_elemwise.h"
#include "hoMatrix_util.h"
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
  comp /= std::abs(r);
  return comp;
}

/**
   Zero out everything except upper triangualar 
*/
void zero_tril(hoNDArray< std::complex<float> >* A)
{
  size_t rows = A->get_size(0);
  size_t cols = A->get_size(1);
  std::complex<float>* d = A->get_data_ptr();
  for (size_t c = 0; c < cols; c++) {
    for (size_t r = 0; r < c; r++) {
      d[r*cols+c] = std::complex<float>(0.0,0.0);
    } 
  }
}

template <typename T> 
bool multiplyOwn(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
{
  try
    {
      GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
      if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
	  r = x;
        }

      long long N = (long long)x.get_number_of_elements();
      long long n;

      const T* pX = x.begin();
      const T* pY = y.begin();
      T* pR = r.begin();

      if ( pR == pX )
        {
#pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
	  for ( n=0; n<(long long)N; n++ )
            {
	      pR[n] *= pY[n];
            }
        }
      else if ( pR == pY )
        {
#pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
	  for ( n=0; n<(long long)N; n++ )
            {
	      pR[n] *= pX[n];
            }
        }
      else
        {
#pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
	  for ( n=0; n<(long long)N; n++ )
            {
	      pR[n] = pX[n] * pY[n];
            }
        }
    }
  catch(...)
    {
      GADGET_ERROR_MSG("Error happened in multiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
      return false;
    }

  return true;
}

template <typename T> void vecMult(size_t N, T* a, T* b, T*c)
{
  size_t i;
#pragma omp parallel for private(i)// shared(num, p, a, n0, n1, r)
  for (i = 0; i < N; i++) {
    c[i] = a[i]*b[i];
  }
}

template <typename T> double hoNDArray_norm1(hoNDArray< std::complex<T> > * a)
{
  size_t N = a->get_number_of_elements();
  hoNDArray<T> ab;
  ab.create(N); 
  std::complex<T>* a_ptr = a->get_data_ptr();
  T* ab_ptr = ab.get_data_ptr();

  size_t i;
#pragma omp parallel for private(i)// shared(num, p, a, n0, n1, r)
  for (i = 0; i < N; i++) {
    ab_ptr[i] = abs(a_ptr[i]);
  }
  
  return hoNDArray_asum(&ab);

}

template<typename T> 
bool fftw_fft2(hoNDArray< std::complex<T> >& a, hoNDArray< std::complex<T> >& r, bool forward)
{
  r = a;

  int n0 = a.get_size(1);
  int n1 = a.get_size(0);

  T fftRatio = 1.0/std::sqrt( T(n0*n1) );

  size_t num = a.get_number_of_elements()/(n0*n1);
  std::cout << "Number of FFTs: " << num << std::endl;
  long long n;

  if ( typeid(T) == typeid(float) )
    {
      fftwf_init_threads();
      fftwf_plan_with_nthreads(omp_get_max_threads());

      fftwf_plan p;

      {
	GadgetronTimer tp("FFTW Planning", true);
	// mutex_.lock();
	if ( forward )
	  {
	    p = fftwf_plan_dft_2d(n0, n1,
				  reinterpret_cast<fftwf_complex*>(a.begin()), 
				  reinterpret_cast<fftwf_complex*>(r.begin()),
				  FFTW_FORWARD, FFTW_ESTIMATE);
	  }
	else
	  {
	    p = fftwf_plan_dft_2d(n0, n1,
				  reinterpret_cast<fftwf_complex*>(a.begin()), 
				  reinterpret_cast<fftwf_complex*>(r.begin()),
				  FFTW_BACKWARD, FFTW_ESTIMATE);
	  }
	// mutex_.unlock();
      }

      {
	GadgetronTimer t("FFT loop time", true); 
	//#pragma omp parallel for private(n) shared(num, p, a, n0, n1, r)
	  for ( n=0; n<num; n++ )
	    {
	      fftwf_execute_dft(p, reinterpret_cast<fftwf_complex*>(a.begin()+n*n0*n1), 
				reinterpret_cast<fftwf_complex*>(r.begin()+n*n0*n1));
	    }
      }

      {
	// mutex_.lock();
	fftwf_destroy_plan(p);
	// mutex_.unlock();
      }
      fftwf_cleanup_threads();
    }
  else if ( typeid(T) == typeid(double) )
    {
      fftw_init_threads();
      fftw_plan_with_nthreads(omp_get_max_threads());
      fftw_plan p;

      {
	// mutex_.lock();
	if ( forward )
	  {
	    p = fftw_plan_dft_2d(n0, n1,
				 reinterpret_cast<fftw_complex*>(a.begin()), 
				 reinterpret_cast<fftw_complex*>(r.begin()),
				 FFTW_FORWARD, FFTW_ESTIMATE);
	  }
	else
	  {
	    p = fftw_plan_dft_2d(n0, n1,
				 reinterpret_cast<fftw_complex*>(a.begin()), 
				 reinterpret_cast<fftw_complex*>(r.begin()),
				 FFTW_BACKWARD, FFTW_ESTIMATE);
	  }
	// mutex_.unlock();
      }

      //z#pragma omp parallel for private(n) shared(num, p, a, n0, n1, r)
        for ( n=0; n<num; n++ )
	  {
            fftw_execute_dft(p, reinterpret_cast<fftw_complex*>(a.begin()+n*n0*n1), 
			     reinterpret_cast<fftw_complex*>(r.begin()+n*n0*n1));
	  }

        {
	  // mutex_.lock();
	  fftw_destroy_plan(p);
	  fftw_cleanup_threads();
	  // mutex_.unlock();
        }
    }

  {
    GadgetronTimer tt("FFT Scaling", true);
    Gadgetron::hoNDArray_scal(std::complex<float>(fftRatio,0.0), &r);
    //r *= fftRatio;
  }
  return true;
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

  {
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
  }

  GADGET_MSG("------------------------------------------------------------------");
  GADGET_MSG("vector add");
  GADGET_MSG("------------------------------------------------------------------");

  hoNDArray<std::complex<float> > a(*A);
  hoNDArray<std::complex<float> > b(*A);

  hoNDArray<std::complex<float> > res;

  {
    GadgetronTimer t("allocate res", true);
    res = a;
  }

  {
    GadgetronTimer t("axpy Time (system)", true);
    Gadgetron::axpy( alpha, &a, &res);
  }

  {
    GadgetronTimer t("linalg (direct BLAS) axpy Time (system)", true);
    Gadgetron::hoNDArray_axpy( &alpha, &a, &res);
  }

  {
    GadgetronTimer t("operator +", true);
    res += a;
  }

  {
    GadgetronTimer t("vzAdd Time (MKL)", true);
    Gadgetron::add( *A.get(), *A.get(), res);
  }

  GADGET_MSG("------------------------------------------------------------------");
  GADGET_MSG("vector multiplication");
  GADGET_MSG("------------------------------------------------------------------");

  {
    GadgetronTimer t("operator *", true);
    res *= a;
  }

  {
    GadgetronTimer t("vzMul Time (MKL)", true);
    Gadgetron::multiply( a, b, res);
  }

  {
    GadgetronTimer t("multiplyOwn Time (openMP)", true);
    multiplyOwn( a, b, res);
  }

  {
    GadgetronTimer t("vecMult (vectorized)", true);
    vecMult( a.get_number_of_elements(), a.get_data_ptr(), b.get_data_ptr(), res.get_data_ptr());
  }

  GADGET_MSG("------------------------------------------------------------------");
  GADGET_MSG("norm2");
  GADGET_MSG("------------------------------------------------------------------");

  float rn;

  {
    GadgetronTimer t("nrm2", true);
    rn = Gadgetron::nrm2(&a);
  }
  std::cout << "nrm2 = " << rn << std::endl;
	
  {
    GadgetronTimer t("Time (MKL)", true);
    Gadgetron::norm2( a, rn);
  }
  std::cout << "nrm2 = " << rn << std::endl;

  {
    GadgetronTimer t("Time (DIRECT BLAS)", true);
    rn = Gadgetron::hoNDArray_norm2(&a);
  }
  std::cout << "nrm2 = " << rn << std::endl;

  GADGET_MSG("------------------------------------------------------------------");
  GADGET_MSG("norm1");
  GADGET_MSG("------------------------------------------------------------------");

  {
    GadgetronTimer t("nrm1", true);
    rn = Gadgetron::nrm1(&a);
  }
  std::cout << "nrm1 = " << rn << std::endl;

  {
    GadgetronTimer t("Time (MKL)", true);
    Gadgetron::norm1( a, rn);
  }
  std::cout << "nrm1 = " << rn << std::endl;

  {
    GadgetronTimer t("Time (DIRECT BLAS)", true);
    rn = hoNDArray_norm1(&a);
  }
  std::cout << "nrm1 = " << rn << std::endl;

  GADGET_MSG("------------------------------------------------------------------");
  GADGET_MSG("axpy");
  GADGET_MSG("------------------------------------------------------------------");

  {
    GadgetronTimer t("axpy Time (system)", true);
    Gadgetron::axpy( alpha, &a, &res);
  }

  {
    GadgetronTimer t("linalg (direct BLAS) axpy Time (system)", true);
    Gadgetron::hoNDArray_axpy( &alpha, &a, &res);
  }

  {
    GadgetronTimer t("axpy Time (MKL)", true);
    Gadgetron::axpy( alpha, a, b, res);
  }

  GADGET_MSG("------------------------------------------------------------------");
  GADGET_MSG("fft 2D");
  GADGET_MSG("------------------------------------------------------------------");

  {
    GadgetronTimer t("fft2 (MKL)", true);
    hoNDFFT<float>::instance()->fft2(a, res);
  }
  Gadgetron::norm2(res, rn); GADGET_MSG("rn = " << rn);

  {
    GadgetronTimer t("fftw_fft2", true);
    fftw_fft2(a, res, true);
  }
  Gadgetron::norm2(res, rn); GADGET_MSG("rn = " << rn);

  return 0;
}




