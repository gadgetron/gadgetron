#include "htgrappa.h"
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_linalg.h"
#include "cuNDArray.h"
#include "CUBLASContextProvider.h"

#ifdef USE_OMP
#include <omp.h>
#endif


/*
  This file is used to hide certain Armadillo calls from the nvcc compiler. If Armadillo functions need to
  be called in a *.cu file, it is preferably to wrap the calls in a function and place that function in
  a *.cpp file so that Armadillo code will not be compiled by nvcc.

  Some error handling may be needed in these functions, but eventually SymmetricHermitianPositiveDefiniteLinearSystem_posv
  will be renamed and made to throw exceptions and then it should be handled.

 */



namespace Gadgetron{


template <class T> void ht_grappa_solve_spd_system(cuNDArray<T>& A, cuNDArray<T>& B) {
    boost::shared_ptr< hoNDArray<T> > A_h = A.to_host();
    boost::shared_ptr< hoNDArray<T> > B_h = B.to_host();

    std::vector<size_t> perm_dim ={1,0};

    permute(*A_h,perm_dim);
    permute(*B_h,perm_dim);

    ht_grappa_solve_spd_system(A_h.get(), B_h.get());

    permute(*B_h,perm_dim);
    B = cuNDArray<T>(*B_h);

}


  template <class T> void ht_grappa_solve_spd_system(hoNDArray<T> *A, hoNDArray<T> *B) {
    /*
      We are swithcing off OpenMP threading before this call to posv. There seems to be a bad interaction between openmp, cuda, and BLAS.
      So far this problem has only been observed from *.cu files (or in functions called from *.cu files) but the problem may be more general.

      This is a temporary fix that we should keep an eye on.
     */

    hoNDArray<T> A_ori;
    A_ori = *A;

    try
    {
        posv(*A, *B);
    }
    catch(...)
    {
        // it is found that if signal is very very high, the posv can throw exceptions due to ill-conditioned matrix of A
        // hesv does not require A to be a positive-definite matrix, but an n-by-n symmetric matrix

        GERROR_STREAM("ht_grappa_solve_spd_system : posv(*A, *B) throws exceptions ... ");
        *A = A_ori;
        hesv(*A, *B);
        GERROR_STREAM("ht_grappa_solve_spd_system : hesv(*A, *B) is called ");
    }
  }


    template <class T> int inverse_clib_matrix(cuNDArray<T>* A,
                                           cuNDArray<T>* b,
                                           cuNDArray<T>* coeff,
                                           double lamda)
{
    // A: M*N
    // b: M*K
    size_t M = A->get_size(0);
    size_t N = A->get_size(1);

    size_t K = b->get_size(1);

    std::vector<size_t> AHA_dims(2,N);
    cuNDArray<T> AHA = cuNDArray<T>(AHA_dims);

    std::vector<size_t> AHrhs_dims;
    AHrhs_dims.push_back(N);
    AHrhs_dims.push_back(K);

    coeff->create(AHrhs_dims);

    cublasHandle_t handle = *CUBLASContextProvider::instance()->getCublasHandle();

    complext<float>  alpha = complext<float>(1);
    complext<float>  beta = complext<float>(0);

    //{
    //    std::string filename = debugFolder+"A.cplx";
    //    write_cuNDArray_to_disk(A, filename.c_str());
    //}

    //{
    //    std::string filename = debugFolder+"b.cplx";
    //    write_cuNDArray_to_disk(b, filename.c_str());
    //}

    {
        //GPUTimer t2("compute AHA ...");
        cublasStatus_t stat = cublasCgemm(handle, CUBLAS_OP_C, CUBLAS_OP_N,
                                          N,N,M,(float2*) &alpha,
                                          (float2*) A->get_data_ptr(), M,
                                          (float2*) A->get_data_ptr(), M,
                                          (float2*) &beta, (float2*) AHA.get_data_ptr(), N);

        if (stat != CUBLAS_STATUS_SUCCESS)
        {
            std::cerr << "inverse_clib_matrix: Failed to form AHA product using cublas gemm" << std::endl;
            std::cerr << "---- cublas error code " << stat << std::endl;
            return -1;
        }
    }

    //{
    //    std::string filename = debugFolder+"AHA.cplx";
    //    write_cuNDArray_to_disk(&AHA, filename.c_str());
    //}

    {
        //GPUTimer t2("compute AHrhs ...");
        cublasStatus_t stat = cublasCgemm(handle, CUBLAS_OP_C, CUBLAS_OP_N,
                                          N,K,M,(float2*) &alpha,
                                          (float2*) A->get_data_ptr(), M,
                                          (float2*) b->get_data_ptr(), M,
                                          (float2*) &beta, (float2*)coeff->get_data_ptr(), N);

        if (stat != CUBLAS_STATUS_SUCCESS)
        {
            std::cerr << "inverse_clib_matrix: Failed to form AHrhs product using cublas gemm" << std::endl;
            std::cerr << "---- cublas error code " << stat << std::endl;
            return -1;
        }
    }

    //{
    //    std::string filename = debugFolder+"AHrhs.cplx";
    //    write_cuNDArray_to_disk(coeff, filename.c_str());
    //}

    // apply the regularization
    if ( lamda > 0 )
    {
        hoNDArray<T> AHA_host(N, N);
        float2* pAHA = (float2*) AHA_host.get_data_ptr();

        //GadgetronTimer timer;

        //timer.start("copy AHA to host");
        if (cudaMemcpy(pAHA, AHA.get_data_ptr(), AHA_host.get_number_of_bytes(), cudaMemcpyDeviceToHost) != cudaSuccess)
        {
            std::cerr << "inverse_clib_matrix: Failed to copy AHA to host" << std::endl;
            return -1;
        }
        //timer.stop();

        //timer.start("apply the regularization");
        // apply the regularization
        double trA = std::sqrt(pAHA[0].x*pAHA[0].x + pAHA[0].y*pAHA[0].y);
        size_t c;
        for ( c=1; c<N; c++ )
        {
            float x = pAHA[c+c*N].x;
            float y = pAHA[c+c*N].y;
            trA += std::sqrt(x*x+y*y);
        }

        double value = trA*lamda/N;
        for ( c=0; c<N; c++ )
        {
            float x = pAHA[c+c*N].x;
            float y = pAHA[c+c*N].y;
            pAHA[c+c*N].x = std::sqrt(x*x+y*y) + value;
            pAHA[c+c*N].y = 0;
        }
        //timer.stop();

        //timer.start("copy the AHA to device");
        if (cudaMemcpy(AHA.get_data_ptr(), pAHA, AHA_host.get_number_of_bytes(), cudaMemcpyHostToDevice) != cudaSuccess)
        {
            std::cerr << "inverse_clib_matrix: Failed to copy regularized AHA to device" << std::endl;
            return -1;
        }
        //timer.stop();
    }

    /*
      culaStatus s;
      s = culaDeviceCgels( 'N', N, N, K,
      (culaDeviceFloatComplex*)AHA.get_data_ptr(), N,
      (culaDeviceFloatComplex*)coeff->get_data_ptr(), N);
    */
    {
        //It actually turns out to be faster to do this inversion on the CPU. Problem is probably too small for GPU to make sense
        //GPUTimer cpu_invert_time("CPU Inversion time");
        boost::shared_ptr< hoNDArray<T> > AHA_h = AHA.to_host();
        boost::shared_ptr< hoNDArray<T> > AHrhs_h = coeff->to_host();

        std::vector<size_t> perm_dim;
        perm_dim.push_back(1);
        perm_dim.push_back(0);

        permute(*AHA_h,perm_dim);
        permute(*AHrhs_h,perm_dim);

        ht_grappa_solve_spd_system(AHA_h.get(), AHrhs_h.get());

        permute(*AHrhs_h,perm_dim);
        *coeff = cuNDArray<T>(*AHrhs_h);
    }


    //{
    //    std::string filename = debugFolder+"coeff.cplx";
    //    write_cuNDArray_to_disk(coeff, filename.c_str());
    //}

    /*
    if (s != culaNoError)
      {
        GDEBUG_STREAM("inverse_clib_matrix: linear solve failed" << std::endl);
        return -1;
      }
    */
    return 0;
}


template EXPORTGPUPMRI int inverse_clib_matrix(cuNDArray<complext<float> >* A,
                                               cuNDArray<complext<float> >* b,
                                               cuNDArray<complext<float> >* coeff,
                                               double lamda);
template void ht_grappa_solve_spd_system< float_complext >(hoNDArray< float_complext > *A, hoNDArray< float_complext > *B);
template void ht_grappa_solve_spd_system< float_complext >(cuNDArray< float_complext >& A, cuNDArray< float_complext >& B);

}
