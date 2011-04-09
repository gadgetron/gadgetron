#include "cuCG.h"

__global__ void clear_array_kernel(float* in, unsigned long int elements)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < elements) {
    in[idx_in] = 0.0;
  }
}

__global__ void clear_array_kernel(float2* in, unsigned long int elements)
{
  unsigned long idx_in = blockIdx.x*blockDim.x+threadIdx.x;
  if (idx_in < elements) {
    in[idx_in].x = 0.0;
    in[idx_in].y = 0.0;
  }
}

template <class T> int clear_cuNDArray(cuNDArray<T>* in)
{
  dim3 blockDim(512,1,1);
  dim3 gridDim((unsigned int) ceil((double)in->get_number_of_elements()/blockDim.x), 1, 1 );

  clear_array_kernel<<< gridDim, blockDim >>>( in->get_data_ptr(), in->get_number_of_elements());

  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    std::cerr << "clear_cuNDArray : Error during kernel call: " << cudaGetErrorString(err) << std::endl;
    return -1;
  }

  return 0;
}

float inner_product(cuNDArray<float2>* arr, cuNDArray<float2>* arr2, cublasHandle_t handle)
{
  float2 ret;
  ret.x = 0.0f;
  ret.y = 0.0f;

   if (cublasCdotc(handle, arr->get_number_of_elements(),
		   arr->get_data_ptr(), 1, 
		   arr2->get_data_ptr(), 1,
		   &ret) != CUBLAS_STATUS_SUCCESS) 
     {
       std::cerr << "cuCG inner product calculating using cublas failed" << std::endl;
       return 0.0f;
     }

  return ret.x;
}


float inner_product(cuNDArray<float>* arr, cuNDArray<float>* arr2, cublasHandle_t handle)
{

  float ret = 0.0f;
  if (cublasSdot(handle, arr->get_number_of_elements(),
		  arr->get_data_ptr(), 1, 
		  arr2->get_data_ptr(), 1,
		  &ret) != CUBLAS_STATUS_SUCCESS) 
    {
      std::cerr << "cuCG inner product calculating using cublas failed" << std::endl;
      return 0.0f;
    }
  return ret;
}

int axpy(float a, cuNDArray<float2>* x, cuNDArray<float2>* y, cublasHandle_t handle)
{
  float2 a_int;
  a_int.x = a;
  a_int.y = 0.0;

  if (x->get_number_of_elements() != y->get_number_of_elements()) {
      std::cerr << "cuCG axpy array dimensions mismatch" << std::endl;
      return -1;
  }
  
  if (cublasCaxpy(handle, x->get_number_of_elements(), &a_int,
		  x->get_data_ptr(), 1, 
		  y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS) 
    {
      std::cerr << "cuCG axpy calculating using cublas failed" << std::endl;
      return -2;
    }
  
  return 0;
} 

int axpy(float a, cuNDArray<float>* x, cuNDArray<float>* y, cublasHandle_t handle)
{
  if (x->get_number_of_elements() != y->get_number_of_elements()) {
      std::cerr << "cuCG axpy array dimensions mismatch" << std::endl;
      return -1;
  }
  
  if (cublasSaxpy(handle, x->get_number_of_elements(), &a,
		  x->get_data_ptr(), 1, 
		  y->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS) 
    {
      std::cerr << "cuCG axpy calculating using cublas failed" << std::endl;
      return -2;
    }
  
  return 0;
} 

int scal(float a, cuNDArray<float2>* x, cublasHandle_t handle) 
{
  float2 a_int;
  a_int.x = a;
  a_int.y = 0.0;

  if (cublasCscal(handle, x->get_number_of_elements(), &a_int,
		  x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS) 
    {
      std::cerr << "cuCG scal calculating using cublas failed" << std::endl;
      return -1;
    }

  return 0;
}


int scal(float a, cuNDArray<float>* x, cublasHandle_t handle) 
{
  if (cublasSscal(handle, x->get_number_of_elements(), &a,
		  x->get_data_ptr(), 1) != CUBLAS_STATUS_SUCCESS) 
    {
      std::cerr << "cuCG scal calculating using cublas failed" << std::endl;
      return -1;
    }

  return 0;
}

template <class T> cuNDArray<T> cuCG<T>::solve(cuNDArray<T>* rhs)
{
  cuNDArray<T> rho;
  if (!rho.create(rhs->get_dimensions())) {
    std::cerr << "cuCG<T>::solve : Unable to allocate temp storage (rho)" << std::endl;
    return rho;
  }
  if (clear_cuNDArray(&rho) < 0) {
    std::cerr << "cuCG<T>::solve : failed to clear rho" << std::endl;
    return rho;
  }

  //Calculate residual r
  cuNDArray<T> r;
  if (precond_) {
    if (!r.create(rhs->get_dimensions())) {
      std::cerr << "cuCG<T>::solve : Unable to allocate storage (r)" << std::endl;
      return rho;
    }
    if (precond_->apply(rhs,&r) < 0) {
      std::cerr << "cuCG<T>::solve : Unable to apply preconditioning to rhs" << std::endl;
      return rho;
    }
  } else {
    r =  *rhs;
  }


  double rr_0    = inner_product(&r, &r, cublas_handle_);
  double rr_1    = rr_0;
  double rr      = 0;
  double rr_last = 1e10;

  cuNDArray<T> p;
  if (!p.create(rhs->get_dimensions())) {
    std::cerr << "cuCG<T>::solve : Unable to allocate temp storage (p)" << std::endl;
    return rho;
  }

  cuNDArray<T> p_precond;
  if (precond_) { //We only need this additional storage if we are using a preconditioner
    if (!p_precond.create(rhs->get_dimensions())) {
      std::cerr << "cuCG<T>::solve : Unable to allocate temp storage (p_precond)" << std::endl;
      return rho;
    }
  }

  cuNDArray<T> q;
  if (!q.create(rhs->get_dimensions())) {
    std::cerr << "cuCG<T>::solve : Unable to allocate temp storage (q)" << std::endl;
    return rho;
  }

  cuNDArray<T> q2;
  if (!q2.create(rhs->get_dimensions())) {
    std::cerr << "cuCG<T>::solve : Unable to allocate temp storage (q2)" << std::endl;
    return rho;
  }

  double alpha, beta, rel_res;

  if (output_mode_ >= OUTPUT_VERBOSE) {
    std::cout << "Iterating..." << std::endl;
  }

  for (unsigned int it = 0; it < iterations_; it++) { //iterations_; it++) {
    rr_1 = rr;
    rr = inner_product(&r, &r, cublas_handle_);
    
    //Update p
    if (it == 0){
      p = r;
    } else {        
      beta = rr/rr_1;
      if (scal(beta,&p,cublas_handle_) < 0) {
	std::cerr << "cuCG<T>::solve : failed to scale p" << std::endl;
	return rho;
      }
      if (axpy(1.0,&r,&p,cublas_handle_) < 0) {
	std::cerr << "cuCG<T>::solve : failed to add r to scaled p" << std::endl;
	return rho;
      }
    }

    //Now we need to multiply with the system matrix
    if (clear_cuNDArray(&q) < 0) {
      std::cerr << "cuCG<T>::solve : failed to clear q" << std::endl;
      return rho;
    }
    

    //Take care of preconditioning
    cuNDArray<T>* cur_p = &p;
    if (precond_) {
      if (!precond_->apply(&p,&p_precond) < 0) {
	std::cerr << "cuCG<T>::solve : failed to apply preconditioner to p" << std::endl;
	return rho;
      }
      cur_p = &p_precond;
    }

    for (unsigned int i = 0; i < operators_.size(); i++) {
      if (operators_[i]->mult_MH_M(cur_p, &q2, false) < 0) {
	std::cerr << "cuCG<T>::solve : failed to apply operator number " << i << std::endl;
	return rho;
      }
      if (axpy(weights_[i],&q2,&q,cublas_handle_) < 0) {
	std::cerr << "cuCG<T>::solve : failed to add q1 to q" << std::endl;
	return rho;
      }
    }

    if (precond_) {
      if (!precond_->apply(&q,&q) < 0) {
	std::cerr << "cuCG<T>::solve : failed to apply preconditioner to q" << std::endl;
	return rho;
      }
    }

    alpha = rr/(inner_product(&p,&q,cublas_handle_));
    
    //Update solution
    if (axpy(alpha,&p,&rho,cublas_handle_) < 0) {
      std::cerr << "cuCG<T>::solve : failed to update solution" << std::endl;
      return rho;
    }
    
    //Update residual
    if (axpy(-alpha,&q,&r,cublas_handle_) < 0) {
      std::cerr << "cuCG<T>::solve : failed to update residual" << std::endl;
      return rho;
    }

    /* Calculate relative residual norm */
    rel_res = rr/rr_0;
    
    if (output_mode_ >= OUTPUT_WARNINGS) {
      if (output_mode_ >= OUTPUT_VERBOSE) {
	std::cout << "Iteration " << it+1 << ". rr/rr_0 = " << rel_res << std::endl;
      }
      if ((rr_last-rel_res) < 0.0) {
	std::cout << "----- Warning: CG residual increase. Stability problem! -----" << std::endl;
      }
    }

    if (rel_res < limit_) {
      break;
    } else {
      rr_last = rel_res;
    }

  }

  if (precond_) {
    if (!precond_->apply(&rho,&rho) < 0) {
      std::cerr << "cuCG<T>::solve : failed to apply preconditioner to rho" << std::endl;
      return rho;
    }
  }

  return rho;
}

template class cuCG<float>;
template class cuCG<float2>;
