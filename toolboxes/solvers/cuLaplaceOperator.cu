#include "cuLaplaceOperator.h"
#include "vector_td.h"
#include "vector_td_utilities.h"
#include "check_CUDA.h"

// Template Power function
template<unsigned int i, unsigned int j>
struct Pow
{
  enum { Value = i*Pow<i,j-1>::Value};
};

template <unsigned int i>
struct Pow<i,1>
{
  enum { Value = i};
};

template<class REAL, class T, unsigned int D> __global__ void
laplace_kernel( typename intd<D>::Type dims, T *in, T *out )
{  
  const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  if( idx < prod(dims) ){
    
    T val = T(0);
    typename intd<D>::Type coN;

    typename intd<D>::Type co = idx_to_co<D>(idx, dims);

    typename intd<D>::Type stride;

    for (int i = 0; i < D; i++) stride.vec[i]=0;

    for (int d1 = -1; d1 < 2; d1++){
      stride.vec[0] = d1;

      if (D > 1){
	for (int d2 = -1; d2 < 2; d2++){
	  stride.vec[1] = d2;
	  if (D > 2){
	    for (int d3 = -1; d3 < 2; d3++){
	      stride.vec[2] = d3;
	      coN = (co+dims+stride)%dims;
	      val -=  in[co_to_idx<D>(coN, dims)];
	    }
	  } else { 
	      coN = (co+dims+stride)%dims;
	      val += T(0) - in[co_to_idx<D>(coN, dims)];
	  }
	}
      } else {
	coN = (co+dims+stride)%dims;
	val -= in[co_to_idx<D>(coN, dims)];
      }
    }
    out[idx] = val+in[co_to_idx<D>(co, dims)]*((REAL) Pow<3,D>::Value);
  }
}

template< class T, unsigned int D> void
cuLaplaceOperator<T,D>::compute_laplace( cuNDArray<T> *in, cuNDArray<T> *out, bool accumulate )
{
  
  if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
  	throw std::runtime_error("laplaceOperator::compute_laplace : array dimensions mismatch.");

  }
  
  typename uintd<D>::Type _dims = vector_to_uintd<D>( *(in->get_dimensions().get()) );
  typename intd<D>::Type dims;
  for( unsigned int i=0; i<D; i++ ){
    dims.vec[i] = (int)_dims.vec[i];
  }  
  


  _set_device();
  
  dim3 dimBlock( dims.vec[0] );
  dim3 dimGrid( 1, dims.vec[D-1] );
  
  for( unsigned int d=1; d<D-1; d++ )
    dimGrid.x *= dims.vec[d];
  
  // Invoke kernel
  laplace_kernel<typename realType<T>::type ,T,D><<< dimGrid, dimBlock >>> (dims, in->get_data_ptr(), out->get_data_ptr() );
  
  CHECK_FOR_CUDA_ERROR();

  _restore_device();


}

// Instantiations

template class EXPORTSOLVERS cuLaplaceOperator<float, 1>;
template class EXPORTSOLVERS cuLaplaceOperator<float, 2>;
template class EXPORTSOLVERS cuLaplaceOperator<float, 3>;

template class EXPORTSOLVERS cuLaplaceOperator<float_complext, 1>;
template class EXPORTSOLVERS cuLaplaceOperator<float_complext, 2>;
template class EXPORTSOLVERS cuLaplaceOperator<float_complext, 3>;


template class EXPORTSOLVERS cuLaplaceOperator<double, 1>;
template class EXPORTSOLVERS cuLaplaceOperator<double, 2>;
template class EXPORTSOLVERS cuLaplaceOperator<double, 3>;

template class EXPORTSOLVERS cuLaplaceOperator<double_complext, 1>;
template class EXPORTSOLVERS cuLaplaceOperator<double_complext, 2>;
template class EXPORTSOLVERS cuLaplaceOperator<double_complext, 3>;

