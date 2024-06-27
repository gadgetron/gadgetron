#include "cuConvolutionOperator.h"
#include "vector_td_utilities.h"
#include "cudaDeviceManager.h"
#include "setup_grid.h"

namespace Gadgetron {

  // Mirror, but keep the origin unchanged
  template<class T, unsigned int D> __global__ void
  origin_mirror_kernel( vector_td<unsigned int,D> matrix_size, vector_td<unsigned int,D> origin, const T * __restrict__ in, T * __restrict__ out, bool zero_fill )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    
    if( idx < prod(matrix_size) ){
      
      vector_td<unsigned int,D> in_co = idx_to_co( idx, matrix_size );
      vector_td<unsigned int,D> out_co = matrix_size-in_co;
    
      bool wrap = false;
      for( unsigned int d=0; d<D; d++ ){
	if( out_co.vec[d] == matrix_size.vec[d] ){
	  out_co.vec[d] = 0;
	  wrap = true;
	}
      }
    
      const unsigned int in_idx = co_to_idx(in_co, matrix_size);
      const unsigned int out_idx = co_to_idx(out_co, matrix_size);

      if( wrap && zero_fill )
	out[out_idx] = T(0);
      else
	out[out_idx] = in[in_idx];
    }
  }
  
  // Mirror around the origin -- !! leaving the origin unchanged !!
  // This creates empty space "on the left" that can be filled by zero (default) or the left-over entry.
  template<class REAL, unsigned int D> void
  cuConvolutionOperator<REAL,D>::origin_mirror( cuNDArray< complext<REAL> > *in, cuNDArray< complext<REAL> > *out )
  {
    if( in == 0x0 || out == 0x0 ){
      throw std::runtime_error( "origin_mirror: 0x0 ndarray provided");
    }
    
    if( !in->dimensions_equal(out) ){
      throw std::runtime_error("origin_mirror: image dimensions mismatch");
    }
    
    if( in->get_number_of_dimensions() != D ){
      std::stringstream ss;
      ss << "origin_mirror: number of image dimensions is not " << D;
      throw std::runtime_error(ss.str());
    }

    typename uint64d<D>::Type matrix_size = from_std_vector<size_t,D>( in->get_dimensions() );
  
    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size), &blockDim, &gridDim );

    // Invoke kernel
    origin_mirror_kernel<complext<REAL>,D><<< gridDim, blockDim >>> 
      ( vector_td<unsigned int,D>(matrix_size), vector_td<unsigned int,D>(matrix_size>>1), in->get_data_ptr(), out->get_data_ptr(), true );
    
    CHECK_FOR_CUDA_ERROR();
  }


  template <class REAL, unsigned int D> void 
  cuConvolutionOperator<REAL,D>::operator_fft( bool forwards_transform, cuNDArray< complext<REAL> > *image )
  {
    if( forwards_transform )
      cuNDFFT<REAL>::instance()->fft(image);
    else
      cuNDFFT<REAL>::instance()->ifft(image);
  }    
  
  template EXPORTGPUOPERATORS class cuConvolutionOperator<float,1>;
  template EXPORTGPUOPERATORS class cuConvolutionOperator<float,2>;
  template EXPORTGPUOPERATORS class cuConvolutionOperator<float,3>;
  template EXPORTGPUOPERATORS class cuConvolutionOperator<float,4>;

  template EXPORTGPUOPERATORS class cuConvolutionOperator<double,1>;
  template EXPORTGPUOPERATORS class cuConvolutionOperator<double,2>;
  template EXPORTGPUOPERATORS class cuConvolutionOperator<double,3>;
  template EXPORTGPUOPERATORS class cuConvolutionOperator<double,4>;
  
}
