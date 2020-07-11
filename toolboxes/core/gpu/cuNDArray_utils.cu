#include "cuNDArray_utils.h"
#include "vector_td_utilities.h"
#include "cudaDeviceManager.h"
#include "setup_grid.h"
#include "cuNDArray_math.h"

#include <cuda_runtime.h>
#include <cmath>

namespace Gadgetron {

namespace {
  template <class T> 
  __global__ void cuNDArray_permute_kernel(const  T*  __restrict__ in, T* __restrict__ out,
                                            unsigned int ndim,
                                            const unsigned int* __restrict__ dims,
                                            const unsigned int* __restrict__ strides_out,
                                            unsigned int elements)
  {
    unsigned int idx_in = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    unsigned int idx_out = 0;
    unsigned int idx_in_tmp = idx_in;

    if (idx_in < elements) {

      unsigned int cur_index;
      for (unsigned int i = 0; i < ndim; i++) {
        unsigned int idx_in_remainder = idx_in_tmp / dims[i];
        cur_index = idx_in_tmp-(idx_in_remainder*dims[i]); //cur_index = idx_in_tmp%dims[i];
        idx_out += cur_index*strides_out[i];
        idx_in_tmp = idx_in_remainder;
      }
      out[idx_in] = in[idx_out];
    }
  }

  template <class T> void cuNDArray_permute( const cuNDArray<T>& in,
                                             cuNDArray<T>& out,
                                             const std::vector<size_t>& order)
  {    


    std::vector<unsigned int> dims(in.get_number_of_dimensions());
    std::vector<unsigned int> strides_out(in.get_number_of_dimensions());


    for (unsigned int i = 0; i < in.get_number_of_dimensions(); i++) {
      dims[i] = in.get_size(order[i]);
      strides_out[i] = 1;    
      for (unsigned int j = 0; j < order[i]; j++) {
        strides_out[i] *= in.get_size(j);
      }
    }



    thrust::device_vector<unsigned int> strides_out_dev(strides_out);
    thrust::device_vector<unsigned int> dims_dev(dims);

    dim3 blockDim(512,1,1);
    dim3 gridDim;
    if( in.get_number_of_dimensions() > 2 ){
      gridDim = dim3((unsigned int) std::ceil((double)in.get_size(0)*in.get_size(1)/blockDim.x), 1, 1 );
      for( unsigned int d=2; d<in.get_number_of_dimensions(); d++ )
        gridDim.y *= in.get_size(d);
    }  else {
      gridDim = dim3((unsigned int) std::ceil((double) in.get_number_of_elements() / blockDim.x), 1, 1);
    }

    cuNDArray_permute_kernel<<< gridDim, blockDim >>>( in.get_data_ptr(),out.get_data_ptr(), in.get_number_of_dimensions(),
                                                       thrust::raw_pointer_cast(dims_dev.data()),
                                                       thrust::raw_pointer_cast(strides_out_dev.data()),
                                                       in.get_number_of_elements());

  }
}

  template <class T> cuNDArray<T>
  permute( const cuNDArray<T>&in, const std::vector<size_t>& dim_order)
  {

    std::vector<size_t> dims;
    for (size_t i = 0; i < dim_order.size(); i++)
      dims.push_back(in.get_size(dim_order[i]));

    auto  out = cuNDArray<T>(dims);
    permute( in, out, dim_order);
    return out;
  }

  template <class T> void
  permute(const cuNDArray<T>& in, cuNDArray<T>&out, const std::vector<size_t>& dim_order)
  {

    if (out.get_number_of_dimensions() != in.get_number_of_dimensions() || out.get_number_of_elements() != in.get_number_of_elements()){
    	throw std::runtime_error("permute(): Input and output have differing dimensions and/or differing number of elements");
    }


    //Check ordering array
    if (dim_order.size() > in.get_number_of_dimensions()) {
      throw std::runtime_error("permute(): invalid length of dimension ordering array");
    }

    std::vector<size_t> dim_count(in.get_number_of_dimensions(),0);
    for (unsigned int i = 0; i < dim_order.size(); i++) {
      if (dim_order[i] >= in.get_number_of_dimensions()) {
        throw std::runtime_error("permute(): invalid dimension order array");
      }
      dim_count[dim_order[i]]++;
    }

    //Create an internal array to store the dimensions
    std::vector<size_t> dim_order_int;

    //Check that there are no duplicate dimensions
    for (unsigned int i = 0; i < dim_order.size(); i++) {
      if (dim_count[dim_order[i]] != 1) {
        throw std::runtime_error("permute(): invalid dimension order array (duplicates)");
      }
      dim_order_int.push_back(dim_order[i]);
    }

    for (unsigned int i = 0; i < dim_order_int.size(); i++) {
      if (in.get_size(dim_order_int[i]) != out.get_size(i)) {
        throw std::runtime_error("permute(): dimensions of output array do not match the input array");
      }
    }

    //Pad dimension order array with dimension not mentioned in order array
    if (dim_order_int.size() < in.get_number_of_dimensions()) {
      for (unsigned int i = 0; i < dim_count.size(); i++) {
        if (dim_count[i] == 0) {
          dim_order_int.push_back(i);
        }
      }
    }


    //Check if permute is needed
    {
    	bool skip_permute = true;
    	for (size_t i = 0; i < dim_order_int.size(); i++)
    		skip_permute &= (i == dim_order_int[i]);

    	if (skip_permute){
    		out = in;
    		return;
    	}
    }
    cuNDArray_permute(in, out, dim_order_int);
  }

  template<class T> cuNDArray<T>
  shift_dim(const  cuNDArray<T> &in, int shift )
  {
    std::vector<size_t> order;
    for (int i = 0; i < in.get_number_of_dimensions(); i++) {
      order.push_back((i+shift)%in.get_number_of_dimensions());
    }
    return permute(in, order);
  }

  template<class T> 
  void shift_dim( const cuNDArray<T>&in, cuNDArray<T>&out, int shift )
  {

    std::vector<size_t> order;
    for (int i = 0; i < in.get_number_of_dimensions(); i++) {
      order.push_back((i+shift)%in.get_number_of_dimensions());
    }
    permute(in,out,order);
  }

  // Expand
  //
  template<class T> 
  __global__ void expand_kernel( 
                                const T * __restrict__ in, T * __restrict__ out,
                                unsigned int number_of_elements_in, unsigned int number_of_elements_out)
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;    
    if( idx < number_of_elements_out ){
      out[idx] = in[idx%number_of_elements_in];
    }
  }

  // Expand
  //
  template<class T> cuNDArray<T>
  expand( const cuNDArray<T>&in, size_t new_dim_size )
  {
    unsigned int number_of_elements_out = in.get_number_of_elements()*new_dim_size;

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( number_of_elements_out, &blockDim, &gridDim );

    // Find element stride
    std::vector<size_t> dims = *in.get_dimensions();
    dims.push_back(new_dim_size);

    // Invoke kernel
    cuNDArray<T> out(dims);

    expand_kernel<T><<< gridDim, blockDim >>>( in.get_data_ptr(), out.get_data_ptr(),
                                               in.get_number_of_elements(), number_of_elements_out );
    CHECK_FOR_CUDA_ERROR();
    return out;
  }

  // Crop
  template<class T, unsigned int D> __global__ void crop_kernel
  ( vector_td<unsigned int,D> offset, vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
    const T * __restrict__ in, T * __restrict__ out, unsigned int num_batches, unsigned int num_elements )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const unsigned int frame_offset = idx/num_elements;
    
    if( idx < num_elements*num_batches ){
      const typename uintd<D>::Type co = idx_to_co( idx-frame_offset*num_elements, matrix_size_out );
      const typename uintd<D>::Type co_os = offset + co;
      const unsigned int in_idx = co_to_idx(co_os, matrix_size_in)+frame_offset*prod(matrix_size_in);
      out[idx] = in[in_idx];
    }
  }

  // Crop
  template<class T, unsigned int D>
  void crop( const vector_td<size_t,D>& offset, const vector_td<size_t,D>& crop_size, const cuNDArray<T>& in, cuNDArray<T>&out )
  {
    auto out_dims = to_std_vector(crop_size);
    for (int i = D; i < in.get_number_of_dimensions(); i++)
      out_dims.push_back(in.get_size(i));

    if (!out.dimensions_equal(out_dims))
      out.create(out_dims);


    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *in.get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *out.get_dimensions() );

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in.get_number_of_dimensions(); d++ ){
        number_of_batches *= in.get_size(d);
      }

           if( weak_greater(offset+matrix_size_out, matrix_size_in) ){
             throw std::runtime_error( "crop: cropping size mismatch");
           }

           // Setup block/grid dimensions
           dim3 blockDim; dim3 gridDim;
         setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

         // Invoke kernel
         crop_kernel<T,D><<< gridDim, blockDim >>>
           ( vector_td<unsigned int,D>(offset), vector_td<unsigned int,D>(matrix_size_in), vector_td<unsigned int,D>(matrix_size_out),
           in.get_data_ptr(), out.get_data_ptr(), number_of_batches, prod(matrix_size_out) );
    
    CHECK_FOR_CUDA_ERROR();
  }

  template<class T, unsigned int D> cuNDArray<T>
  crop<T,D>( const vector_td<size_t,D>& offset, const vector_td<size_t,D>& size, const cuNDArray<T>& in )
  {
    cuNDArray<T>  result;
    crop<T,D>(offset,size, in, result);
    return result;
  }  

  // Expand and zero fill
  template<class T, unsigned int D> 
  __global__ void pad_kernel( vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
                              const T * __restrict__ in, T * __restrict__ out, unsigned int number_of_batches, unsigned int num_elements, T val )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const unsigned int frame_offset = idx/num_elements;

    if( idx < num_elements*number_of_batches ){

      const typename uintd<D>::Type co_out = idx_to_co( idx-frame_offset*num_elements, matrix_size_out );
      typename uintd<D>::Type offset;
      for(unsigned int d=0; d<D; d++)
      {
          offset[d] = matrix_size_out[d]/2 - matrix_size_in[d]/2;
      }

      T _out;
      bool inside = (co_out>=offset) && (co_out<(matrix_size_in+offset));

      if( inside )
        _out = in[co_to_idx(co_out-offset, matrix_size_in)+frame_offset*prod(matrix_size_in)];
      else{      
        _out = val;
      }

      out[idx] = _out;
    }
  }

  template<class T, unsigned int D> 
  void pad( const cuNDArray<T>&in, cuNDArray<T>&out, T val )
  { 
    if( in.get_number_of_dimensions() != out.get_number_of_dimensions() ){
      throw std::runtime_error("pad: image dimensions mismatch");
    }

    if( in.get_number_of_dimensions() < D ){
      std::stringstream ss;
      ss << "pad: number of image dimensions should be at least " << D;
      throw std::runtime_error(ss.str());
    }

    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *in.get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *out.get_dimensions() );

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in.get_number_of_dimensions(); d++ ){
      number_of_batches *= in.get_size(d);
    }

    if( weak_greater(matrix_size_in,matrix_size_out) ){
      throw std::runtime_error("pad: size mismatch, cannot expand");
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

    // Invoke kernel
    pad_kernel<T,D><<< gridDim, blockDim >>> 
      ( vector_td<unsigned int,D>(matrix_size_in), vector_td<unsigned int,D>(matrix_size_out),
        in.get_data_ptr(), out.get_data_ptr(), number_of_batches, prod(matrix_size_out), val );

    CHECK_FOR_CUDA_ERROR();
  }

  template<class T, unsigned int D> cuNDArray<T>
  pad<T,D>( const vector_td<size_t,D>& size, const cuNDArray<T>&in, T val )
  {
    std::vector<size_t> dims = to_std_vector(size);
    for( unsigned int d=D; d<in.get_number_of_dimensions(); d++ ){
      dims.push_back(in.get_size(d));
    }
    cuNDArray<T>  result( &dims );
    pad<T,D>(in, result, val);
    return result;
  }

  template<class T, unsigned int D> 
  __global__ void fill_border_kernel( vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
                                      T *image, unsigned int number_of_batches, unsigned int number_of_elements, T val )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

    if( idx < number_of_elements ){
      const vector_td<unsigned int,D> co_out = idx_to_co( idx, matrix_size_out );
      const vector_td<unsigned int,D> offset = (matrix_size_out-matrix_size_in)>>1;
      if( weak_less( co_out, offset ) || weak_greater_equal( co_out, matrix_size_in+offset ) ){
	      for( unsigned int batch=0; batch<number_of_batches; batch++ ){
          image[idx+batch*number_of_elements] = val;
        }
      }
      else
	      ; // do nothing
    }
  }

  // Zero fill border (rectangular)
  template<class T, unsigned int D> 
  void fill_border( const vector_td<size_t,D>& matrix_size_in, cuNDArray<T>&in_out, T val )
  { 
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *in_out.get_dimensions() );

    if( weak_greater(matrix_size_in, matrix_size_out) ){
      throw std::runtime_error("fill_border: size mismatch, cannot zero fill");
    }

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in_out.get_number_of_dimensions(); d++ ){
      number_of_batches *= in_out.get_size(d);
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim );

    // Invoke kernel
    fill_border_kernel<T,D><<< gridDim, blockDim >>>
      ( vector_td<unsigned int,D>(matrix_size_in), vector_td<unsigned int,D>(matrix_size_out),
        in_out.get_data_ptr(), number_of_batches, prod(matrix_size_out), val );

    CHECK_FOR_CUDA_ERROR();
  }


  template<class T, unsigned int D>
  __global__ void fill_border_kernel( typename realType<T>::Type radius, vector_td<int,D> matrix_size,
                                      T *image, unsigned int number_of_batches, unsigned int number_of_elements, T val )
  {
    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

    if( idx < number_of_elements ){
      const vector_td<typename realType<T>::Type,D> co_out( (matrix_size>>1) - idx_to_co( idx, matrix_size ));
      if(  norm(co_out) > radius ){
	      for( unsigned int batch=0; batch<number_of_batches; batch++ ){
          image[idx+batch*number_of_elements] = val;
        }
      }
      else
	      ; // do nothing
    }
  }

  // Zero fill border (radial)
  template<class T, unsigned int D>
  void fill_border( typename realType<T>::Type radius, cuNDArray<T>&in_out, T val )
  {
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *in_out.get_dimensions() );


    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in_out.get_number_of_dimensions(); d++ ){
      number_of_batches *= in_out.get_size(d);
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim );

    // Invoke kernel
    fill_border_kernel<T,D><<< gridDim, blockDim >>>
      (radius, vector_td<int,D>(matrix_size_out),
        in_out.get_data_ptr(), number_of_batches, prod(matrix_size_out), val );

    CHECK_FOR_CUDA_ERROR();
  }
  template<class T, unsigned int D> __global__ void 
  upsample_kernel( typename uintd<D>::Type matrix_size_in,
                   typename uintd<D>::Type matrix_size_out,
                   unsigned int num_batches,
                   const T * __restrict__ image_in,
                   T * __restrict__ image_out )
  {
    typedef typename realType<T>::Type REAL;
    
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const unsigned int num_elements_out = prod(matrix_size_out);
    
    if( idx < num_elements_out*num_batches ){
      
      const unsigned int batch = idx/num_elements_out;
      const unsigned int batch_offset_in = batch*prod(matrix_size_in);
      
      const typename uintd<D>::Type co_out = idx_to_co( idx-batch*num_elements_out, matrix_size_out );
      const typename uintd<D>::Type co_in = co_out >> 1;
      const typename uintd<D>::Type ones(1);
      const typename uintd<D>::Type twos(2);
      const typename uintd<D>::Type offset = co_out%twos;
      
      const unsigned int num_cells = 1 << D;
      
      T cellsum(0);
      unsigned int count = 0;
      
      for( unsigned int i=0; i<num_cells; i++ ){
        
        const typename uintd<D>::Type stride = idx_to_co( i, twos );
        
        if( offset >= stride ){
          cellsum += image_in[batch_offset_in+co_to_idx(amin(co_in+stride, matrix_size_in-ones), matrix_size_in)];
          count++;
        }
      }

      image_out[idx] = cellsum / REAL(count);
    }
  }

  //
  // Linear upsampling by a factor of two (on a D-dimensional grid) 
  // Note that this operator is the transpose of the downsampling operator below by design
  // - based on Briggs et al, A Multigrid Tutorial 2nd edition, pp. 34-35
  // 
  
  template<class T, unsigned int D>  cuNDArray<T> upsample( const cuNDArray<T>& in )
	{

    std::vector<size_t> dims_out = *in.get_dimensions();
    for( unsigned int i=0; i<D; i++ ) dims_out[i] <<= 1;
    cuNDArray<T> out(&dims_out);
    upsample<T,D>( in, out );
    return out;
	}

  template<class T, unsigned int D> void upsample( const cuNDArray<T>&in, cuNDArray<T>&out )
  {
    typename uint64d<D>::Type matrix_size_in  = from_std_vector<size_t,D>( *in.get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *out.get_dimensions() );

    if( (matrix_size_in<<1) != matrix_size_out ){
      throw std::runtime_error("upsample: arrays do not correspond to upsampling by a factor of two");
    }

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<out.get_number_of_dimensions(); d++ ){
      number_of_batches *= out.get_size(d);
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

    // Invoke kernel
    upsample_kernel<T,D><<< gridDim, blockDim >>>
      ( vector_td<unsigned int,D>(matrix_size_in), vector_td<unsigned int,D>(matrix_size_out),
        number_of_batches, in.get_data_ptr(), out.get_data_ptr() );

    CHECK_FOR_CUDA_ERROR();    
  }
  //
  // Linear downsampling by a factor of two (on a D-dimensional grid)
  // Note that this operator is the transpose of the upsampling operator above by design
  // - based on Briggs et al, A Multigrid Tutorial 2nd edition, pp. 36.
  // 

  template<class T, unsigned int D> __global__ void 
  downsample_kernel( typename intd<D>::Type matrix_size_in,
                     typename intd<D>::Type matrix_size_out,
                     int num_batches,
                     const T * __restrict__ image_in,
                     T * __restrict__ image_out )
  {
    typedef typename realType<T>::Type REAL;
    
    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const int num_elements_out = prod(matrix_size_out);
    
    if( idx < num_elements_out*num_batches ){
      
      const int batch = idx/num_elements_out;
      const int batch_offset_in = batch*prod(matrix_size_in);
      
      const typename intd<D>::Type co_out = idx_to_co( idx-batch*num_elements_out, matrix_size_out );
      const typename intd<D>::Type co_in = co_out << 1;
      
      T cellsum[D+1];
      for( unsigned int d=0; d<D+1; d++ ){
        cellsum[d] = T(0);
      }
      
      //const int num_cells = pow(3,D); // no pow for integers on device
      int num_cells = 1; 
      for( int i=0; i<D; i++ ) num_cells *=3;

      const REAL denominator = pow(REAL(4),REAL(D));
      
      for( int i=0; i<num_cells; i++ ){
        
        const typename intd<D>::Type zeros(0);
        const typename intd<D>::Type ones(1);
        const typename intd<D>::Type threes(3);
        const typename intd<D>::Type stride = idx_to_co(i,threes)-ones; // in the range [-1;1]^D
        
        int distance = 0;
        for( int d=0; d<D; d++ ){
          if( abs(stride[d])>0 )
            distance++;
        }
        
        cellsum[distance] += image_in[batch_offset_in+co_to_idx(amax(zeros, amin(matrix_size_in-ones,co_in+stride)), matrix_size_in)];
      }
      
      T res = T(0);
      
      for( unsigned int d=0; d<D+1; d++ ){
        res += (REAL(1<<(D-d))*cellsum[d]);
      }
      
      image_out[idx] = res / denominator;
    }
  }

  template<class T, unsigned int D>  cuNDArray<T> downsample( const cuNDArray<T>& in )
  {
    
    std::vector<size_t> dims_out = *in.get_dimensions();
    for( unsigned int i=0; i<D; i++ ) dims_out[i] >>= 1;
    cuNDArray<T> out(&dims_out);
    downsample<T,D>( in, out );
    return out;
  }

  template<class T, unsigned int D> void downsample(const  cuNDArray<T>& in, cuNDArray<T>& out )
  {

    typename uint64d<D>::Type matrix_size_in  = from_std_vector<size_t,D>( *in.get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *out.get_dimensions() );

    if( (matrix_size_in>>1) != matrix_size_out ){
      throw std::runtime_error("downsample: arrays do not correspond to downsampling by a factor of two");
    }

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<out.get_number_of_dimensions(); d++ ){
      number_of_batches *= out.get_size(d);
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

    // Invoke kernel
    downsample_kernel<T,D><<< gridDim, blockDim >>>
      ( vector_td<int,D>(matrix_size_in), vector_td<int,D>(matrix_size_out),
        (int)number_of_batches, in.get_data_ptr(), out.get_data_ptr() );

    CHECK_FOR_CUDA_ERROR();    
  }

  //
  // Instantiation
  //

  template EXPORTGPUCORE  cuNDArray<float > permute( const cuNDArray<float>&, const std::vector<size_t>&);
  template EXPORTGPUCORE  cuNDArray<double > permute( const cuNDArray<double>&, const std::vector<size_t>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > permute(const cuNDArray<float_complext>&, const std::vector<size_t>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > permute(const cuNDArray<double_complext>&, const std::vector<size_t>& );

  template EXPORTGPUCORE void permute( const cuNDArray<float>&, cuNDArray<float>&, const std::vector<size_t>&);
  template EXPORTGPUCORE void permute( const cuNDArray<double>&, cuNDArray<double>&, const std::vector<size_t>&);
  template EXPORTGPUCORE void permute( const cuNDArray<float_complext>&, cuNDArray<float_complext>&, const std::vector<size_t>&);
  template EXPORTGPUCORE void permute( const cuNDArray<double_complext>&, cuNDArray<double_complext>&, const std::vector<size_t>&);

  template EXPORTGPUCORE  cuNDArray<float > shift_dim( const cuNDArray<float>&, int );
  template EXPORTGPUCORE  cuNDArray<double > shift_dim( const cuNDArray<double>&, int );
  template EXPORTGPUCORE  cuNDArray<float_complext > shift_dim( const cuNDArray<float_complext>&, int );
  template EXPORTGPUCORE  cuNDArray<double_complext > shift_dim( const cuNDArray<double_complext>&, int );

  template EXPORTGPUCORE void shift_dim( const cuNDArray<float>&, cuNDArray<float>&, int shift );
  template EXPORTGPUCORE void shift_dim( const cuNDArray<double>&, cuNDArray<double>&, int shift );
  template EXPORTGPUCORE void shift_dim( const cuNDArray<float_complext>&, cuNDArray<float_complext>&, int shift );
  template EXPORTGPUCORE void shift_dim( const cuNDArray<double_complext>&, cuNDArray<double_complext>&, int shift );

  template EXPORTGPUCORE  cuNDArray<float > expand<float>( const cuNDArray<float>&, size_t);
  template EXPORTGPUCORE  cuNDArray<double > expand<double>( const cuNDArray<double>&, size_t);
  template EXPORTGPUCORE  cuNDArray<float_complext > expand<float_complext>( const cuNDArray<float_complext>&, size_t);
  template EXPORTGPUCORE  cuNDArray<double_complext > expand<double_complext>( const cuNDArray<double_complext>&, size_t);

  template EXPORTGPUCORE  cuNDArray<float > crop<float,1>(const vector_td<size_t,1>&, const vector_td<size_t,1>&, const cuNDArray<float>&);
  template EXPORTGPUCORE  cuNDArray<float > crop<float,2>(const vector_td<size_t,2>&, const vector_td<size_t,2>&, const cuNDArray<float>&);
  template EXPORTGPUCORE  cuNDArray<float > crop<float,3>(const vector_td<size_t,3>&, const vector_td<size_t,3>&, const cuNDArray<float>&);
  template EXPORTGPUCORE  cuNDArray<float > crop<float,4>(const vector_td<size_t,4>&, const vector_td<size_t,4>&, const cuNDArray<float>&);

  template EXPORTGPUCORE  cuNDArray<float_complext > crop<float_complext,1>( const typename uint64d<1>::Type&, const typename uint64d<1>::Type&, const cuNDArray<float_complext>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > crop<float_complext,2>( const typename uint64d<2>::Type&, const typename uint64d<2>::Type&, const cuNDArray<float_complext>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > crop<float_complext,3>( const typename uint64d<3>::Type&, const typename uint64d<3>::Type&, const cuNDArray<float_complext>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > crop<float_complext,4>( const typename uint64d<4>::Type&, const typename uint64d<4>::Type&, const cuNDArray<float_complext>&);

  template EXPORTGPUCORE void crop<float,1>( const vector_td<size_t,1>&, const vector_td<size_t,1>&, const cuNDArray<float>&, cuNDArray<float>&);
  template EXPORTGPUCORE void crop<float,2>( const vector_td<size_t,2>&, const vector_td<size_t,2>&, const cuNDArray<float>&, cuNDArray<float>&);
  template EXPORTGPUCORE void crop<float,3>( const vector_td<size_t,3>&, const vector_td<size_t,3>&, const cuNDArray<float>&, cuNDArray<float>&);
  template EXPORTGPUCORE void crop<float,4>( const vector_td<size_t,4>&, const vector_td<size_t,4>&, const cuNDArray<float>&, cuNDArray<float>&);

  template EXPORTGPUCORE void crop<complext<float>,1>( const uint64d1&, const uint64d1&, const cuNDArray<complext<float> >&, cuNDArray< complext<float> >&);
  template EXPORTGPUCORE void crop<complext<float>,2>( const uint64d2&, const uint64d2&, const cuNDArray<complext<float> >&, cuNDArray< complext<float> >&);
  template EXPORTGPUCORE void crop<complext<float>,3>( const uint64d3&, const uint64d3&, const cuNDArray<complext<float> >&, cuNDArray< complext<float> >&);
  template EXPORTGPUCORE void crop<complext<float>,4>( const uint64d4&, const uint64d4&, const cuNDArray<complext<float> >&, cuNDArray< complext<float> >&);

  template EXPORTGPUCORE  cuNDArray<float > pad<float,1>( const uint64d1&, const cuNDArray<float>&, float );
  template EXPORTGPUCORE  cuNDArray<float > pad<float,2>( const uint64d2&, const cuNDArray<float>&, float );
  template EXPORTGPUCORE  cuNDArray<float > pad<float,3>( const uint64d3&, const cuNDArray<float>&, float );
  template EXPORTGPUCORE  cuNDArray<float > pad<float,4>( const uint64d4&, const cuNDArray<float>&, float );

  template EXPORTGPUCORE  cuNDArray<float_complext > pad<float_complext,1>( const uint64d1&, const cuNDArray<float_complext>&, float_complext );
  template EXPORTGPUCORE  cuNDArray<float_complext > pad<float_complext,2>( const uint64d2&,  const cuNDArray<float_complext>&, float_complext );
  template EXPORTGPUCORE  cuNDArray<float_complext > pad<float_complext,3>( const uint64d3&, const cuNDArray<float_complext>&, float_complext );
  template EXPORTGPUCORE  cuNDArray<float_complext > pad<float_complext,4>( const uint64d4&, const cuNDArray<float_complext>&, float_complext );

  template EXPORTGPUCORE void pad<float,1>( const cuNDArray<float>&, cuNDArray<float>&, float);
  template EXPORTGPUCORE void pad<float,2>( const cuNDArray<float>&, cuNDArray<float>&, float);
  template EXPORTGPUCORE void pad<float,3>( const cuNDArray<float>&, cuNDArray<float>&, float);
  template EXPORTGPUCORE void pad<float,4>( const cuNDArray<float>&, cuNDArray<float>&, float);

  template EXPORTGPUCORE void pad<float_complext,1>( const cuNDArray<float_complext>&, cuNDArray<float_complext>&, float_complext);
  template EXPORTGPUCORE void pad<float_complext,2>( const cuNDArray<float_complext>&, cuNDArray<float_complext>&, float_complext);
  template EXPORTGPUCORE void pad<float_complext,3>( const cuNDArray<float_complext>&, cuNDArray<float_complext>&, float_complext);
  template EXPORTGPUCORE void pad<float_complext,4>( const cuNDArray<float_complext>&, cuNDArray<float_complext>&, float_complext);

  template EXPORTGPUCORE void fill_border<float,1>(const uint64d1&, cuNDArray<float>&,float);
  template EXPORTGPUCORE void fill_border<float,2>(const uint64d2&, cuNDArray<float>&,float);
  template EXPORTGPUCORE void fill_border<float,3>(const uint64d3&, cuNDArray<float>&,float);
  template EXPORTGPUCORE void fill_border<float,4>(const uint64d4&, cuNDArray<float>&,float);
  template EXPORTGPUCORE void fill_border<float,1>(float, cuNDArray<float>&,float);
	template EXPORTGPUCORE void fill_border<float,2>(float, cuNDArray<float>&,float);
	template EXPORTGPUCORE void fill_border<float,3>(float, cuNDArray<float>&,float);
	template EXPORTGPUCORE void fill_border<float,4>(float, cuNDArray<float>&,float);

  template EXPORTGPUCORE void fill_border<float_complext,1>(const uint64d1&, cuNDArray<float_complext>&,float_complext);
  template EXPORTGPUCORE void fill_border<float_complext,2>(const uint64d2&, cuNDArray<float_complext>&,float_complext);
  template EXPORTGPUCORE void fill_border<float_complext,3>(const uint64d3&, cuNDArray<float_complext>&,float_complext);
  template EXPORTGPUCORE void fill_border<float_complext,4>(const uint64d4&, cuNDArray<float_complext>&,float_complext);
  template EXPORTGPUCORE void fill_border<float_complext,1>(float, cuNDArray<float_complext>&,float_complext);
	template EXPORTGPUCORE void fill_border<float_complext,2>(float, cuNDArray<float_complext>&,float_complext);
	template EXPORTGPUCORE void fill_border<float_complext,3>(float, cuNDArray<float_complext>&,float_complext);
	template EXPORTGPUCORE void fill_border<float_complext,4>(float, cuNDArray<float_complext>&,float_complext);



  template EXPORTGPUCORE  cuNDArray<double > crop<double,1>(const vector_td<size_t,1>&, const vector_td<size_t,1>&, const cuNDArray<double>&);
  template EXPORTGPUCORE  cuNDArray<double > crop<double,2>(const vector_td<size_t,2>&, const vector_td<size_t,2>&, const cuNDArray<double>&);
  template EXPORTGPUCORE  cuNDArray<double > crop<double,3>(const vector_td<size_t,3>&, const vector_td<size_t,3>&, const cuNDArray<double>&);
  template EXPORTGPUCORE  cuNDArray<double > crop<double,4>(const vector_td<size_t,4>&, const vector_td<size_t,4>&, const cuNDArray<double>&);

  template EXPORTGPUCORE  cuNDArray<double_complext > crop<double_complext,1>(const vector_td<size_t,1>&, const vector_td<size_t,1>&, const cuNDArray<double_complext>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > crop<double_complext,2>(const vector_td<size_t,2>&, const vector_td<size_t,2>&, const cuNDArray<double_complext>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > crop<double_complext,3>(const vector_td<size_t,3>&, const vector_td<size_t,3>&, const cuNDArray<double_complext>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > crop<double_complext,4>(const vector_td<size_t,4>&, const vector_td<size_t,4>&, const cuNDArray<double_complext>&);

  template EXPORTGPUCORE void crop<double,1>(const uint64d1&, const uint64d1&, const cuNDArray<double>&, cuNDArray<double>&);
  template EXPORTGPUCORE void crop<double,2>(const uint64d2&, const uint64d2&, const cuNDArray<double>&, cuNDArray<double>&);
  template EXPORTGPUCORE void crop<double,3>(const uint64d3&, const uint64d3&, const cuNDArray<double>&, cuNDArray<double>&);
  template EXPORTGPUCORE void crop<double,4>(const uint64d4&, const uint64d4&, const cuNDArray<double>&, cuNDArray<double>&);

  template EXPORTGPUCORE void crop<complext<double>,1>(const uint64d1&, const uint64d1&, const cuNDArray<complext<double> >&, cuNDArray< complext<double> >&);
  template EXPORTGPUCORE void crop<complext<double>,2>(const uint64d2&, const uint64d2&, const cuNDArray<complext<double> >&, cuNDArray< complext<double> >&);
  template EXPORTGPUCORE void crop<complext<double>,3>(const uint64d3&, const uint64d3&, const cuNDArray<complext<double> >&, cuNDArray< complext<double> >&);
  template EXPORTGPUCORE void crop<complext<double>,4>(const uint64d4&, const uint64d4&, const cuNDArray<complext<double> >&, cuNDArray< complext<double> >&);

  template EXPORTGPUCORE  cuNDArray<double > pad<double,1>(const uint64d1&, const cuNDArray<double>&, double );
  template EXPORTGPUCORE  cuNDArray<double > pad<double,2>(const uint64d2&, const cuNDArray<double>&, double );
  template EXPORTGPUCORE  cuNDArray<double > pad<double,3>(const uint64d3&, const cuNDArray<double>&, double );
  template EXPORTGPUCORE  cuNDArray<double > pad<double,4>(const uint64d4&, const cuNDArray<double>&, double );

  template EXPORTGPUCORE  cuNDArray<double_complext > pad<double_complext,1>(const uint64d1&, const cuNDArray<double_complext>&, double_complext );
  template EXPORTGPUCORE  cuNDArray<double_complext > pad<double_complext,2>( const uint64d2&,const cuNDArray<double_complext>&, double_complext );
  template EXPORTGPUCORE  cuNDArray<double_complext > pad<double_complext,3>( const uint64d3&,const cuNDArray<double_complext>&, double_complext );
  template EXPORTGPUCORE  cuNDArray<double_complext > pad<double_complext,4>( const uint64d4&,const cuNDArray<double_complext>&, double_complext );

  template EXPORTGPUCORE void pad<double,1>( const cuNDArray<double>&, cuNDArray<double>&, double);
  template EXPORTGPUCORE void pad<double,2>( const cuNDArray<double>&, cuNDArray<double>&, double);
  template EXPORTGPUCORE void pad<double,3>( const cuNDArray<double>&, cuNDArray<double>&, double);
  template EXPORTGPUCORE void pad<double,4>( const cuNDArray<double>&, cuNDArray<double>&, double);

  template EXPORTGPUCORE void pad<double_complext,1>( const cuNDArray<double_complext>&, cuNDArray<double_complext>&, double_complext);
  template EXPORTGPUCORE void pad<double_complext,2>( const cuNDArray<double_complext>&, cuNDArray<double_complext>&, double_complext);
  template EXPORTGPUCORE void pad<double_complext,3>( const cuNDArray<double_complext>&, cuNDArray<double_complext>&, double_complext);
  template EXPORTGPUCORE void pad<double_complext,4>( const cuNDArray<double_complext>&, cuNDArray<double_complext>&, double_complext);

  template EXPORTGPUCORE void fill_border<double,1>(const uint64d1&, cuNDArray<double>&,double);
  template EXPORTGPUCORE void fill_border<double,2>(const uint64d2&, cuNDArray<double>&,double);
  template EXPORTGPUCORE void fill_border<double,3>(const uint64d3&, cuNDArray<double>&,double);
  template EXPORTGPUCORE void fill_border<double,4>(const uint64d4&, cuNDArray<double>&,double);
  template EXPORTGPUCORE void fill_border<double,1>(double, cuNDArray<double>&,double);
	template EXPORTGPUCORE void fill_border<double,2>(double, cuNDArray<double>&,double);
	template EXPORTGPUCORE void fill_border<double,3>(double, cuNDArray<double>&,double);
	template EXPORTGPUCORE void fill_border<double,4>(double, cuNDArray<double>&,double);

  template EXPORTGPUCORE void fill_border<double_complext,1>(const uint64d1&, cuNDArray<double_complext>&,double_complext);
  template EXPORTGPUCORE void fill_border<double_complext,2>(const uint64d2&, cuNDArray<double_complext>&,double_complext);
  template EXPORTGPUCORE void fill_border<double_complext,3>(const uint64d3&, cuNDArray<double_complext>&,double_complext);
  template EXPORTGPUCORE void fill_border<double_complext,4>(const uint64d4&, cuNDArray<double_complext>&,double_complext);
  template EXPORTGPUCORE void fill_border<double_complext,1>(double, cuNDArray<double_complext>&,double_complext);
	template EXPORTGPUCORE void fill_border<double_complext,2>(double, cuNDArray<double_complext>&,double_complext);
	template EXPORTGPUCORE void fill_border<double_complext,3>(double, cuNDArray<double_complext>&,double_complext);
	template EXPORTGPUCORE void fill_border<double_complext,4>(double, cuNDArray<double_complext>&,double_complext);


  template EXPORTGPUCORE  cuNDArray<float > upsample<float,1>(const cuNDArray<float>&);
  template EXPORTGPUCORE  cuNDArray<float > upsample<float,2>(const cuNDArray<float>&);
  template EXPORTGPUCORE  cuNDArray<float > upsample<float,3>(const cuNDArray<float>&);
  template EXPORTGPUCORE  cuNDArray<float > upsample<float,4>(const cuNDArray<float>&);

  template EXPORTGPUCORE  cuNDArray<float_complext > upsample<float_complext,1>(const cuNDArray<float_complext>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > upsample<float_complext,2>(const cuNDArray<float_complext>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > upsample<float_complext,3>(const cuNDArray<float_complext>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > upsample<float_complext,4>(const cuNDArray<float_complext>&);

  template EXPORTGPUCORE  cuNDArray<double > upsample<double,1>(const cuNDArray<double>&);
  template EXPORTGPUCORE  cuNDArray<double > upsample<double,2>(const cuNDArray<double>&);
  template EXPORTGPUCORE  cuNDArray<double > upsample<double,3>(const cuNDArray<double>&);
  template EXPORTGPUCORE  cuNDArray<double > upsample<double,4>(const cuNDArray<double>&);

  template EXPORTGPUCORE  cuNDArray<double_complext > upsample<double_complext,1>(const cuNDArray<double_complext>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > upsample<double_complext,2>(const cuNDArray<double_complext>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > upsample<double_complext,3>(const cuNDArray<double_complext>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > upsample<double_complext,4>(const cuNDArray<double_complext>&);

  template EXPORTGPUCORE void upsample<float,1>(const cuNDArray<float>&, cuNDArray<float>&);
  template EXPORTGPUCORE void upsample<float,2>(const cuNDArray<float>&, cuNDArray<float>&);
  template EXPORTGPUCORE void upsample<float,3>(const cuNDArray<float>&, cuNDArray<float>&);
  template EXPORTGPUCORE void upsample<float,4>(const cuNDArray<float>&, cuNDArray<float>&);

  template EXPORTGPUCORE void upsample<float_complext,1>(const cuNDArray<float_complext>&, cuNDArray<float_complext>&);
  template EXPORTGPUCORE void upsample<float_complext,2>(const cuNDArray<float_complext>&, cuNDArray<float_complext>&);
  template EXPORTGPUCORE void upsample<float_complext,3>(const cuNDArray<float_complext>&, cuNDArray<float_complext>&);
  template EXPORTGPUCORE void upsample<float_complext,4>(const cuNDArray<float_complext>&, cuNDArray<float_complext>&);

  template EXPORTGPUCORE void upsample<double,1>(const cuNDArray<double>&, cuNDArray<double>&);
  template EXPORTGPUCORE void upsample<double,2>(const cuNDArray<double>&, cuNDArray<double>&);
  template EXPORTGPUCORE void upsample<double,3>(const cuNDArray<double>&, cuNDArray<double>&);
  template EXPORTGPUCORE void upsample<double,4>(const cuNDArray<double>&, cuNDArray<double>&);

  template EXPORTGPUCORE void upsample<double_complext,1>(const cuNDArray<double_complext>&, cuNDArray<double_complext>&);
  template EXPORTGPUCORE void upsample<double_complext,2>(const cuNDArray<double_complext>&, cuNDArray<double_complext>&);
  template EXPORTGPUCORE void upsample<double_complext,3>(const cuNDArray<double_complext>&, cuNDArray<double_complext>&);
  template EXPORTGPUCORE void upsample<double_complext,4>(const cuNDArray<double_complext>&, cuNDArray<double_complext>&);

  template EXPORTGPUCORE  cuNDArray<float > downsample<float,1>(const cuNDArray<float>&);
  template EXPORTGPUCORE  cuNDArray<float > downsample<float,2>(const cuNDArray<float>&);
  template EXPORTGPUCORE  cuNDArray<float > downsample<float,3>(const cuNDArray<float>&);
  template EXPORTGPUCORE  cuNDArray<float > downsample<float,4>(const cuNDArray<float>&);

  template EXPORTGPUCORE  cuNDArray<float_complext > downsample<float_complext,1>(const cuNDArray<float_complext>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > downsample<float_complext,2>(const cuNDArray<float_complext>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > downsample<float_complext,3>(const cuNDArray<float_complext>&);
  template EXPORTGPUCORE  cuNDArray<float_complext > downsample<float_complext,4>(const cuNDArray<float_complext>&);

  template EXPORTGPUCORE  cuNDArray<double > downsample<double,1>(const cuNDArray<double>&);
  template EXPORTGPUCORE  cuNDArray<double > downsample<double,2>(const cuNDArray<double>&);
  template EXPORTGPUCORE  cuNDArray<double > downsample<double,3>(const cuNDArray<double>&);
  template EXPORTGPUCORE  cuNDArray<double > downsample<double,4>(const cuNDArray<double>&);

  template EXPORTGPUCORE  cuNDArray<double_complext > downsample<double_complext,1>(const cuNDArray<double_complext>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > downsample<double_complext,2>(const cuNDArray<double_complext>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > downsample<double_complext,3>(const cuNDArray<double_complext>&);
  template EXPORTGPUCORE  cuNDArray<double_complext > downsample<double_complext,4>(const cuNDArray<double_complext>&);

  template EXPORTGPUCORE void downsample<float,1>(const cuNDArray<float>&, cuNDArray<float>&);
  template EXPORTGPUCORE void downsample<float,2>(const cuNDArray<float>&, cuNDArray<float>&);
  template EXPORTGPUCORE void downsample<float,3>(const cuNDArray<float>&, cuNDArray<float>&);
  template EXPORTGPUCORE void downsample<float,4>(const cuNDArray<float>&, cuNDArray<float>&);

  template EXPORTGPUCORE void downsample<float_complext,1>(const cuNDArray<float_complext>&, cuNDArray<float_complext>&);
  template EXPORTGPUCORE void downsample<float_complext,2>(const cuNDArray<float_complext>&, cuNDArray<float_complext>&);
  template EXPORTGPUCORE void downsample<float_complext,3>(const cuNDArray<float_complext>&, cuNDArray<float_complext>&);
  template EXPORTGPUCORE void downsample<float_complext,4>(const cuNDArray<float_complext>&, cuNDArray<float_complext>&);

  template EXPORTGPUCORE void downsample<double,1>(const cuNDArray<double>&, cuNDArray<double>&);
  template EXPORTGPUCORE void downsample<double,2>(const cuNDArray<double>&, cuNDArray<double>&);
  template EXPORTGPUCORE void downsample<double,3>(const cuNDArray<double>&, cuNDArray<double>&);
  template EXPORTGPUCORE void downsample<double,4>(const cuNDArray<double>&, cuNDArray<double>&);

  template EXPORTGPUCORE void downsample<double_complext,1>(const cuNDArray<double_complext>&, cuNDArray<double_complext>&);
  template EXPORTGPUCORE void downsample<double_complext,2>(const cuNDArray<double_complext>&, cuNDArray<double_complext>&);
  template EXPORTGPUCORE void downsample<double_complext,3>(const cuNDArray<double_complext>&, cuNDArray<double_complext>&);
  template EXPORTGPUCORE void downsample<double_complext,4>(const cuNDArray<double_complext>&, cuNDArray<double_complext>&);


  // We can probably instantiate the functions below functionsfor many more types? E.g. arrays of floatd2. 
  // For now we just introduce what we have needed...
  //

  template EXPORTGPUCORE  cuNDArray<floatd2 > expand<floatd2>( const cuNDArray<floatd2>&, size_t);
}
