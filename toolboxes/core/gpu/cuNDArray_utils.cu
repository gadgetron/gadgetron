#include "cuNDArray_utils.h"
#include "vector_td_utilities.h"
#include "cudaDeviceManager.h"
#include "setup_grid.h"
#include "cuNDArray_math.h"

#include <math_functions.h>
#include <cmath>

namespace Gadgetron {

  template <class T> 
  __global__ void cuNDArray_permute_kernel(const  T*  __restrict__ in, T* __restrict__ out,
                                            unsigned int ndim,
                                            const unsigned int* __restrict__ dims,
                                            const unsigned int* __restrict__ strides_out,
                                            unsigned int elements,
                                            int shift_mode)
  {
    unsigned int idx_in = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    unsigned int idx_out = 0;
    unsigned int idx_in_tmp = idx_in;

    if (idx_in < elements) {

      unsigned int cur_index;
      for (unsigned int i = 0; i < ndim; i++) {
        unsigned int idx_in_remainder = idx_in_tmp / dims[i];
        cur_index = idx_in_tmp-(idx_in_remainder*dims[i]); //cur_index = idx_in_tmp%dims[i];
        if (shift_mode < 0) { //IFFTSHIFT
          idx_out += ((cur_index+(dims[i]>>1))%dims[i])*strides_out[i];
        } else if (shift_mode > 0) { //FFTSHIFT
          idx_out += ((cur_index+((dims[i]+1)>>1))%dims[i])*strides_out[i];
        } else {
          idx_out += cur_index*strides_out[i];
        }
        idx_in_tmp = idx_in_remainder;
      }
      out[idx_in] = in[idx_out];
    }
  }

  template <class T> void cuNDArray_permute( cuNDArray<T>* in,
                                             cuNDArray<T>* out,
                                             std::vector<size_t> *order,
                                             int shift_mode)
  {    
    if( out == 0x0 ){
      throw cuda_error("cuNDArray_permute(internal): 0x0 output");
    }


    cudaError_t err;

    T* in_ptr = in->get_data_ptr();
    T* out_ptr = 0;

    if (out) {
      out_ptr = out->get_data_ptr();
    } else {
      if (cudaMalloc((void**) &out_ptr, in->get_number_of_elements()*sizeof(T)) != cudaSuccess) {
        throw cuda_error("cuNDArray_permute : Error allocating CUDA memory");
      }
    }

    unsigned int* dims        = new unsigned int[in->get_number_of_dimensions()];
    unsigned int* strides_out = new unsigned int[in->get_number_of_dimensions()];

    if (!dims || !strides_out) {
      throw cuda_error("cuNDArray_permute: failed to allocate temporary storage for arrays");
    }

    for (unsigned int i = 0; i < in->get_number_of_dimensions(); i++) {
      dims[i] = (*in->get_dimensions())[(*order)[i]];
      strides_out[i] = 1;    
      for (unsigned int j = 0; j < (*order)[i]; j++) {
        strides_out[i] *= (*in->get_dimensions())[j];
      }
    }

    unsigned int* dims_dev        = 0;
    unsigned int* strides_out_dev = 0;

    if (cudaMalloc((void**) &dims_dev, in->get_number_of_dimensions()*sizeof(unsigned int)) != cudaSuccess) {
      throw cuda_error("cuNDArray_permute : Error allocating CUDA dims memory");
    }

    if (cudaMalloc((void**) &strides_out_dev, in->get_number_of_dimensions()*sizeof(unsigned int)) != cudaSuccess) {
      throw cuda_error("cuNDArray_permute : Error allocating CUDA strides_out memory");
    }

    if (cudaMemcpy(dims_dev, dims, in->get_number_of_dimensions()*sizeof(unsigned int), cudaMemcpyHostToDevice) != cudaSuccess) {
      err = cudaGetLastError();
      std::stringstream ss;
      ss << "cuNDArray_permute : Error uploading dimensions to device, " << cudaGetErrorString(err);
      throw cuda_error(ss.str());
    }

    if (cudaMemcpy(strides_out_dev, strides_out, in->get_number_of_dimensions()*sizeof(unsigned int), cudaMemcpyHostToDevice) != cudaSuccess) {
      throw cuda_error("cuNDArray_permute : Error uploading strides to device");
    }

    dim3 blockDim(512,1,1);
    dim3 gridDim;
    if( in->get_number_of_dimensions() > 2 ){
      gridDim = dim3((unsigned int) std::ceil((double)in->get_size(0)*in->get_size(1)/blockDim.x), 1, 1 );
      for( unsigned int d=2; d<in->get_number_of_dimensions(); d++ )
        gridDim.y *= in->get_size(d);
    }
    else
      gridDim = dim3((unsigned int) std::ceil((double)in->get_number_of_elements()/blockDim.x), 1, 1 );

    cuNDArray_permute_kernel<<< gridDim, blockDim >>>( in_ptr, out_ptr, in->get_number_of_dimensions(), 
                                                       dims_dev, strides_out_dev, in->get_number_of_elements(), shift_mode);

    err = cudaGetLastError();
    if( err != cudaSuccess ){
      std::stringstream ss;
      ss <<"cuNDArray_permute : Error during kernel call: " << cudaGetErrorString(err);
      throw cuda_error(ss.str());
    }

    if (cudaFree(dims_dev) != cudaSuccess) {
      err = cudaGetLastError();
      std::stringstream ss;
      ss << "cuNDArray_permute: failed to delete device memory (dims_dev) " << cudaGetErrorString(err);
      throw cuda_error(ss.str());
    }

    if (cudaFree(strides_out_dev) != cudaSuccess) {
      err = cudaGetLastError();
      std::stringstream ss;
      ss << "cuNDArray_permute: failed to delete device memory (strides_out_dev) "<< cudaGetErrorString(err);
      throw cuda_error(ss.str());
    }    
    delete [] dims;
    delete [] strides_out;    
  }  

  template <class T> boost::shared_ptr< cuNDArray<T> >
  permute( cuNDArray<T> *in, std::vector<size_t> *dim_order, int shift_mode )
  {
    if( in == 0x0 || dim_order == 0x0 ) {
      throw std::runtime_error("permute(): invalid pointer provided");
    }    

    std::vector<size_t> dims;
    for (size_t i = 0; i < dim_order->size(); i++)
      dims.push_back(in->get_size(dim_order->at(i)));

    boost::shared_ptr< cuNDArray<T> > out( new cuNDArray<T>(dims) );
    permute( in, out.get(), dim_order, shift_mode );
    return out;
  }

  template <class T> void
  permute( cuNDArray<T> *in, cuNDArray<T> *out, std::vector<size_t> *dim_order, int shift_mode )
  {
    if( in == 0x0 || out == 0x0 || dim_order == 0x0 ) {
      throw std::runtime_error("permute(): invalid pointer provided");
    }    
    if (out->get_number_of_dimensions() != in->get_number_of_dimensions() || out->get_number_of_elements() != in->get_number_of_elements()){
    	throw std::runtime_error("permute(): Input and output have differing dimensions and/or differing number of elements");
    }

    if( in == out ){
      throw std::runtime_error("permute(): in-place permutation not supported");
    }   

    //Check ordering array
    if (dim_order->size() > in->get_number_of_dimensions()) {
      throw std::runtime_error("permute(): invalid length of dimension ordering array");
    }

    std::vector<size_t> dim_count(in->get_number_of_dimensions(),0);
    for (unsigned int i = 0; i < dim_order->size(); i++) {
      if ((*dim_order)[i] >= in->get_number_of_dimensions()) {
        throw std::runtime_error("permute(): invalid dimension order array");
      }
      dim_count[(*dim_order)[i]]++;
    }

    //Create an internal array to store the dimensions
    std::vector<size_t> dim_order_int;

    //Check that there are no duplicate dimensions
    for (unsigned int i = 0; i < dim_order->size(); i++) {
      if (dim_count[(*dim_order)[i]] != 1) {
        throw std::runtime_error("permute(): invalid dimension order array (duplicates)");
      }
      dim_order_int.push_back((*dim_order)[i]);
    }

    for (unsigned int i = 0; i < dim_order_int.size(); i++) {
      if ((*in->get_dimensions())[dim_order_int[i]] != out->get_size(i)) {
        throw std::runtime_error("permute(): dimensions of output array do not match the input array");
      }
    }

    //Pad dimension order array with dimension not mentioned in order array
    if (dim_order_int.size() < in->get_number_of_dimensions()) {
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
    		*out = *in;
    		return;
    	}
    }
    cuNDArray_permute(in, out, &dim_order_int, shift_mode);
  }

  template<class T> boost::shared_ptr< cuNDArray<T> >
  shift_dim( cuNDArray<T> *in, int shift )
  {
    if( in == 0x0 ) {
      throw std::runtime_error("shift_dim(): invalid input pointer provided");
    }    

    std::vector<size_t> order;
    for (int i = 0; i < in->get_number_of_dimensions(); i++) {
      order.push_back(static_cast<unsigned int>((i+shift)%in->get_number_of_dimensions()));
    }
    return permute(in, &order);
  }

  template<class T> 
  void shift_dim( cuNDArray<T> *in, cuNDArray<T> *out, int shift )
  {
    if( in == 0x0 || out == 0x0 ) {
      throw std::runtime_error("shift_dim(): invalid pointer provided");
    }    

    std::vector<size_t> order;
    for (int i = 0; i < in->get_number_of_dimensions(); i++) {
      order.push_back(static_cast<unsigned int>((i+shift)%in->get_number_of_dimensions()));
    }
    permute(in,out,&order);
  }

  // Expand
  //
  template<class T> 
  __global__ void expand_kernel( 
                                const T * __restrict__ in, T * __restrict__ out,
                                unsigned int number_of_elements_in, unsigned int number_of_elements_out, unsigned int new_dim_size )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;    
    if( idx < number_of_elements_out ){
      out[idx] = in[idx%number_of_elements_in];
    }
  }

  // Expand
  //
  template<class T> boost::shared_ptr< cuNDArray<T> > 
  expand( cuNDArray<T> *in, size_t new_dim_size )
  {
    unsigned int number_of_elements_out = in->get_number_of_elements()*new_dim_size;

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( number_of_elements_out, &blockDim, &gridDim );

    // Find element stride
    std::vector<size_t> dims = *in->get_dimensions();
    dims.push_back(new_dim_size);

    // Invoke kernel
    boost::shared_ptr< cuNDArray<T> > out( new cuNDArray<T>());
    out->create(&dims);

    expand_kernel<T><<< gridDim, blockDim >>>( in->get_data_ptr(), out->get_data_ptr(), 
                                               in->get_number_of_elements(), number_of_elements_out, new_dim_size );

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
      const typename uintd<D>::Type co = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
      const typename uintd<D>::Type co_os = offset + co;
      const unsigned int in_idx = co_to_idx<D>(co_os, matrix_size_in)+frame_offset*prod(matrix_size_in);
      out[idx] = in[in_idx];
    }
  }

  // Crop
  template<class T, unsigned int D>
  void crop( typename uint64d<D>::Type offset, cuNDArray<T> *in, cuNDArray<T> *out )
  {
    if( in == 0x0 || out == 0x0 ){
      throw std::runtime_error("crop: 0x0 ndarray provided");
    }

    if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
      throw std::runtime_error("crop: image dimensions mismatch");
    }

    if( in->get_number_of_dimensions() < D ){
      std::stringstream ss;
      ss << "crop: number of image dimensions should be at least " << D;
      throw std::runtime_error(ss.str());
    }

    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *in->get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *out->get_dimensions() );

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
        number_of_batches *= in->get_size(d);
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
           in->get_data_ptr(), out->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
    
    CHECK_FOR_CUDA_ERROR();
  }

  template<class T, unsigned int D> boost::shared_ptr< cuNDArray<T> > 
  crop( typename uint64d<D>::Type offset, typename uint64d<D>::Type size, cuNDArray<T> *in )
  {
    if( in == 0x0 ){
      throw std::runtime_error("crop: 0x0 array provided");
    }
    std::vector<size_t> dims = to_std_vector(size);
    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
      dims.push_back(in->get_size(d));
    }
    boost::shared_ptr< cuNDArray<T> > result( new cuNDArray<T>(&dims) );
    crop<T,D>(offset, in, result.get());
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

      const typename uintd<D>::Type co_out = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
      typename uintd<D>::Type offset;
      for(unsigned int d=0; d<D; d++)
      {
          offset[d] = matrix_size_out[d]/2 - matrix_size_in[d]/2;
      }

      T _out;
      bool inside = (co_out>=offset) && (co_out<(matrix_size_in+offset));

      if( inside )
        _out = in[co_to_idx<D>(co_out-offset, matrix_size_in)+frame_offset*prod(matrix_size_in)];
      else{      
        _out = val;
      }

      out[idx] = _out;
    }
  }

  template<class T, unsigned int D> 
  void pad( cuNDArray<T> *in, cuNDArray<T> *out, T val )
  { 
    if( in == 0x0 || out == 0x0 ){
      throw std::runtime_error("pad: 0x0 ndarray provided");
    }

    if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
      throw std::runtime_error("pad: image dimensions mismatch");
    }

    if( in->get_number_of_dimensions() < D ){
      std::stringstream ss;
      ss << "pad: number of image dimensions should be at least " << D;
      throw std::runtime_error(ss.str());
    }

    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *in->get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *out->get_dimensions() );

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
      number_of_batches *= in->get_size(d);
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
        in->get_data_ptr(), out->get_data_ptr(), number_of_batches, prod(matrix_size_out), val );

    CHECK_FOR_CUDA_ERROR();
  }

  template<class T, unsigned int D> boost::shared_ptr< cuNDArray<T> >
  pad( typename uint64d<D>::Type size, cuNDArray<T> *in, T val )
  {
    if( in == 0x0 ){
      throw std::runtime_error("pad: 0x0 array provided");
    }
    std::vector<size_t> dims = to_std_vector(size);
    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
      dims.push_back(in->get_size(d));
    }
    boost::shared_ptr< cuNDArray<T> > result( new cuNDArray<T>(&dims) );
    pad<T,D>(in, result.get(), val);
    return result;
  }

  template<class T, unsigned int D> 
  __global__ void fill_border_kernel( vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
                                      T *image, unsigned int number_of_batches, unsigned int number_of_elements, T val )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

    if( idx < number_of_elements ){
      const vector_td<unsigned int,D> co_out = idx_to_co<D>( idx, matrix_size_out );
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
  void fill_border( typename uint64d<D>::Type matrix_size_in, cuNDArray<T> *in_out, T val )
  { 
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *in_out->get_dimensions() );

    if( weak_greater(matrix_size_in, matrix_size_out) ){
      throw std::runtime_error("fill_border: size mismatch, cannot zero fill");
    }

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in_out->get_number_of_dimensions(); d++ ){
      number_of_batches *= in_out->get_size(d);
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim );

    // Invoke kernel
    fill_border_kernel<T,D><<< gridDim, blockDim >>>
      ( vector_td<unsigned int,D>(matrix_size_in), vector_td<unsigned int,D>(matrix_size_out),
        in_out->get_data_ptr(), number_of_batches, prod(matrix_size_out), val );

    CHECK_FOR_CUDA_ERROR();
  }


  template<class T, unsigned int D>
  __global__ void fill_border_kernel( typename realType<T>::Type radius, vector_td<int,D> matrix_size,
                                      T *image, unsigned int number_of_batches, unsigned int number_of_elements, T val )
  {
    const int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

    if( idx < number_of_elements ){
      const vector_td<typename realType<T>::Type,D> co_out( (matrix_size>>1) - idx_to_co<D>( idx, matrix_size ));
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
  void fill_border( typename realType<T>::Type radius, cuNDArray<T> *in_out, T val )
  {
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *in_out->get_dimensions() );


    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in_out->get_number_of_dimensions(); d++ ){
      number_of_batches *= in_out->get_size(d);
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim );

    // Invoke kernel
    fill_border_kernel<T,D><<< gridDim, blockDim >>>
      (radius, vector_td<int,D>(matrix_size_out),
        in_out->get_data_ptr(), number_of_batches, prod(matrix_size_out), val );

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
      
      const typename uintd<D>::Type co_out = idx_to_co<D>( idx-batch*num_elements_out, matrix_size_out );
      const typename uintd<D>::Type co_in = co_out >> 1;
      const typename uintd<D>::Type ones(1);
      const typename uintd<D>::Type twos(2);
      const typename uintd<D>::Type offset = co_out%twos;
      
      const unsigned int num_cells = 1 << D;
      
      T cellsum(0);
      unsigned int count = 0;
      
      for( unsigned int i=0; i<num_cells; i++ ){
        
        const typename uintd<D>::Type stride = idx_to_co<D>( i, twos );
        
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
  
  template<class T, unsigned int D> boost::shared_ptr< cuNDArray<T> > upsample( cuNDArray<T>* in )
	{
    if( in == 0x0 )
      throw std::runtime_error("upsample: illegal input pointer");

    std::vector<size_t> dims_out = *in->get_dimensions();
    for( unsigned int i=0; i<D; i++ ) dims_out[i] <<= 1;
    boost::shared_ptr< cuNDArray<T> > out(new cuNDArray<T>(&dims_out));
    upsample<T,D>( in, out.get() );
    return out;
	}

  template<class T, unsigned int D> void upsample( cuNDArray<T> *in, cuNDArray<T> *out )
  {
    if( in == 0x0 || out == 0x0 )
      throw std::runtime_error("upsample: illegal input pointer");

    typename uint64d<D>::Type matrix_size_in  = from_std_vector<size_t,D>( *in->get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *out->get_dimensions() );

    if( (matrix_size_in<<1) != matrix_size_out ){
      throw std::runtime_error("upsample: arrays do not correspond to upsampling by a factor of two");
    }

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<out->get_number_of_dimensions(); d++ ){
      number_of_batches *= out->get_size(d);
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

    // Invoke kernel
    upsample_kernel<T,D><<< gridDim, blockDim >>>
      ( vector_td<unsigned int,D>(matrix_size_in), vector_td<unsigned int,D>(matrix_size_out),
        number_of_batches, in->get_data_ptr(), out->get_data_ptr() );

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
      
      const typename intd<D>::Type co_out = idx_to_co<D>( idx-batch*num_elements_out, matrix_size_out );
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
        const typename intd<D>::Type stride = idx_to_co<D>(i,threes)-ones; // in the range [-1;1]^D
        
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

  template<class T, unsigned int D> boost::shared_ptr< cuNDArray<T> > downsample( cuNDArray<T>* in )
  {
    if( in == 0x0 )
      throw std::runtime_error("downsample: illegal input pointer");
    
    std::vector<size_t> dims_out = *in->get_dimensions();
    for( unsigned int i=0; i<D; i++ ) dims_out[i] >>= 1;
    boost::shared_ptr< cuNDArray<T> > out(new cuNDArray<T>(&dims_out));
    downsample<T,D>( in, out.get() );
    return out;
  }

  template<class T, unsigned int D> void downsample( cuNDArray<T> *in, cuNDArray<T> *out )
  {
    if( in == 0x0 || out == 0x0 )
      throw std::runtime_error("downsample: illegal input pointer");

    typename uint64d<D>::Type matrix_size_in  = from_std_vector<size_t,D>( *in->get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = from_std_vector<size_t,D>( *out->get_dimensions() );

    if( (matrix_size_in>>1) != matrix_size_out ){
      throw std::runtime_error("downsample: arrays do not correspond to downsampling by a factor of two");
    }

    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<out->get_number_of_dimensions(); d++ ){
      number_of_batches *= out->get_size(d);
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

    // Invoke kernel
    downsample_kernel<T,D><<< gridDim, blockDim >>>
      ( vector_td<int,D>(matrix_size_in), vector_td<int,D>(matrix_size_out),
        (int)number_of_batches, in->get_data_ptr(), out->get_data_ptr() );

    CHECK_FOR_CUDA_ERROR();    
  }

  //
  // Instantiation
  //

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > permute( cuNDArray<float>*, std::vector<size_t>*, int );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > permute( cuNDArray<double>*, std::vector<size_t>*, int );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > permute( cuNDArray<float_complext>*, std::vector<size_t>*, int );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > permute( cuNDArray<double_complext>*, std::vector<size_t>*, int );  

  template EXPORTGPUCORE void permute( cuNDArray<float>*, cuNDArray<float>*, std::vector<size_t>*, int);
  template EXPORTGPUCORE void permute( cuNDArray<double>*, cuNDArray<double>*, std::vector<size_t>*, int);
  template EXPORTGPUCORE void permute( cuNDArray<float_complext>*, cuNDArray<float_complext>*, std::vector<size_t>*, int);
  template EXPORTGPUCORE void permute( cuNDArray<double_complext>*, cuNDArray<double_complext>*, std::vector<size_t>*, int);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > shift_dim( cuNDArray<float>*, int );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > shift_dim( cuNDArray<double>*, int );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > shift_dim( cuNDArray<float_complext>*, int );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > shift_dim( cuNDArray<double_complext>*, int );

  template EXPORTGPUCORE void shift_dim( cuNDArray<float>*, cuNDArray<float>*, int shift );
  template EXPORTGPUCORE void shift_dim( cuNDArray<double>*, cuNDArray<double>*, int shift );
  template EXPORTGPUCORE void shift_dim( cuNDArray<float_complext>*, cuNDArray<float_complext>*, int shift );
  template EXPORTGPUCORE void shift_dim( cuNDArray<double_complext>*, cuNDArray<double_complext>*, int shift );

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > expand<float>( cuNDArray<float>*, size_t);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > expand<double>( cuNDArray<double>*, size_t);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > expand<float_complext>( cuNDArray<float_complext>*, size_t);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > expand<double_complext>( cuNDArray<double_complext>*, size_t);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > crop<float,1>( typename uint64d<1>::Type, typename uint64d<1>::Type, cuNDArray<float>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > crop<float,2>( typename uint64d<2>::Type, typename uint64d<2>::Type, cuNDArray<float>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > crop<float,3>( typename uint64d<3>::Type, typename uint64d<3>::Type, cuNDArray<float>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > crop<float,4>( typename uint64d<4>::Type, typename uint64d<4>::Type, cuNDArray<float>*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > crop<float_complext,1>( typename uint64d<1>::Type, typename uint64d<1>::Type, cuNDArray<float_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > crop<float_complext,2>( typename uint64d<2>::Type, typename uint64d<2>::Type, cuNDArray<float_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > crop<float_complext,3>( typename uint64d<3>::Type, typename uint64d<3>::Type, cuNDArray<float_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > crop<float_complext,4>( typename uint64d<4>::Type, typename uint64d<4>::Type, cuNDArray<float_complext>*);

  template EXPORTGPUCORE void crop<float,1>( uint64d1, cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUCORE void crop<float,2>( uint64d2, cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUCORE void crop<float,3>( uint64d3, cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUCORE void crop<float,4>( uint64d4, cuNDArray<float>*, cuNDArray<float>*);

  template EXPORTGPUCORE void crop<complext<float>,1>( uint64d1, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
  template EXPORTGPUCORE void crop<complext<float>,2>( uint64d2, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
  template EXPORTGPUCORE void crop<complext<float>,3>( uint64d3, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
  template EXPORTGPUCORE void crop<complext<float>,4>( uint64d4, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > pad<float,1>( typename uint64d<1>::Type, cuNDArray<float>*, float );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > pad<float,2>( typename uint64d<2>::Type, cuNDArray<float>*, float );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > pad<float,3>( typename uint64d<3>::Type, cuNDArray<float>*, float );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > pad<float,4>( typename uint64d<4>::Type, cuNDArray<float>*, float );

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > pad<float_complext,1>( typename uint64d<1>::Type, cuNDArray<float_complext>*, float_complext );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > pad<float_complext,2>( typename uint64d<2>::Type, cuNDArray<float_complext>*, float_complext );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > pad<float_complext,3>( typename uint64d<3>::Type, cuNDArray<float_complext>*, float_complext );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > pad<float_complext,4>( typename uint64d<4>::Type, cuNDArray<float_complext>*, float_complext );

  template EXPORTGPUCORE void pad<float,1>( cuNDArray<float>*, cuNDArray<float>*, float);
  template EXPORTGPUCORE void pad<float,2>( cuNDArray<float>*, cuNDArray<float>*, float);
  template EXPORTGPUCORE void pad<float,3>( cuNDArray<float>*, cuNDArray<float>*, float);
  template EXPORTGPUCORE void pad<float,4>( cuNDArray<float>*, cuNDArray<float>*, float);

  template EXPORTGPUCORE void pad<float_complext,1>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, float_complext);
  template EXPORTGPUCORE void pad<float_complext,2>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, float_complext);  
  template EXPORTGPUCORE void pad<float_complext,3>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, float_complext);
  template EXPORTGPUCORE void pad<float_complext,4>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, float_complext);

  template EXPORTGPUCORE void fill_border<float,1>(uint64d1, cuNDArray<float>*,float);
  template EXPORTGPUCORE void fill_border<float,2>(uint64d2, cuNDArray<float>*,float);
  template EXPORTGPUCORE void fill_border<float,3>(uint64d3, cuNDArray<float>*,float);
  template EXPORTGPUCORE void fill_border<float,4>(uint64d4, cuNDArray<float>*,float);
  template EXPORTGPUCORE void fill_border<float,1>(float, cuNDArray<float>*,float);
	template EXPORTGPUCORE void fill_border<float,2>(float, cuNDArray<float>*,float);
	template EXPORTGPUCORE void fill_border<float,3>(float, cuNDArray<float>*,float);
	template EXPORTGPUCORE void fill_border<float,4>(float, cuNDArray<float>*,float);

  template EXPORTGPUCORE void fill_border<float_complext,1>(uint64d1, cuNDArray<float_complext>*,float_complext);
  template EXPORTGPUCORE void fill_border<float_complext,2>(uint64d2, cuNDArray<float_complext>*,float_complext);
  template EXPORTGPUCORE void fill_border<float_complext,3>(uint64d3, cuNDArray<float_complext>*,float_complext);
  template EXPORTGPUCORE void fill_border<float_complext,4>(uint64d4, cuNDArray<float_complext>*,float_complext);
  template EXPORTGPUCORE void fill_border<float_complext,1>(float, cuNDArray<float_complext>*,float_complext);
	template EXPORTGPUCORE void fill_border<float_complext,2>(float, cuNDArray<float_complext>*,float_complext);
	template EXPORTGPUCORE void fill_border<float_complext,3>(float, cuNDArray<float_complext>*,float_complext);
	template EXPORTGPUCORE void fill_border<float_complext,4>(float, cuNDArray<float_complext>*,float_complext);


  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > crop<double,1>( typename uint64d<1>::Type, typename uint64d<1>::Type, cuNDArray<double>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > crop<double,2>( typename uint64d<2>::Type, typename uint64d<2>::Type, cuNDArray<double>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > crop<double,3>( typename uint64d<3>::Type, typename uint64d<3>::Type, cuNDArray<double>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > crop<double,4>( typename uint64d<4>::Type, typename uint64d<4>::Type, cuNDArray<double>*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > crop<double_complext,1>( typename uint64d<1>::Type, typename uint64d<1>::Type, cuNDArray<double_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > crop<double_complext,2>( typename uint64d<2>::Type, typename uint64d<2>::Type, cuNDArray<double_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > crop<double_complext,3>( typename uint64d<3>::Type, typename uint64d<3>::Type, cuNDArray<double_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > crop<double_complext,4>( typename uint64d<4>::Type, typename uint64d<4>::Type, cuNDArray<double_complext>*);

  template EXPORTGPUCORE void crop<double,1>( uint64d1, cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUCORE void crop<double,2>( uint64d2, cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUCORE void crop<double,3>( uint64d3, cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUCORE void crop<double,4>( uint64d4, cuNDArray<double>*, cuNDArray<double>*);

  template EXPORTGPUCORE void crop<complext<double>,1>( uint64d1, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
  template EXPORTGPUCORE void crop<complext<double>,2>( uint64d2, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
  template EXPORTGPUCORE void crop<complext<double>,3>( uint64d3, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
  template EXPORTGPUCORE void crop<complext<double>,4>( uint64d4, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > pad<double,1>( typename uint64d<1>::Type, cuNDArray<double>*, double );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > pad<double,2>( typename uint64d<2>::Type, cuNDArray<double>*, double );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > pad<double,3>( typename uint64d<3>::Type, cuNDArray<double>*, double );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > pad<double,4>( typename uint64d<4>::Type, cuNDArray<double>*, double );

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > pad<double_complext,1>( typename uint64d<1>::Type, cuNDArray<double_complext>*, double_complext );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > pad<double_complext,2>( typename uint64d<2>::Type, cuNDArray<double_complext>*, double_complext );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > pad<double_complext,3>( typename uint64d<3>::Type, cuNDArray<double_complext>*, double_complext );
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > pad<double_complext,4>( typename uint64d<4>::Type, cuNDArray<double_complext>*, double_complext );

  template EXPORTGPUCORE void pad<double,1>( cuNDArray<double>*, cuNDArray<double>*, double);
  template EXPORTGPUCORE void pad<double,2>( cuNDArray<double>*, cuNDArray<double>*, double);
  template EXPORTGPUCORE void pad<double,3>( cuNDArray<double>*, cuNDArray<double>*, double);
  template EXPORTGPUCORE void pad<double,4>( cuNDArray<double>*, cuNDArray<double>*, double);

  template EXPORTGPUCORE void pad<double_complext,1>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, double_complext);
  template EXPORTGPUCORE void pad<double_complext,2>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, double_complext);  
  template EXPORTGPUCORE void pad<double_complext,3>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, double_complext);
  template EXPORTGPUCORE void pad<double_complext,4>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, double_complext);

  template EXPORTGPUCORE void fill_border<double,1>(uint64d1, cuNDArray<double>*,double);
  template EXPORTGPUCORE void fill_border<double,2>(uint64d2, cuNDArray<double>*,double);
  template EXPORTGPUCORE void fill_border<double,3>(uint64d3, cuNDArray<double>*,double);
  template EXPORTGPUCORE void fill_border<double,4>(uint64d4, cuNDArray<double>*,double);
  template EXPORTGPUCORE void fill_border<double,1>(double, cuNDArray<double>*,double);
	template EXPORTGPUCORE void fill_border<double,2>(double, cuNDArray<double>*,double);
	template EXPORTGPUCORE void fill_border<double,3>(double, cuNDArray<double>*,double);
	template EXPORTGPUCORE void fill_border<double,4>(double, cuNDArray<double>*,double);

  template EXPORTGPUCORE void fill_border<double_complext,1>(uint64d1, cuNDArray<double_complext>*,double_complext);
  template EXPORTGPUCORE void fill_border<double_complext,2>(uint64d2, cuNDArray<double_complext>*,double_complext);
  template EXPORTGPUCORE void fill_border<double_complext,3>(uint64d3, cuNDArray<double_complext>*,double_complext);
  template EXPORTGPUCORE void fill_border<double_complext,4>(uint64d4, cuNDArray<double_complext>*,double_complext);
  template EXPORTGPUCORE void fill_border<double_complext,1>(double, cuNDArray<double_complext>*,double_complext);
	template EXPORTGPUCORE void fill_border<double_complext,2>(double, cuNDArray<double_complext>*,double_complext);
	template EXPORTGPUCORE void fill_border<double_complext,3>(double, cuNDArray<double_complext>*,double_complext);
	template EXPORTGPUCORE void fill_border<double_complext,4>(double, cuNDArray<double_complext>*,double_complext);


  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > upsample<float,1>(cuNDArray<float>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > upsample<float,2>(cuNDArray<float>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > upsample<float,3>(cuNDArray<float>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > upsample<float,4>(cuNDArray<float>*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > upsample<float_complext,1>(cuNDArray<float_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > upsample<float_complext,2>(cuNDArray<float_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > upsample<float_complext,3>(cuNDArray<float_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > upsample<float_complext,4>(cuNDArray<float_complext>*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > upsample<double,1>(cuNDArray<double>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > upsample<double,2>(cuNDArray<double>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > upsample<double,3>(cuNDArray<double>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > upsample<double,4>(cuNDArray<double>*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > upsample<double_complext,1>(cuNDArray<double_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > upsample<double_complext,2>(cuNDArray<double_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > upsample<double_complext,3>(cuNDArray<double_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > upsample<double_complext,4>(cuNDArray<double_complext>*);

  template EXPORTGPUCORE void upsample<float,1>(cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUCORE void upsample<float,2>(cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUCORE void upsample<float,3>(cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUCORE void upsample<float,4>(cuNDArray<float>*, cuNDArray<float>*);

  template EXPORTGPUCORE void upsample<float_complext,1>(cuNDArray<float_complext>*, cuNDArray<float_complext>*);
  template EXPORTGPUCORE void upsample<float_complext,2>(cuNDArray<float_complext>*, cuNDArray<float_complext>*);
  template EXPORTGPUCORE void upsample<float_complext,3>(cuNDArray<float_complext>*, cuNDArray<float_complext>*);
  template EXPORTGPUCORE void upsample<float_complext,4>(cuNDArray<float_complext>*, cuNDArray<float_complext>*);

  template EXPORTGPUCORE void upsample<double,1>(cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUCORE void upsample<double,2>(cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUCORE void upsample<double,3>(cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUCORE void upsample<double,4>(cuNDArray<double>*, cuNDArray<double>*);

  template EXPORTGPUCORE void upsample<double_complext,1>(cuNDArray<double_complext>*, cuNDArray<double_complext>*);
  template EXPORTGPUCORE void upsample<double_complext,2>(cuNDArray<double_complext>*, cuNDArray<double_complext>*);
  template EXPORTGPUCORE void upsample<double_complext,3>(cuNDArray<double_complext>*, cuNDArray<double_complext>*);
  template EXPORTGPUCORE void upsample<double_complext,4>(cuNDArray<double_complext>*, cuNDArray<double_complext>*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > downsample<float,1>(cuNDArray<float>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > downsample<float,2>(cuNDArray<float>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > downsample<float,3>(cuNDArray<float>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> > downsample<float,4>(cuNDArray<float>*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > downsample<float_complext,1>(cuNDArray<float_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > downsample<float_complext,2>(cuNDArray<float_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > downsample<float_complext,3>(cuNDArray<float_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> > downsample<float_complext,4>(cuNDArray<float_complext>*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > downsample<double,1>(cuNDArray<double>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > downsample<double,2>(cuNDArray<double>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > downsample<double,3>(cuNDArray<double>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> > downsample<double,4>(cuNDArray<double>*);

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > downsample<double_complext,1>(cuNDArray<double_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > downsample<double_complext,2>(cuNDArray<double_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > downsample<double_complext,3>(cuNDArray<double_complext>*);
  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> > downsample<double_complext,4>(cuNDArray<double_complext>*);

  template EXPORTGPUCORE void downsample<float,1>(cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUCORE void downsample<float,2>(cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUCORE void downsample<float,3>(cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUCORE void downsample<float,4>(cuNDArray<float>*, cuNDArray<float>*);

  template EXPORTGPUCORE void downsample<float_complext,1>(cuNDArray<float_complext>*, cuNDArray<float_complext>*);
  template EXPORTGPUCORE void downsample<float_complext,2>(cuNDArray<float_complext>*, cuNDArray<float_complext>*);
  template EXPORTGPUCORE void downsample<float_complext,3>(cuNDArray<float_complext>*, cuNDArray<float_complext>*);
  template EXPORTGPUCORE void downsample<float_complext,4>(cuNDArray<float_complext>*, cuNDArray<float_complext>*);

  template EXPORTGPUCORE void downsample<double,1>(cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUCORE void downsample<double,2>(cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUCORE void downsample<double,3>(cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUCORE void downsample<double,4>(cuNDArray<double>*, cuNDArray<double>*);

  template EXPORTGPUCORE void downsample<double_complext,1>(cuNDArray<double_complext>*, cuNDArray<double_complext>*);
  template EXPORTGPUCORE void downsample<double_complext,2>(cuNDArray<double_complext>*, cuNDArray<double_complext>*);
  template EXPORTGPUCORE void downsample<double_complext,3>(cuNDArray<double_complext>*, cuNDArray<double_complext>*);
  template EXPORTGPUCORE void downsample<double_complext,4>(cuNDArray<double_complext>*, cuNDArray<double_complext>*);


  // We can probably instantiate the functions below functionsfor many more types? E.g. arrays of floatd2. 
  // For now we just introduce what we have needed...
  //

  template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd2> > expand<floatd2>( cuNDArray<floatd2>*, size_t);  
}
