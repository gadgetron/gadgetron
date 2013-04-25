#include "NFFT_utils.h"
#include "vector_td_utilities.h"
#include "check_CUDA.h"
#include "cudaDeviceManager.h"

namespace Gadgetron {

  void setup_grid( unsigned int number_of_elements, dim3 *blockDim, dim3* gridDim, unsigned int num_batches = 1 )
  {    
    int cur_device = cudaDeviceManager::Instance()->getCurrentDevice();
    int maxGridDim = cudaDeviceManager::Instance()->max_griddim(cur_device);
    
    // For small arrays we keep the block dimension fairly small
    *blockDim = dim3(256);
    *gridDim = dim3((number_of_elements+blockDim->x-1)/blockDim->x, num_batches);

    // Extend block/grid dimensions for large arrays
    if( gridDim->x > maxGridDim){
      blockDim->x = maxGridDim;
      gridDim->x = (number_of_elements+blockDim->x-1)/blockDim->x;
    }
    
    if( gridDim->x > maxGridDim ){
      gridDim->x = ((unsigned int)std::sqrt((float)number_of_elements)+blockDim->x-1)/blockDim->x;
      gridDim->y *= ((number_of_elements+blockDim->x*gridDim->x-1)/(blockDim->x*gridDim->x));
    }
    
    if( gridDim->x >maxGridDim || gridDim->y >maxGridDim){      
      BOOST_THROW_EXCEPTION(cuda_error("Grid dimension larger than supported by device"));
    }
  }

  // Crop
  template<class T, unsigned int D> __global__ void
  crop_kernel( vector_td<unsigned int,D> offset, vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
	       T *in, T *out, unsigned int num_batches, unsigned int num_elements )
  {
    typedef vector_td<unsigned int,D> uintd;
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const unsigned int frame_offset = idx/num_elements;

    if( idx < num_elements*num_batches ){
      const uintd co = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
      const uintd co_os = offset + co;
      const unsigned int in_idx = co_to_idx<D>(co_os, matrix_size_in)+frame_offset*prod(matrix_size_in);
      out[idx] = in[in_idx];
    }
  }

  // Crop
  template<class T, unsigned int D> EXPORTGPUNFFT
  void crop( typename uintd<D>::Type offset, cuNDArray<T> *in, cuNDArray<T> *out )
  {
    if( in == 0x0 || out == 0x0 ){
      BOOST_THROW_EXCEPTION(runtime_error("crop: 0x0 ndarray provided"));
    }

    if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
      BOOST_THROW_EXCEPTION(runtime_error("crop: image dimensions mismatch"));
    }

    if( in->get_number_of_dimensions() < D ){
      std::stringstream ss;
      ss << "crop: number of image dimensions should be at least " << D;
      BOOST_THROW_EXCEPTION(runtime_error(ss.str()));
    }

    typename uintd<D>::Type matrix_size_in = from_std_vector<unsigned int,D>( *in->get_dimensions() );
    typename uintd<D>::Type matrix_size_out = from_std_vector<unsigned int,D>( *out->get_dimensions() );
 
    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
      number_of_batches *= in->get_size(d);
    }

    if( weak_greater(offset+matrix_size_out, matrix_size_in) ){
      BOOST_THROW_EXCEPTION(runtime_error( "crop: cropping size mismatch"));
    }
  
    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

    // Invoke kernel
    crop_kernel<T,D><<< gridDim, blockDim >>>
      ( offset, matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
    CHECK_FOR_CUDA_ERROR();
  }

  // Expand and zero fill
  template<class T, unsigned int D> __global__ void
  expand_with_zero_fill_kernel( vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
				T *in, T *out, unsigned int number_of_batches, unsigned int num_elements )
  {
    typedef vector_td<unsigned int,D> uintd;
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const unsigned int frame_offset = idx/num_elements;

    if( idx < num_elements*number_of_batches ){

      const uintd co_out = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
      const uintd offset = (matrix_size_out-matrix_size_in)>>1;
      T _out;
      bool inside = (co_out>=offset) && (co_out<(matrix_size_in+offset));

      if( inside )
	_out = in[co_to_idx<D>(co_out-offset, matrix_size_in)+frame_offset*prod(matrix_size_in)];
      else{      
	_out = T(0);
      }

      out[idx] = _out;
    }
  }

  // Expand and zero fill
  template<class T, unsigned int D> 
  void expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out )
  { 
    if( in == 0x0 || out == 0x0 ){
      BOOST_THROW_EXCEPTION(runtime_error("zero_fill: 0x0 ndarray provided"));
    }

    if( in->get_number_of_dimensions() != out->get_number_of_dimensions() ){
      BOOST_THROW_EXCEPTION(runtime_error("zero_fill: image dimensions mismatch"));
    }

    if( in->get_number_of_dimensions() < D ){
      std::stringstream ss;
      ss << "zero_fill: number of image dimensions should be at least " << D;
      BOOST_THROW_EXCEPTION(runtime_error(ss.str()));
    }

    typename uintd<D>::Type matrix_size_in = from_std_vector<unsigned int,D>( *in->get_dimensions() );
    typename uintd<D>::Type matrix_size_out = from_std_vector<unsigned int,D>( *out->get_dimensions() );
  
    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
      number_of_batches *= in->get_size(d);
    }

    if( weak_greater(matrix_size_in,matrix_size_out) ){
      std::runtime_error("expand: size mismatch, cannot expand");
    }
 
    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );
 
    // Invoke kernel
    expand_with_zero_fill_kernel<T,D><<< gridDim, blockDim >>> 
      ( matrix_size_in, matrix_size_out, in->get_data_ptr(), out->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
    CHECK_FOR_CUDA_ERROR();
  }

  // Zero fill border (rectangular)
  template<class T, unsigned int D> __global__ void
  zero_fill_border_kernel( vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
			   T *image, unsigned int number_of_batches, unsigned int number_of_elements )
  {
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

    if( idx < number_of_elements ){
      const vector_td<unsigned int,D> co_out = idx_to_co<D>( idx, matrix_size_out );
      const vector_td<unsigned int,D> offset = (matrix_size_out-matrix_size_in)>>1;
      if( weak_less( co_out, offset ) || weak_greater_equal( co_out, matrix_size_in+offset ) ){
	T zero = T(0);
	for( unsigned int batch=0; batch<number_of_batches; batch++ ){
	  image[idx+batch*number_of_elements] = zero;
	}
      }
      else
	; // do nothing
    }
  }

  // Zero fill border (rectangular)
  template<class T, unsigned int D> 
  void zero_fill_border( typename uintd<D>::Type matrix_size_in, cuNDArray<T> *in_out )
  { 
    typename uintd<D>::Type matrix_size_out = from_std_vector<unsigned int,D>( *in_out->get_dimensions() );
 
    if( weak_greater(matrix_size_in, matrix_size_out) ){
      BOOST_THROW_EXCEPTION(runtime_error("zero_fill: size mismatch, cannot zero fill"));
    }
 
    unsigned int number_of_batches = 1;
    for( unsigned int d=D; d<in_out->get_number_of_dimensions(); d++ ){
      number_of_batches *= in_out->get_size(d);
    }

    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( prod(matrix_size_out), &blockDim, &gridDim );
 
    // Invoke kernel
    zero_fill_border_kernel<T,D><<< gridDim, blockDim >>>
      ( matrix_size_in, matrix_size_out, in_out->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
    CHECK_FOR_CUDA_ERROR();
  }
 
  //
  // Instantiation
  //

  template EXPORTGPUNFFT void
  crop<float,1>( uintd1, cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUNFFT void
  crop<float,2>( uintd2, cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUNFFT void
  crop<float,3>( uintd3, cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUNFFT void
  crop<float,4>( uintd4, cuNDArray<float>*, cuNDArray<float>*);

  template EXPORTGPUNFFT void
  crop<complext<float>,1>( uintd1, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
  template EXPORTGPUNFFT void
  crop<complext<float>,2>( uintd2, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
  template EXPORTGPUNFFT void
  crop<complext<float>,3>( uintd3, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
  template EXPORTGPUNFFT void
  crop<complext<float>,4>( uintd4, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*);
  
  template EXPORTGPUNFFT void
  expand_with_zero_fill<float,1>( cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<float,2>( cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<float,3>( cuNDArray<float>*, cuNDArray<float>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<float,4>( cuNDArray<float>*, cuNDArray<float>*);
  
  template EXPORTGPUNFFT void
  expand_with_zero_fill<float_complext,1>( cuNDArray<float_complext>*, cuNDArray<float_complext>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<float_complext,2>( cuNDArray<float_complext>*, cuNDArray<float_complext>*);  
  template EXPORTGPUNFFT void
  expand_with_zero_fill<float_complext,3>( cuNDArray<float_complext>*, cuNDArray<float_complext>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<float_complext,4>( cuNDArray<float_complext>*, cuNDArray<float_complext>*);
  
  template EXPORTGPUNFFT void zero_fill_border<float,1>(uintd1, cuNDArray<float>*);
  template EXPORTGPUNFFT void zero_fill_border<float,2>(uintd2, cuNDArray<float>*);
  template EXPORTGPUNFFT void zero_fill_border<float,3>(uintd3, cuNDArray<float>*);
  template EXPORTGPUNFFT void zero_fill_border<float,4>(uintd4, cuNDArray<float>*);
  
  template EXPORTGPUNFFT void zero_fill_border<float_complext,1>(uintd1, cuNDArray<float_complext>*);
  template EXPORTGPUNFFT void zero_fill_border<float_complext,2>(uintd2, cuNDArray<float_complext>*);
  template EXPORTGPUNFFT void zero_fill_border<float_complext,3>(uintd3, cuNDArray<float_complext>*);
  template EXPORTGPUNFFT void zero_fill_border<float_complext,4>(uintd4, cuNDArray<float_complext>*);

  template EXPORTGPUNFFT void
  crop<double,1>( uintd1, cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUNFFT void
  crop<double,2>( uintd2, cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUNFFT void
  crop<double,3>( uintd3, cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUNFFT void
  crop<double,4>( uintd4, cuNDArray<double>*, cuNDArray<double>*);

  template EXPORTGPUNFFT void
  crop<complext<double>,1>( uintd1, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
  template EXPORTGPUNFFT void
  crop<complext<double>,2>( uintd2, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
  template EXPORTGPUNFFT void
  crop<complext<double>,3>( uintd3, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
  template EXPORTGPUNFFT void
  crop<complext<double>,4>( uintd4, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*);
  
  template EXPORTGPUNFFT void
  expand_with_zero_fill<double,1>( cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<double,2>( cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<double,3>( cuNDArray<double>*, cuNDArray<double>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<double,4>( cuNDArray<double>*, cuNDArray<double>*);
  
  template EXPORTGPUNFFT void
  expand_with_zero_fill<double_complext,1>( cuNDArray<double_complext>*, cuNDArray<double_complext>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<double_complext,2>( cuNDArray<double_complext>*, cuNDArray<double_complext>*);  
  template EXPORTGPUNFFT void
  expand_with_zero_fill<double_complext,3>( cuNDArray<double_complext>*, cuNDArray<double_complext>*);
  template EXPORTGPUNFFT void
  expand_with_zero_fill<double_complext,4>( cuNDArray<double_complext>*, cuNDArray<double_complext>*);
  
  template EXPORTGPUNFFT void zero_fill_border<double,1>(uintd1, cuNDArray<double>*);
  template EXPORTGPUNFFT void zero_fill_border<double,2>(uintd2, cuNDArray<double>*);
  template EXPORTGPUNFFT void zero_fill_border<double,3>(uintd3, cuNDArray<double>*);
  template EXPORTGPUNFFT void zero_fill_border<double,4>(uintd4, cuNDArray<double>*);
  
  template EXPORTGPUNFFT void zero_fill_border<double_complext,1>(uintd1, cuNDArray<double_complext>*);
  template EXPORTGPUNFFT void zero_fill_border<double_complext,2>(uintd2, cuNDArray<double_complext>*);
  template EXPORTGPUNFFT void zero_fill_border<double_complext,3>(uintd3, cuNDArray<double_complext>*);
  template EXPORTGPUNFFT void zero_fill_border<double_complext,4>(uintd4, cuNDArray<double_complext>*);
}
