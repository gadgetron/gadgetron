#include "cuRegistration_utils.h"
#include "cudaDeviceManager.h"
#include "vector_td_utilities.h"

namespace Gadgetron{

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
      BOOST_THROW_EXCEPTION(runtime_error("Grid dimension larger than supported by device"));
    }
  }

  // Downsample
  template<class REAL, unsigned int D> __global__ void
  downsample_kernel( REAL *in, REAL *out, 
		     vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
		     unsigned int num_elements, unsigned int num_batches )
  {
    // We have started a thread for each output element
    typedef vector_td<unsigned int,D> uintd;
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    const unsigned int frame_offset = idx/num_elements;
  
    if( idx < num_elements*num_batches ){

      const uintd co_out = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
      const uintd co_in = co_out << 1;

      const uintd twos = to_vector_td<unsigned int,D>(2);
      const unsigned int num_adds = 1 << D;
      unsigned int actual_adds = 0;

      REAL res = REAL(0);

      for( unsigned int i=0; i<num_adds; i++ ){
	const uintd local_co = idx_to_co<D>( i, twos );
	if( weak_greater_equal( local_co, matrix_size_out ) ) continue; // To allow array dimensions of size 1
	const unsigned int in_idx = co_to_idx<D>(co_in+local_co, matrix_size_in)+frame_offset*prod(matrix_size_in);
	actual_adds++;
	res += in[in_idx];
      }    
      out[idx] = res/REAL(actual_adds);
    }
  }

  // Downsample
  template<class REAL, unsigned int D> 
  boost::shared_ptr< cuNDArray<REAL> > downsample( cuNDArray<REAL> *in )
  {
    // A few sanity checks 

    if( in == 0x0 ){
      BOOST_THROW_EXCEPTION(runtime_error("downsample(): illegal input provided."));
    }
    
    if( in->get_number_of_dimensions() < D ){
      BOOST_THROW_EXCEPTION(runtime_error( "downsample(): the number of array dimensions should be at least D"));      
    }
    
    for( unsigned int d=0; d<D; d++ ){
      if( (in->get_size(d)%2) == 1 && in->get_size(d) != 1 ){
	BOOST_THROW_EXCEPTION(runtime_error( "downsample(): uneven array dimensions larger than one not accepted"));
      }
    }
    
    typename uintd<D>::Type matrix_size_in = from_std_vector<unsigned int,D>( *in->get_dimensions() );
    typename uintd<D>::Type matrix_size_out = matrix_size_in >> 1;

    for( unsigned int d=0; d<D; d++ ){
      if( matrix_size_out[d] == 0 ) 
	matrix_size_out[d] = 1;
    }
  
    unsigned int number_of_elements = prod(matrix_size_out);
    unsigned int number_of_batches = 1;

    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
      number_of_batches *= in->get_size(d);
    }
  
    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;

    setup_grid( number_of_elements, &blockDim, &gridDim, number_of_batches );
    
    std::vector<unsigned int> dims = to_std_vector(matrix_size_out);
    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
      dims.push_back(in->get_size(d));
    }
  
    boost::shared_ptr< cuNDArray<REAL> > out( new cuNDArray<REAL>(&dims) );
    
    // Invoke kernel
    downsample_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
    
    CHECK_FOR_CUDA_ERROR();

    return out;
  }

  // Utility to check if all neighbors required for the linear interpolation exists
  // ... do not include dimensions of size 1

  template<class REAL, unsigned int D> __device__ 
  bool is_border_pixel( vector_td<unsigned int,D> co, vector_td<unsigned int,D> dims )
  {
    for( unsigned int dim=0; dim<D; dim++ ){
      if( dims[dim] > 1 && ( co[dim] == 0 || co[dim] == (dims[dim]-1) ) )
	return true;
    }
    return false;
  }

  // Linear upsampling
  template<class REAL, unsigned int D> __global__ void
  upsample_lin_kernel( REAL *in, REAL *out,
		       vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
		       unsigned int num_elements, unsigned int num_batches )
  {
    // We have started a thread for each output element
    typedef vector_td<unsigned int,D> uintd;
    const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
    
    if( idx < num_elements*num_batches ){

      REAL res = REAL(0);
      const unsigned int num_neighbors = 1 << D;
      const unsigned int frame_idx = idx/num_elements;
      const uintd co_out = idx_to_co<D>( idx-frame_idx*num_elements, matrix_size_out );

      // We will only proceed if all neighbours exist (this adds a zero-boundary to the upsampled image/vector field)
      //
    
      if( !is_border_pixel<REAL,D>(co_out, matrix_size_out) ){
      
	for( unsigned int i=0; i<num_neighbors; i++ ){
	
	  // Determine coordinate of neighbor in input
	  //

	  const uintd twos = to_vector_td<unsigned int,D>(2);
	  const uintd stride = idx_to_co<D>( i, twos );

	  if( weak_greater_equal( stride, matrix_size_out ) ) continue; // To allow array dimensions of 1

	  // Be careful about dimensions of size 1
	  uintd ones = to_vector_td<unsigned int,D>(1);
	  for( unsigned int d=0; d<D; d++ ){
	    if( matrix_size_out[d] == 1 )
	      ones[d] = 0;
	  }
	  uintd co_in = ((co_out-ones)>>1)+stride;
	
	  // Read corresponding pixel value
	  //
	
	  const unsigned int in_idx = co_to_idx<D>(co_in, matrix_size_in)+frame_idx*prod(matrix_size_in);
	  REAL value = in[in_idx];
	
	  // Determine weight
	  //
	
	  REAL weight = REAL(1);
	
	  for( unsigned int dim=0; dim<D; dim++ ){	  
	    if( matrix_size_in[dim] > 1 ){
	      if( stride.vec[dim] == (co_out.vec[dim]%2) ) {
		weight *= REAL(0.25);
	      }
	      else{
		weight *= REAL(0.75);
	      }
	    }
	  }
	
	  // Accumulate result
	  //
	
	  res += weight*value;
	}
      }
      out[idx] = res;
    }
  }

  // Linear interpolation upsampling
  template<class REAL, unsigned int D> boost::shared_ptr< cuNDArray<REAL> >
  upsample( cuNDArray<REAL> *in )
  {
    // A few sanity checks 

    if( in == 0x0 ){
      BOOST_THROW_EXCEPTION(runtime_error("upsample(): illegal input provided."));
    }

    if( in->get_number_of_dimensions() < D ){
      BOOST_THROW_EXCEPTION(runtime_error( "upsample(): the number of array dimensions should be at least D"));
    }
    
    typename uintd<D>::Type matrix_size_in = from_std_vector<unsigned int,D>( *in->get_dimensions() );
    typename uintd<D>::Type matrix_size_out = matrix_size_in << 1;

    for( unsigned int d=0; d<D; d++ ){
      if( matrix_size_in[d] == 1 )
	matrix_size_out[d] = 1;
    }
  
    unsigned int number_of_elements = prod(matrix_size_out);
    unsigned int number_of_batches = 1;

    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
      number_of_batches *= in->get_size(d);
    }
  
    // Setup block/grid dimensions
    dim3 blockDim; dim3 gridDim;
    setup_grid( number_of_elements, &blockDim, &gridDim, number_of_batches );
    
    std::vector<unsigned int> dims = to_std_vector(matrix_size_out);
    for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
      dims.push_back(in->get_size(d));
    }

    boost::shared_ptr< cuNDArray<REAL> > out( new cuNDArray<REAL>(&dims) );

    // Invoke kernel
    upsample_lin_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
    
    CHECK_FOR_CUDA_ERROR();
    
    return out;
  }

  //
  // Instantiation
  //
  
  template EXPORTGPUREG boost::shared_ptr< cuNDArray<float> > downsample<float,1>(cuNDArray<float>*);
  template EXPORTGPUREG boost::shared_ptr< cuNDArray<float> > upsample<float,1>(cuNDArray<float>*);

  template EXPORTGPUREG boost::shared_ptr< cuNDArray<float> > downsample<float,2>(cuNDArray<float>*);
  template EXPORTGPUREG boost::shared_ptr< cuNDArray<float> > upsample<float,2>(cuNDArray<float>*);

  template EXPORTGPUREG boost::shared_ptr< cuNDArray<float> > downsample<float,3>(cuNDArray<float>*);
  template EXPORTGPUREG boost::shared_ptr< cuNDArray<float> > upsample<float,3>(cuNDArray<float>*);

  template EXPORTGPUREG boost::shared_ptr< cuNDArray<float> > downsample<float,4>(cuNDArray<float>*);
  template EXPORTGPUREG boost::shared_ptr< cuNDArray<float> > upsample<float,4>(cuNDArray<float>*);

  template EXPORTGPUREG boost::shared_ptr< cuNDArray<double> > downsample<double,1>(cuNDArray<double>*);
  template EXPORTGPUREG boost::shared_ptr< cuNDArray<double> > upsample<double,1>(cuNDArray<double>*);

  template EXPORTGPUREG boost::shared_ptr< cuNDArray<double> > downsample<double,2>(cuNDArray<double>*);
  template EXPORTGPUREG boost::shared_ptr< cuNDArray<double> > upsample<double,2>(cuNDArray<double>*);

  template EXPORTGPUREG boost::shared_ptr< cuNDArray<double> > downsample<double,3>(cuNDArray<double>*);
  template EXPORTGPUREG boost::shared_ptr< cuNDArray<double> > upsample<double,3>(cuNDArray<double>*);

  template EXPORTGPUREG boost::shared_ptr< cuNDArray<double> > downsample<double,4>(cuNDArray<double>*);
  template EXPORTGPUREG boost::shared_ptr< cuNDArray<double> > upsample<double,4>(cuNDArray<double>*);
}
