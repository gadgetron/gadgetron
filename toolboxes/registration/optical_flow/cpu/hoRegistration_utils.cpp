#include "hoRegistration_utils.h"
#include "vector_td_utilities.h"

#ifdef USE_OMP
#include <omp.h>
#endif 

namespace Gadgetron{

  // Utility to check if all neighbors required for the linear interpolation exists
  // ... do not include dimensions of size 1

  template<class REAL, unsigned int D> inline bool
  is_border_pixel( vector_td<size_t,D> co, vector_td<size_t,D> dims )
  {
    for( size_t dim=0; dim<D; dim++ ){
      if( dims[dim] > 1 && ( co[dim] == 0 || co[dim] == (dims[dim]-1) ) )
	return true;
    }
    return false;
  }

  // Downsample
  template<class REAL, unsigned int D> 
  boost::shared_ptr< hoNDArray<REAL> > downsample( hoNDArray<REAL> *_in )
  {
    // A few sanity checks 

    if( _in == 0x0 ){
      throw std::runtime_error( "downsample(): illegal input provided.");
    }
    
    if( _in->get_number_of_dimensions() < D ){
      throw std::runtime_error( "downsample(): the number of array dimensions should be at least D");
    }
    
    for( size_t d=0; d<D; d++ ){
      if( (_in->get_size(d)%2) == 1 && _in->get_size(d) != 1 ){
	throw std::runtime_error( "downsample(): uneven array dimensions larger than one not accepted");
      }
    }
    
    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *_in->get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = matrix_size_in >> 1;

    for( size_t d=0; d<D; d++ ){
      if( matrix_size_out[d] == 0 ) 
	matrix_size_out[d] = 1;
    }
  
    size_t num_elements = prod(matrix_size_out);
    size_t num_batches = 1;

    for( size_t d=D; d<_in->get_number_of_dimensions(); d++ ){
      num_batches *= _in->get_size(d);
    }
  
    std::vector<size_t> dims = to_std_vector(matrix_size_out);
    for( size_t d=D; d<_in->get_number_of_dimensions(); d++ ){
      dims.push_back(_in->get_size(d));
    }
  
    REAL *in = _in->get_data_ptr();

    boost::shared_ptr< hoNDArray<REAL> > _out( new hoNDArray<REAL>(&dims) );
    REAL *out = _out->get_data_ptr();
    
    typedef vector_td<size_t,D> uint64d;

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( long long idx=0; idx < num_elements*num_batches; idx++ ){

      const size_t frame_offset = idx/num_elements;
      const uint64d co_out = idx_to_co<D>( idx-frame_offset*num_elements, matrix_size_out );
      const uint64d co_in = co_out << 1;
      const uint64d twos = to_vector_td<size_t,D>(2);
      const size_t num_adds = 1 << D;

      size_t actual_adds = 0;
      REAL res = REAL(0);

      for( size_t i=0; i<num_adds; i++ ){
	const uint64d local_co = idx_to_co<D>( i, twos );
	if( weak_greater_equal( local_co, matrix_size_out ) ) continue; // To allow array dimensions of size 1
	const size_t in_idx = co_to_idx<D>(co_in+local_co, matrix_size_in)+frame_offset*prod(matrix_size_in);
	actual_adds++;
	res += in[in_idx];
      }    
      out[idx] = res/REAL(actual_adds);
    }

    return _out;
  }

  // Linear interpolation upsampling
  template<class REAL, unsigned int D> boost::shared_ptr< hoNDArray<REAL> >
  upsample( hoNDArray<REAL> *_in )
  {
    // A few sanity checks 

    if( _in == 0x0 ){
      throw std::runtime_error("upsample(): illegal input provided.");
    }

    if( _in->get_number_of_dimensions() < D ){
      throw std::runtime_error( "upsample(): the number of array dimensions should be at least D");
    }
    
    typename uint64d<D>::Type matrix_size_in = from_std_vector<size_t,D>( *_in->get_dimensions() );
    typename uint64d<D>::Type matrix_size_out = matrix_size_in << 1;

    for( size_t d=0; d<D; d++ ){
      if( matrix_size_in[d] == 1 )
	matrix_size_out[d] = 1;
    }
  
    size_t num_elements = prod(matrix_size_out);
    size_t num_batches = 1;

    for( size_t d=D; d<_in->get_number_of_dimensions(); d++ ){
      num_batches *= _in->get_size(d);
    }
  
    std::vector<size_t> dims = to_std_vector(matrix_size_out);
    for( size_t d=D; d<_in->get_number_of_dimensions(); d++ ){
      dims.push_back(_in->get_size(d));
    }

    REAL *in = _in->get_data_ptr();

    boost::shared_ptr< hoNDArray<REAL> > _out( new hoNDArray<REAL>(&dims) );
    REAL *out = _out->get_data_ptr();
    
    typedef vector_td<size_t,D> uint64d;

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for( long long idx=0; idx < num_elements*num_batches; idx++ ){
      
      REAL res = REAL(0);

      const size_t num_neighbors = 1 << D;
      const size_t frame_idx = idx/num_elements;
      const uint64d co_out = idx_to_co<D>( idx-frame_idx*num_elements, matrix_size_out );

      // We will only proceed if all neighbours exist (this adds a zero-boundary to the upsampled image/vector field)
      //
    
      if( !is_border_pixel<REAL,D>(co_out, matrix_size_out) ){
      
	for( size_t i=0; i<num_neighbors; i++ ){
	
	  // Determine coordinate of neighbor in input
	  //

	  const uint64d twos = to_vector_td<size_t,D>(2);
	  const uint64d stride = idx_to_co<D>( i, twos );

	  if( weak_greater_equal( stride, matrix_size_out ) ) continue; // To allow array dimensions of 1

	  // Be careful about dimensions of size 1
	  uint64d ones = to_vector_td<size_t,D>(1);
	  for( size_t d=0; d<D; d++ ){
	    if( matrix_size_out[d] == 1 )
	      ones[d] = 0;
	  }
	  uint64d co_in = ((co_out-ones)>>1)+stride;
	
	  // Read corresponding pixel value
	  //
	
	  const size_t in_idx = co_to_idx<D>(co_in, matrix_size_in)+frame_idx*prod(matrix_size_in);
	  REAL value = in[in_idx];
	
	  // Determine weight
	  //
	
	  REAL weight = REAL(1);
	
	  for( size_t dim=0; dim<D; dim++ ){	  
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
    
    return _out;
  }

  //
  // Instantiation
  //
  
  template EXPORTCPUREG boost::shared_ptr< hoNDArray<float> > downsample<float,1>(hoNDArray<float>*);
  template EXPORTCPUREG boost::shared_ptr< hoNDArray<float> > upsample<float,1>(hoNDArray<float>*);

  template EXPORTCPUREG boost::shared_ptr< hoNDArray<float> > downsample<float,2>(hoNDArray<float>*);
  template EXPORTCPUREG boost::shared_ptr< hoNDArray<float> > upsample<float,2>(hoNDArray<float>*);

  template EXPORTCPUREG boost::shared_ptr< hoNDArray<float> > downsample<float,3>(hoNDArray<float>*);
  template EXPORTCPUREG boost::shared_ptr< hoNDArray<float> > upsample<float,3>(hoNDArray<float>*);

  template EXPORTCPUREG boost::shared_ptr< hoNDArray<float> > downsample<float,4>(hoNDArray<float>*);
  template EXPORTCPUREG boost::shared_ptr< hoNDArray<float> > upsample<float,4>(hoNDArray<float>*);

  template EXPORTCPUREG boost::shared_ptr< hoNDArray<double> > downsample<double,1>(hoNDArray<double>*);
  template EXPORTCPUREG boost::shared_ptr< hoNDArray<double> > upsample<double,1>(hoNDArray<double>*);

  template EXPORTCPUREG boost::shared_ptr< hoNDArray<double> > downsample<double,2>(hoNDArray<double>*);
  template EXPORTCPUREG boost::shared_ptr< hoNDArray<double> > upsample<double,2>(hoNDArray<double>*);

  template EXPORTCPUREG boost::shared_ptr< hoNDArray<double> > downsample<double,3>(hoNDArray<double>*);
  template EXPORTCPUREG boost::shared_ptr< hoNDArray<double> > upsample<double,3>(hoNDArray<double>*);

  template EXPORTCPUREG boost::shared_ptr< hoNDArray<double> > downsample<double,4>(hoNDArray<double>*);
  template EXPORTCPUREG boost::shared_ptr< hoNDArray<double> > upsample<double,4>(hoNDArray<double>*);
}
