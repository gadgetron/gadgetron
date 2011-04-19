#include "b1_map.hcu"
#include "vectord_operators.hcu"
#include "vectord_utilities.hcu"
#include "ndarray_device_utilities.hcu"
#include "check_CUDA.h"
#include "cuNDFFT.h"

#include <math_functions.h>
#include <iostream>

using namespace std;

template< class REAL, unsigned int D> __host__ 
auto_ptr< cuNDArray<real_complex<REAL> > > set_box_convkernel( uintd<D> dims, uintd<D> box );

template<class REAL> __host__ 
auto_ptr< cuNDArray<real_complex<REAL> > > extract_csm( cuNDArray<real_complex<REAL> > *corrm_in, unsigned int number_of_batches, unsigned int number_of_elements );

template<class REAL> __host__ 
void set_phase_reference( cuNDArray<real_complex<REAL> > *csm, unsigned int number_of_batches, unsigned int number_of_elements );

//
// Main method:
//

template<class REAL, unsigned int D> auto_ptr< cuNDArray<real_complex<REAL> > >
estimate_b1_map( cuNDArray<real_complex<REAL> >* data_in )
{
  if( data_in->get_number_of_dimensions() < 2 ){
    cout << endl << "estimate_b1_map:: dimensionality mismatch." << endl; 
    return auto_ptr< cuNDArray<real_complex<REAL> > >(0x0);
  }

  if( data_in->get_number_of_dimensions()-1 != D ){
    cout << endl << "estimate_b1_map:: dimensionality mismatch." << endl; 
    return auto_ptr< cuNDArray<real_complex<REAL> > >(0x0);
  }

  vector<unsigned int> image_dims, dims_to_xform;
  unsigned int pixels_per_coil = 1;
  
  for( unsigned int i=0; i<D; i++ ){
    image_dims.push_back(data_in->get_size(i));
    dims_to_xform.push_back(i);
    pixels_per_coil *= data_in->get_size(i);
  }
  
  unsigned int ncoils = data_in->get_size(D);

  // Make a copy of input data (use aligned structs)
  cuNDArray<real_complex<REAL> > *_data_out = (sizeof(REAL)==sizeof(float)) ?
    (cuNDArray<real_complex<REAL> >*) new cuNDArray<float_complex>((cuNDArray<float_complex>*) data_in) : 
    (cuNDArray<real_complex<REAL> >*) new cuNDArray<double_complex>((cuNDArray<double_complex>*) data_in);
  auto_ptr< cuNDArray<real_complex<REAL> > > data_out(_data_out);
  
  // Normalize by the RSS of the coils
  if( !cuNDA_rss_normalize<REAL, real_complex<REAL> >( data_out.get(), D ) ){
    cout << endl << "estimate_b1_map:: error in rss_normalize" << endl;
    return auto_ptr< cuNDArray<real_complex<REAL> > >(0x0);
  }
  
  // Now calculate the correlation matrices
  auto_ptr< cuNDArray<real_complex<REAL> > > corrm = cuNDA_correlation( data_out.get() );
  data_out.reset();
  
  // Compute smoothing kernel for convolution
  uintd<D> dims; cuNDA_fromVec<D>( image_dims, dims );

  uintd<D> box; to_vectord<unsigned int,D>(box,7);
  auto_ptr< cuNDArray<real_complex<REAL> > > conv_kernel = set_box_convkernel<REAL,D>( dims, box );

  // Perform convolution by multiplication in image space
  cuNDFFT().fft( (cuNDArray<cuFloatComplex>*) conv_kernel.get() );           // TODO: fixme (requires new cuNDFFT interface)
  cuNDFFT().fft( (cuNDArray<cuFloatComplex>*) corrm.get(), dims_to_xform );  // TODO: fixme (requires new cuNDFFT interface)
  cuNDA_scale( conv_kernel.get(), corrm.get() );
  cuNDFFT().ifft( (cuNDArray<cuFloatComplex>*) corrm.get(), dims_to_xform ); // TODO: fixme (requires new cuNDFFT interface)
  conv_kernel.reset();

  // Get the dominant eigenvector for each correlation matrix.
  auto_ptr< cuNDArray<real_complex<REAL> > > csm = extract_csm<REAL>( corrm.get(), ncoils, pixels_per_coil );
  corrm.reset();
  
  // Set phase according to reference (coil 0)
  set_phase_reference<REAL>( csm.get(), ncoils, pixels_per_coil );
  
  return csm;
}

template<class REAL, unsigned int D> __global__ void
set_box_convkernel_kernel( real_complex<REAL> *out, uintd<D> dims, uintd<D> box )
{
  unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < prod(dims) ){
    
    vectord<unsigned int,D> co = idx_to_co(idx,dims);
    vectord<unsigned int,D> offset_dim = (dims>>1);
    vectord<unsigned int,D> offset_box = (box>>1);
    
    if( weak_less(co, offset_dim-offset_box ) || weak_greater_equal(co, offset_dim+offset_box ))
      out[idx] = get_zero<real_complex<REAL> >();
    else{
      real_complex<REAL> _out = get_one<real_complex<REAL> >();
      real_complex<REAL> scale; scale.vec[0] = (REAL)prod(box); scale.vec[1] = get_zero<REAL>();
      scale = reciprocal<real_complex<REAL> >(scale);
      _out *= scale;
      out[idx] = _out;
    }
  }
}

template<class REAL, unsigned int D> auto_ptr< cuNDArray<real_complex<REAL> > >
set_box_convkernel( uintd<D> dims, uintd<D> box )
{
  // Allocate output (aligned)
  cuNDArray<real_complex<REAL> > *out = (sizeof(REAL)==sizeof(float)) ?
    (cuNDArray<real_complex<REAL> >*) cuNDArray<float_complex>::allocate(cuNDA_toVec(dims)) :
    (cuNDArray<real_complex<REAL> >*) cuNDArray<double_complex>::allocate(cuNDA_toVec(dims));
  
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)prod(dims)/blockDim.x));
  
  if( out != 0x0 )
    set_box_convkernel_kernel<REAL,D><<< gridDim, blockDim >>>( out->get_data_ptr(), dims, box );
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr< cuNDArray<real_complex<REAL> > >(out);
}
/*
extern __shared__ char shared_mem[];

// Extract CSM
template<class REAL> __global__ void
extract_csm_kernel( real_complex<REAL> *corrm, real_complex<REAL> *csm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int i = threadIdx.x;

  if( idx < num_elements ){    
    
    // Get the dominant eigenvector for each correlation matrix.
    // Copying Peter Kellman's approach we use the power method:
    //  b_k+1 = A*b_k / ||A*b_k||
    
    real_complex<REAL> *data_out = (real_complex<REAL>*) shared_mem;
    real_complex<REAL> *tmp_v = &(((real_complex<REAL>*) shared_mem)[num_batches*blockDim.x]);

    const unsigned int iterations = 2;

    for( unsigned int c=0; c<num_batches; c++){
      data_out[c*blockDim.x+i] = get_one<real_complex<REAL> >();
    }
    
    for( unsigned int it=0; it<iterations; it++ ){

      for( unsigned int c=0; c<num_batches; c++){
	tmp_v[c*blockDim.x+i] = get_zero<real_complex<REAL> >();
      }
      
      for( unsigned j=0; j<num_batches; j++){
	for( unsigned int k=0; k<num_batches; k++){
	  tmp_v[j*blockDim.x+i] += corrm[(k*num_batches+j)*num_elements+idx]*data_out[k*blockDim.x+i];
	}
      }

      REAL tmp = get_zero<REAL>();
      
      for (unsigned int c=0; c<num_batches; c++){
	tmp += norm_squared(tmp_v[c*blockDim.x+i]);
      }
      
      tmp = sqrt(tmp);
      tmp = reciprocal(tmp);
      
      for (unsigned int c=0; c<num_batches; c++){
	vectord<REAL,2> _res = tmp*tmp_v[c*blockDim.x+i];
	real_complex<REAL> res; res.vec[0]=_res.vec[0]; res.vec[1]=_res.vec[1]; // TODO: do this assignment elegantly
	data_out[c*blockDim.x+i] = res;
      }
    }

    for (unsigned int c=0; c<num_batches; c++){
      csm[c*num_elements+idx] = data_out[c*blockDim.x+i];
    }
  }
}

*/
// Extract CSM
template<class REAL> __global__ void
extract_csm_kernel( real_complex<REAL> *corrm, real_complex<REAL> *csm, unsigned int num_batches, unsigned int num_elements, real_complex<REAL> *tmp_v )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < num_elements ){    
    
    // Get the dominant eigenvector for each correlation matrix.
    // Copying Peter Kellman's approach we use the power method:
    //  b_k+1 = A*b_k / ||A*b_k||
    
    const unsigned int iterations = 2;

    for( unsigned int c=0; c<num_batches; c++){
      csm[c*num_elements+idx] = get_one<real_complex<REAL> >();
    }
    
    for( unsigned int it=0; it<iterations; it++ ){

      for( unsigned int c=0; c<num_batches; c++){
	tmp_v[c*num_elements+idx] = get_zero<real_complex<REAL> >();
      }
      
      for( unsigned j=0; j<num_batches; j++){
	for( unsigned int k=0; k<num_batches; k++){
	  tmp_v[j*num_elements+idx] += corrm[(k*num_batches+j)*num_elements+idx]*csm[k*num_elements+idx];
	}
      }

      REAL tmp = get_zero<REAL>();
      
      for (unsigned int c=0; c<num_batches; c++){
	tmp += norm_squared(tmp_v[c*num_elements+idx]);
      }
      
      tmp = sqrt(tmp);
      tmp = reciprocal(tmp);
      
      for (unsigned int c=0; c<num_batches; c++){
	vectord<REAL,2> _res = tmp*tmp_v[c*num_elements+idx];
	real_complex<REAL> res; res.vec[0]=_res.vec[0]; res.vec[1]=_res.vec[1]; // TODO: do this assignment elegantly
	csm[c*num_elements+idx] = res;
      }
    }
  }
}

// Extract CSM
template<class REAL> __host__ 
auto_ptr< cuNDArray<real_complex<REAL> > > extract_csm( cuNDArray<real_complex<REAL> > *corrm_in, unsigned int number_of_batches, unsigned int number_of_elements )
{
  vector<unsigned int> image_dims;

  for( unsigned int i=0; i<corrm_in->get_number_of_dimensions()-1; i++ ){
    image_dims.push_back(corrm_in->get_size(i));
  }
  
  // Allocate output (aligned)
  cuNDArray<real_complex<REAL> > *out = (sizeof(REAL)==sizeof(float)) ?
    (cuNDArray<real_complex<REAL> >*) cuNDArray<float_complex>::allocate(image_dims) :
    (cuNDArray<real_complex<REAL> >*) cuNDArray<double_complex>::allocate(image_dims);

  dim3 blockDim(128);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  /*  
  if( out != 0x0 )
    extract_csm_kernel<REAL><<< gridDim, blockDim, number_of_batches*blockDim.x*2*sizeof(real_complex<REAL>) >>>
      ( corrm_in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
  */

  // Temporary buffer. TODO: use shared memory
  cuNDArray<real_complex<REAL> > *tmp_v = (sizeof(REAL)==sizeof(float)) ?
    (cuNDArray<real_complex<REAL> >*) cuNDArray<float_complex>::allocate(image_dims) :
    (cuNDArray<real_complex<REAL> >*) cuNDArray<double_complex>::allocate(image_dims);

  if( out != 0x0 )
    extract_csm_kernel<REAL><<< gridDim, blockDim >>>
      ( corrm_in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements, tmp_v->get_data_ptr() );

  CHECK_FOR_CUDA_ERROR();
  
  delete tmp_v;
  return auto_ptr< cuNDArray<real_complex<REAL> > >(out);
}

// Set refence phase
template<class REAL> __global__ void
set_phase_reference_kernel( real_complex<REAL> *csm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < num_elements ){    
    REAL angle = arg(csm[idx]); //Phase of the first coil
    REAL sin_a, cos_a; sin_cos( angle, &sin_a, &cos_a );

    real_complex<REAL> tmp;
    tmp.vec[0] = cos_a; tmp.vec[1] = sin_a;
    tmp = conj(tmp);

    for( unsigned int c=0; c<num_batches; c++ ){
      real_complex<REAL> val =  csm[c*num_elements+idx];
      val *= tmp;
      csm[c*num_elements+idx] = val;
    }
  }
}
  
// Extract CSM
template<class REAL> __host__ 
void set_phase_reference( cuNDArray<real_complex<REAL> > *csm, unsigned int number_of_batches, unsigned int number_of_elements )
{
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  set_phase_reference_kernel<REAL><<< gridDim, blockDim >>>( csm->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

//
// Template instantiation
//

template auto_ptr< cuNDArray<real_complex<float> > > estimate_b1_map<float,2>(cuNDArray<real_complex<float> >*);
template auto_ptr< cuNDArray<real_complex<float> > > estimate_b1_map<float,3>(cuNDArray<real_complex<float> >*);
template auto_ptr< cuNDArray<real_complex<float> > > estimate_b1_map<float,4>(cuNDArray<real_complex<float> >*);
