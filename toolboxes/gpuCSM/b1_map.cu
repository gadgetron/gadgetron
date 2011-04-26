#include "b1_map.hcu"
#include "vector_td_operators.hcu"
#include "vector_td_utilities.hcu"
#include "ndarray_vector_td_utilities.hcu"
#include "check_CUDA.h"
#include "cuNDFFT.h"

#include <math_functions.h>
#include <iostream>

using namespace std;

template< class REAL, unsigned int D> __host__ 
auto_ptr< cuNDArray<typename complext<REAL>::Type> > set_box_convkernel( typename uintd<D>::Type dims, typename uintd<D>::Type box );

template<class REAL> __host__ 
auto_ptr< cuNDArray<typename complext<REAL>::Type> > extract_csm( cuNDArray<typename complext<REAL>::Type> *corrm_in, unsigned int number_of_batches, unsigned int number_of_elements );

template<class REAL> __host__ 
void set_phase_reference( cuNDArray<typename complext<REAL>::Type> *csm, unsigned int number_of_batches, unsigned int number_of_elements );

//
// Main method:
//

template<class REAL, unsigned int D> auto_ptr< cuNDArray<typename complext<REAL>::Type> >
estimate_b1_map( cuNDArray<typename complext<REAL>::Type> *data_in )
{
  if( data_in->get_number_of_dimensions() < 2 ){
    cout << endl << "estimate_b1_map:: dimensionality mismatch." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type > >(0x0);
  }

  if( data_in->get_number_of_dimensions()-1 != D ){
    cout << endl << "estimate_b1_map:: dimensionality mismatch." << endl; 
    return auto_ptr< cuNDArray<typename complext<REAL>::Type > >(0x0);
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
  cuNDArray<typename complext<REAL>::Type > *_data_out = new cuNDArray<typename complext<REAL>::Type>(*data_in);
  auto_ptr< cuNDArray<typename complext<REAL>::Type> > data_out(_data_out);
  
  // Normalize by the RSS of the coils
  if( !cuNDA_rss_normalize<REAL, typename complext<REAL>::Type>( data_out.get(), D ) ){
    cout << endl << "estimate_b1_map:: error in rss_normalize" << endl;
    return auto_ptr< cuNDArray<typename complext<REAL>::Type> >(0x0);
  }
  
  // Now calculate the correlation matrices
  auto_ptr<cuNDArray<typename complext<REAL>::Type> > corrm = cuNDA_correlation( data_out.get() );
  data_out.reset();
  
  // Compute smoothing kernel for convolution
  typename uintd<D>::Type dims; cuNDA_fromVec<D>( image_dims, dims );

  typename uintd<D>::Type box; to_vector_td<unsigned int,D>(box,7);
  auto_ptr<cuNDArray<typename complext<REAL>::Type> > conv_kernel = set_box_convkernel<REAL,D>( dims, box );

  // Perform convolution by multiplication in image space
  cuNDFFT().fft( (cuNDArray<cuFloatComplex>*) conv_kernel.get() );           // TODO: fixme (requires new cuNDFFT interface)
  cuNDFFT().fft( (cuNDArray<cuFloatComplex>*) corrm.get(), dims_to_xform );  // TODO: fixme (requires new cuNDFFT interface)
  cuNDA_scale( conv_kernel.get(), corrm.get() );
  cuNDFFT().ifft( (cuNDArray<cuFloatComplex>*) corrm.get(), dims_to_xform ); // TODO: fixme (requires new cuNDFFT interface)
  conv_kernel.reset();

  // Get the dominant eigenvector for each correlation matrix.
  auto_ptr<cuNDArray<typename complext<REAL>::Type> > csm = extract_csm<REAL>( corrm.get(), ncoils, pixels_per_coil );
  corrm.reset();
  
  // Set phase according to reference (coil 0)
  set_phase_reference<REAL>( csm.get(), ncoils, pixels_per_coil );
  
  return csm;
}

template<class REAL, unsigned int D> __global__ void
set_box_convkernel_kernel( typename complext<REAL>::Type *out, typename uintd<D>::Type dims, typename uintd<D>::Type box )
{
  unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < prod(dims) ){
    
    typename uintd<D>::Type co = idx_to_co<D>(idx,dims);
    typename uintd<D>::Type offset_dim = (dims>>1);
    typename uintd<D>::Type offset_box = (box>>1);
    
    if( weak_less(co, offset_dim-offset_box ) || weak_greater_equal(co, offset_dim+offset_box ))
      out[idx] = get_zero<typename complext<REAL>::Type >();
    else{
      typename complext<REAL>::Type _out = get_one<typename complext<REAL>::Type >();
      typename complext<REAL>::Type scale; scale.vec[0] = (REAL)prod(box); scale.vec[1] = get_zero<REAL>();
      scale = reciprocal<typename complext<REAL>::Type >(scale);
      _out *= scale;
      out[idx] = _out;
    }
  }
}

template<class REAL, unsigned int D> auto_ptr<cuNDArray<typename complext<REAL>::Type> >
set_box_convkernel( typename uintd<D>::Type dims, typename uintd<D>::Type box )
{
  // Allocate output (aligned)
  cuNDArray<typename complext<REAL>::Type> *out = cuNDArray<typename complext<REAL>::Type>::allocate(cuNDA_toVec<D>(dims));
  
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)prod(dims)/blockDim.x));
  
  if( out != 0x0 )
    set_box_convkernel_kernel<REAL,D><<< gridDim, blockDim >>>( out->get_data_ptr(), dims, box );
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr<cuNDArray<typename complext<REAL>::Type> >(out);
}
/*
extern __shared__ char shared_mem[];

// Extract CSM
template<class REAL> __global__ void
extract_csm_kernel( typename complext<REAL>::Type *corrm, typename complext<REAL>::Type *csm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int i = threadIdx.x;

  if( idx < num_elements ){    
    
    // Get the dominant eigenvector for each correlation matrix.
    // Copying Peter Kellman's approach we use the power method:
    //  b_k+1 = A*b_k / ||A*b_k||
    
    typename complext<REAL>::Type *data_out = (typename complext<REAL>::Type*) shared_mem;
    typename complext<REAL>::Type *tmp_v = &(((typename complext<REAL>::Type*) shared_mem)[num_batches*blockDim.x]);

    const unsigned int iterations = 2;

    for( unsigned int c=0; c<num_batches; c++){
      data_out[c*blockDim.x+i] = get_one<typename complext<REAL>::Type >();
    }
    
    for( unsigned int it=0; it<iterations; it++ ){

      for( unsigned int c=0; c<num_batches; c++){
	tmp_v[c*blockDim.x+i] = get_zero<typename complext<REAL>::Type >();
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
	typename complext<REAL>::Type res; res.vec[0]=_res.vec[0]; res.vec[1]=_res.vec[1]; // TODO: do this assignment elegantly
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
extract_csm_kernel( typename complext<REAL>::Type *corrm, typename complext<REAL>::Type *csm, unsigned int num_batches, unsigned int num_elements, typename complext<REAL>::Type *tmp_v )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < num_elements ){    
    
    // Get the dominant eigenvector for each correlation matrix.
    // Copying Peter Kellman's approach we use the power method:
    //  b_k+1 = A*b_k / ||A*b_k||
    
    const unsigned int iterations = 2;

    for( unsigned int c=0; c<num_batches; c++){
      csm[c*num_elements+idx] = get_one<typename complext<REAL>::Type >();
    }
    
    for( unsigned int it=0; it<iterations; it++ ){

      for( unsigned int c=0; c<num_batches; c++){
	tmp_v[c*num_elements+idx] = get_zero<typename complext<REAL>::Type >();
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
	typename complext<REAL>::Type res = tmp*tmp_v[c*num_elements+idx];
	csm[c*num_elements+idx] = res;
      }
    }
  }
}

// Extract CSM
template<class REAL> __host__ 
auto_ptr<cuNDArray<typename complext<REAL>::Type> > extract_csm(cuNDArray<typename complext<REAL>::Type> *corrm_in, unsigned int number_of_batches, unsigned int number_of_elements )
{
  vector<unsigned int> image_dims;

  for( unsigned int i=0; i<corrm_in->get_number_of_dimensions()-1; i++ ){
    image_dims.push_back(corrm_in->get_size(i));
  }
  
  // Allocate output (aligned)
  cuNDArray<typename complext<REAL>::Type> *out = cuNDArray<typename complext<REAL>::Type>::allocate(image_dims);

  dim3 blockDim(128);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));

  /*  
  if( out != 0x0 )
    extract_csm_kernel<REAL><<< gridDim, blockDim, number_of_batches*blockDim.x*2*sizeof(typename complext<REAL>::Type) >>>
      ( corrm_in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
  */

  // Temporary buffer. TODO: use shared memory
  cuNDArray<typename complext<REAL>::Type> *tmp_v = cuNDArray<typename complext<REAL>::Type>::allocate(image_dims);

  if( out != 0x0 )
    extract_csm_kernel<REAL><<< gridDim, blockDim >>>
      ( corrm_in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements, tmp_v->get_data_ptr() );

  CHECK_FOR_CUDA_ERROR();
  
  delete tmp_v;
  return auto_ptr<cuNDArray<typename complext<REAL>::Type> >(out);
}

// Set refence phase
template<class REAL> __global__ void
set_phase_reference_kernel( typename complext<REAL>::Type *csm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < num_elements ){    
    REAL angle = arg<REAL>(csm[idx]); //Phase of the first coil
    REAL sin_a, cos_a; sin_cos( angle, &sin_a, &cos_a );

    typename complext<REAL>::Type tmp;
    tmp.vec[0] = cos_a; tmp.vec[1] = sin_a;
    tmp = conj(tmp);

    for( unsigned int c=0; c<num_batches; c++ ){
      typename complext<REAL>::Type val =  csm[c*num_elements+idx];
      val *= tmp;
      csm[c*num_elements+idx] = val;
    }
  }
}
  
// Extract CSM
template<class REAL> __host__ 
void set_phase_reference(cuNDArray<typename complext<REAL>::Type> *csm, unsigned int number_of_batches, unsigned int number_of_elements )
{
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  set_phase_reference_kernel<REAL><<< gridDim, blockDim >>>( csm->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

//
// Template instantiation
//

template auto_ptr< cuNDArray<typename complext<float>::Type > > estimate_b1_map<float,1>(cuNDArray<typename complext<float>::Type >*);
template auto_ptr< cuNDArray<typename complext<float>::Type > > estimate_b1_map<float,2>(cuNDArray<typename complext<float>::Type >*);
template auto_ptr< cuNDArray<typename complext<float>::Type > > estimate_b1_map<float,3>(cuNDArray<typename complext<float>::Type >*);
template auto_ptr< cuNDArray<typename complext<float>::Type > > estimate_b1_map<float,4>(cuNDArray<typename complext<float>::Type >*);

template auto_ptr< cuNDArray<typename complext<double>::Type > > estimate_b1_map<double,1>(cuNDArray<typename complext<double>::Type >*);
template auto_ptr< cuNDArray<typename complext<double>::Type > > estimate_b1_map<double,2>(cuNDArray<typename complext<double>::Type >*);
template auto_ptr< cuNDArray<typename complext<double>::Type > > estimate_b1_map<double,3>(cuNDArray<typename complext<double>::Type >*);
template auto_ptr< cuNDArray<typename complext<double>::Type > > estimate_b1_map<double,4>(cuNDArray<typename complext<double>::Type >*);
