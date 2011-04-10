#include "b1_map.hcu"
#include "uintd_utilities.hcu"
#include "ndarray_device_utilities.hcu"
#include "check_CUDA.h"
#include "cuNDFFT.h"

#include <math_functions.h>
#include <iostream>

using namespace std;

template<class REAL, class T> __host__ 
auto_ptr< cuNDArray<T> > extract_csm( cuNDArray<T> *corrm_in, unsigned int number_of_batches, unsigned int number_of_elements );

template<class REAL, class T> __host__ 
void set_phase_reference( cuNDArray<T> *csm, unsigned int number_of_batches, unsigned int number_of_elements );

template<class REAL, class T> auto_ptr< cuNDArray<T> >
estimate_b1_map( cuNDArray<T>* data_in )
{
  if( data_in->get_number_of_dimensions() < 2 ){
    cout << endl << "estimate_b1_map:: dimensionality mismatch." << endl; 
    auto_ptr< cuNDArray<T> >(0x0);
  }

  vector<unsigned int> image_dims, dims_to_xform;
  unsigned int pixels_per_coil = 1;
  
  unsigned int d = data_in->get_number_of_dimensions()-1;

  for( unsigned int i=0; i<d; i++ ){
    image_dims.push_back(data_in->get_size(i));
    dims_to_xform.push_back(i);
    pixels_per_coil *= data_in->get_size(i);
  }
  
  unsigned int ncoils = data_in->get_size(d);

  // Make a copy of input data
  cuNDArray<T> data_out = *data_in;

  // Normalize by the RSS of the coils
  if( !cuNDA_rss_normalize<REAL, T>( &data_out, d ) ){
    cout << endl << "estimate_b1_map:: error in rss_normalize" << endl;
    auto_ptr< cuNDArray<T> >(0x0);
  }

  // Now calculate the correlation matrices
  auto_ptr< cuNDArray<T> > corrm = cuNDA_correlation( &data_out );
  
  // Compute smoothing kernel for convolution
  auto_ptr< cuNDArray<T> > conv_kernel; // TODO: compute/set kernel
  
  // Perform convolution by multiplication in image space
  cuNDFFT().fft( conv_kernel.get() );
  cuNDFFT().fft( corrm.get(), dims_to_xform );
  cuNDA_scale( conv_kernel.get(), corrm.get() );
  cuNDFFT().ifft( corrm.get(), dims_to_xform );
  
  // Get the dominant eigenvector for each correlation matrix.
  auto_ptr< cuNDArray<T> > csm = extract_csm<REAL, T>( corrm.get(), ncoils, pixels_per_coil );
  
  // Set phase according to reference (coil 0)
  set_phase_reference<REAL, T>( csm.get(), ncoils, pixels_per_coil );
  
  return auto_ptr< cuNDArray<T> >(csm);
}

extern __shared__ char shared_mem[];

// Extract CSM
template<class REAL, class T> __global__ void
extract_csm_kernel( T *corrm, T *csm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int i = threadIdx.x;

  if( idx < num_elements ){    
    
    // Get the dominant eigenvector for each correlation matrix.
    // Copying Peter Kellman's approach we use the power method:
    //  b_k+1 = A*b_k / ||A*b_k||
    
    T *data_out = (T*) shared_mem;
    T *tmp_v = &(((T*) shared_mem)[num_batches*blockDim.x]);

    const unsigned int iterations = 2;

    for( unsigned int c=0; c<num_batches; c++){
      data_out[c*blockDim.x+i] = get_one<T>();
    }
    
    for( unsigned int it=0; it<iterations; it++ ){

      for( unsigned int c=0; c<num_batches; c++){
	tmp_v[c*blockDim.x+i] = get_zero<T>();
      }
      
      for( unsigned j=0; j<num_batches; j++){
	for( unsigned int k=0; k<num_batches; k++){
	  tmp_v[j*blockDim.x+i] += corrm[(k*num_batches+j)*num_elements+idx]*data_out[k*blockDim.x+i];
	}
      }

      REAL tmp = get_zero<REAL>();
      
      for (unsigned int c=0; c<num_batches; c++){
	tmp += norm_sq(tmp_v[c*blockDim.x+i]);
      }
      
      tmp = sqrt(tmp);
      tmp = reciprocal(tmp);
      
      for (unsigned int c=0; c<num_batches; c++){
	data_out[c*blockDim.x+i] = tmp*tmp_v[c*blockDim.x+i];
      }
    }
      
    for (unsigned int c=0; c<num_batches; c++){
      csm[c*num_elements+idx] = data_out[c*blockDim.x+i];
    }
  }
}

// Extract CSM
template<class REAL, class T> __host__ 
auto_ptr< cuNDArray<T> > extract_csm( cuNDArray<T> *corrm_in, unsigned int number_of_batches, unsigned int number_of_elements )
{
  unsigned int d = corrm_in->get_number_of_dimensions()-2;
  vector<unsigned int> image_dims;

  for( unsigned int i=0; i<d; i++ ){
    image_dims.push_back(corrm_in->get_size(i));
  }
  
  cuNDArray<T> *out = cuNDArray<T>::allocate(image_dims);
  
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  if( out != 0x0 )
    extract_csm_kernel<REAL, T><<< gridDim, blockDim, number_of_batches*blockDim.x*2*sizeof(T) >>>( corrm_in->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
  
  return auto_ptr< cuNDArray<T> >(out);
}

// Set refence phase
template<class REAL, class T> __global__ void
set_phase_reference_kernel( T *csm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;

  if( idx < num_elements ){    
    // TODO: find arg in cuComplex.h eller...
    REAL angle;// = arg(csm[idx]); //Phase of the first coil
    REAL cos_a, sin_a; sincosf( angle, &sin_a, &cos_a ); // TODO overload!

    T tmp = conj(make_reald(cos_a,sin_a));

    for( unsigned int c=0; c<num_batches; c++ ){
      T val =  csm[c*num_elements+idx];
      val *= tmp;
      csm[c*num_elements+idx] = val;
    }
  }
}
  
// Extract CSM
template<class REAL, class T> __host__ 
void set_phase_reference( cuNDArray<T> *csm, unsigned int number_of_batches, unsigned int number_of_elements )
{
  dim3 blockDim(512);
  dim3 gridDim((unsigned int) ceil((double)number_of_elements/blockDim.x));
  
  set_phase_reference_kernel<REAL, T><<< gridDim, blockDim >>>( csm->get_data_ptr(), number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
}

//
// Template instantiation
//

template auto_ptr< cuNDArray<cuFloatComplex> > estimate_b1_map<float, cuFloatComplex>(cuNDArray<cuFloatComplex>*);
