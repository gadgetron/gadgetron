#include "ndarray_vector_td_utilities.h"
#include "real_utilities.h"
#include "real_utilities_device.h"
#include "check_CUDA.h"

#include <cublas_v2.h>

#include <vector>
#include <cmath>
#include <sstream>
#include <boost/throw_exception.hpp>
#include "GadgetronCuException.h"
#include "cuGTBLAS.h"
#include "cudaDeviceManager.h"


// Default template arguments seems to require c++-0x, which we can't assume. 
// We use a dummy type instead...
typedef float dummy;

//
// Some internal utilities
//

// Prepare array for processing:
//
// Sets the device context to the denoted compute device and makes device copies if necessary.
// If compute mode is CUNDA_NDARRAY_DEVICE the mandatory array 'in1' determines the compute device
//

using namespace Gadgetron;
template< unsigned int D, typename I1, typename I2, typename I3 > 
static void prepare( int compute_device, int *cur_device, int *old_device,
		     cuNDArray<I1> *in1,       cuNDArray<I1> **in1_int,
		     cuNDArray<I2> *in2 = 0x0, cuNDArray<I2> **in2_int = 0x0,
		     cuNDArray<I3> *in3 = 0x0, cuNDArray<I3> **in3_int = 0x0 )
{
  // Test validity of D
  if( D==0 || D>3 ){
    BOOST_THROW_EXCEPTION(runtime_error( ">>>Internal error<<< :prepare: D out of range"));

  }

  if( !cur_device || !old_device ){
    BOOST_THROW_EXCEPTION(runtime_error( ">>>Internal error<<< :prepare: device ids 0x0"));

  }

  // Test validity of input pointer
  if( !in1 || !in1_int ){
    BOOST_THROW_EXCEPTION(runtime_error( "unable to process 0x0 input"));

  }
  if( D>1 && (!in2 || !in2_int) ){
    BOOST_THROW_EXCEPTION(runtime_error( "unable to process 0x0 input"));

  }
  if( D>2 && (!in3 || !in3_int) ){
    BOOST_THROW_EXCEPTION(runtime_error( "unable to process 0x0 input"));

  }
  
  // Get current Cuda device
  if( cudaGetDevice(old_device) != cudaSuccess ) {
    BOOST_THROW_EXCEPTION(runtime_error( "unable to get device no"));

  }

  // Set the cuda device to use for computation
  if( compute_device == CUNDA_CURRENT_DEVICE ){
    *cur_device = *old_device; 
  }
  else if( compute_device == CUNDA_NDARRAY_DEVICE ){
    // Let D indicate which ndarray that determines the device
    // D denotes the output array (the latter ndarray in the list), if any, otherwise a sole input ndarray
    if( D == 1)
      *cur_device = in1->get_device();
    else if( D == 2 )
      *cur_device = in2->get_device();
    else if( D == 3 )
      *cur_device = in3->get_device();
  }
  else{
    BOOST_THROW_EXCEPTION(runtime_error( ">>>Internal error<<< :prepare: unknown compute mode"));

  }

  if( *cur_device != *old_device && cudaSetDevice(*cur_device) != cudaSuccess) {
    BOOST_THROW_EXCEPTION(runtime_error( "unable to set device no"));

  }
  
  // Transfer arrays to compute device if necessary
  if( *cur_device != in1->get_device() )
    *in1_int = new cuNDArray<I1>(*in1); // device transfer
  else
    *in1_int = in1;
  
  if( D>1 ){
    if( *cur_device != in2->get_device() )
      *in2_int = new cuNDArray<I2>(*in2); // device transfer
    else
      *in2_int = in2;
  }
  
  if( D>2 ){
    if( *cur_device != in3->get_device() )
      *in3_int = new cuNDArray<I3>(*in3); // device transfer
    else
      *in3_int = in3;
  }
  

}  

// Restore active device and free internal memory 
//
template< unsigned int D, typename I1, typename O, typename I2, typename I3 > 
static void restore( int old_device,
		     cuNDArray<I1> *in1, cuNDArray<I1> *in1_int,
		     unsigned int out_idx = 0, cuNDA_device alloc_device = CUNDA_NDARRAY_DEVICE, cuNDArray<O>  *out = 0x0,
		     cuNDArray<I2> *in2 = 0x0, cuNDArray<I2> *in2_int = 0x0,
		     cuNDArray<I3> *in3 = 0x0, cuNDArray<I3> *in3_int = 0x0 )
{
  // Test validity of D
  if( D==0 || D>3 ){
    BOOST_THROW_EXCEPTION(runtime_error( ">>>Internal error<<< :prepare: D out of range"));

  }

  // Test validity of input pointer
  if( !in1 || !in1_int ){
    BOOST_THROW_EXCEPTION(runtime_error( "unable to process 0x0 input"));

  }
  if( D>1 && (!in2 || !in2_int) ){
    BOOST_THROW_EXCEPTION(runtime_error( "unable to process 0x0 input"));

  }
  if( D>2 && (!in3 || !in3_int) ){
    BOOST_THROW_EXCEPTION(runtime_error( "unable to process 0x0 input"));

  }

  // Check if output ndarray resides on the desired device
  //
  if( out ){
    if( alloc_device == CUNDA_CURRENT_DEVICE && out->get_device() != old_device ){
      out->set_device( old_device ); } // device copy
    else if( alloc_device == CUNDA_NDARRAY_DEVICE && out->get_device() != in1->get_device() ){
      out->set_device( in1->get_device() ); } // device copy
  }

  // Check if in_out ndarray resides on the desired device
  //
  if( out_idx > 0 && out_idx < 4 ){ 

   if( out_idx > D ){
      BOOST_THROW_EXCEPTION(runtime_error( ">>>Internal error<<< :restore: array index out of range"));

    }
   
   if( D == 1 ){
     if( in1->get_device() != in1_int->get_device() ){ 
       *in1 = *in1_int; } // device transfer by assignment
   }
   if( D == 2 ){
     if( out_idx == 1 && in1->get_device() != in1_int->get_device() ){ 
       *in1 = *in1_int; } // device transfer by assignment
     else if( out_idx == 2 && in2->get_device() != in2_int->get_device() ){ 
       *in2 = *in2_int; } // device transfer by assignment
   }
   if( D == 3 ){
     if( out_idx == 1 && in1->get_device() != in1_int->get_device() ){ 
       *in1 = *in1_int; } // device transfer by assignment
     else if( out_idx == 2 && in2->get_device() != in2_int->get_device() ){ 
       *in2 = *in2_int; } // device transfer by assignment
     else if( out_idx == 3 && in3->get_device() != in3_int->get_device() ){ 
       *in3 = *in3_int; } // device transfer by assignment
   }
  }
  else if( out_idx != 0 ){
    BOOST_THROW_EXCEPTION(runtime_error( ">>>Internal error<<< :restore: illegal device specified"));

  }

  // Check if internal array needs deletion (they do only if they were created in ::prepare()
  //
  if( in1->get_device() != in1_int->get_device() ){
    delete in1_int;
  }
  if( D>1 && in2->get_device() != in2_int->get_device() ){
    delete in2_int;
  }
  if( D>2 && in3->get_device() != in3_int->get_device() ){
    delete in3_int;
  }

  // Get current Cuda device
  int device;
  if( cudaGetDevice(&device) != cudaSuccess ) {
    BOOST_THROW_EXCEPTION(cuda_error( "unable to get device no"));

  }

  // Restore old device
  if( device != old_device && cudaSetDevice(old_device) != cudaSuccess) {
    BOOST_THROW_EXCEPTION(cuda_error( "unable to restore device no"));

  }
    

}

// Common block/grid configuration utility
//
static void setup_grid( unsigned int cur_device, unsigned int number_of_elements,
			dim3 *blockDim, dim3* gridDim, unsigned int num_batches=1 )
{

  // For small arrays we keep the block dimension fairly small
  *blockDim = dim3(256);
  *gridDim = dim3((number_of_elements+blockDim->x-1)/blockDim->x, num_batches);
  int maxGridDim = cudaDeviceManager::Instance()->max_griddim(cur_device);
  // Extend block/grid dimensions for large arrays
  if( gridDim->x > maxGridDim){
    blockDim->x = maxGridDim;
    gridDim->x = (number_of_elements+blockDim->x-1)/blockDim->x;
  }

  if( gridDim->x > maxGridDim ){
    gridDim->x = ((unsigned int)sqrt((float)number_of_elements)+blockDim->x-1)/blockDim->x;
    gridDim->y *= ((number_of_elements+blockDim->x*gridDim->x-1)/(blockDim->x*gridDim->x));
  }
   
  if( gridDim->x >maxGridDim || gridDim->y >maxGridDim){

    BOOST_THROW_EXCEPTION(cuda_error("Grid dimension larger than supported by device"));
  }


}

// Common stride setup utility
//
template<class T> static void find_stride( cuNDArray<T> *in, unsigned int dim,
					   unsigned int *stride, std::vector<unsigned int> *dims )
{
  *stride = 1;
  for( unsigned int i=0; i<in->get_number_of_dimensions(); i++ ){
    if( i != dim )
      dims->push_back(in->get_size(i));
    if( i < dim )
      *stride *= in->get_size(i);
  }
}

//
// Implementation of public utilities
//

// cAbs
//
template<class REAL, class T> __global__ void
cAbs_kernel( T *in, REAL *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx<number_of_elements ){
    T val = in[idx];
    out[idx] = Gadgetron::abs(val);
  }
}

// Abs
//
template<class T>
boost::shared_ptr< cuNDArray<typename realType<T>::type > >
Gadgetron::abs( cuNDArray<T> *in,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
	typedef typename realType<T>::type REAL;
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Prepare 
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) cAbs_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}


template<typename T>
class minmax_clamp_functor : public thrust::unary_function<T,T>
{
public:
	minmax_clamp_functor(T _min,T _max):min(_min),max(_max) {};
	 __inline__ __host__ __device__ T operator()(const T &y) const {
		 if (y < min) return min;
		 else if (y > max) return max;
		 else return y;
	 }
private:
	const T min,max;
};


template<typename T>
class minmax_clamp_functor<complext<T> > : public thrust::unary_function<complext<T>, complext<T> >
{
public:
	minmax_clamp_functor(T _min,T _max):min(_min),max(_max) {};
	 __inline__ __host__ __device__ complext<T>  operator()(const complext<T>  &y) const {
		 if (real(y) < min) return  complext<T>(min);
		 else if (real(y) > max) return  complext<T>(max);
		 else return  complext<T>(real(y));
	 }
private:
	const T min,max;
};

// CLAMP functions
template<class  T> EXPORTGPUCORE
void Gadgetron::clamp(cuNDArray<T> *in_out, typename realType<T>::type min, typename realType<T>::type max){
	thrust::transform(in_out->begin(),in_out->end(),in_out->begin(),minmax_clamp_functor<T>(min,max));
}

template<typename T>
class max_clamp_functor : public thrust::unary_function<T,T>
{
public:

	max_clamp_functor(T _max):max(_max) {};
	 __inline__ __host__ __device__ T operator()(const T &y) const {
		 if (y > max) return max;
		 return y;
	 }
private:
	const T max;
};

template<typename T>
class max_clamp_functor<complext<T> > : public thrust::unary_function<complext<T>, complext<T> >
{
public:
		max_clamp_functor(T _max):max(_max) {};
	 __inline__ __host__ __device__ complext<T>  operator()(const complext<T>  &y) const {
		 if (real(y) > max) return  complext<T>(max);
		 else return  complext<T>(real(y));
	 }
private:
	const T max;
};
// CLAMP functions
template<class  T> EXPORTGPUCORE
void Gadgetron::clamp_max(cuNDArray<T> *in_out, typename realType<T>::type max){
	thrust::transform(in_out->begin(),in_out->end(),in_out->begin(),max_clamp_functor<T>(max));
}

template<typename T>
class min_clamp_functor : public thrust::unary_function<T,T>
{
public:

	min_clamp_functor(T _min):min(_min) {};
	 __inline__ __host__ __device__ T operator()(const T &y) const {
		 if (y < min) return min;
		 return y;
	 }
private:
	const T min;
};

template<typename T>
class min_clamp_functor<complext<T> > : public thrust::unary_function<complext<T>, complext<T> >
{
public:
	min_clamp_functor(T _min):min(_min) {};
	 __inline__ __host__ __device__ complext<T>  operator()(const complext<T>  &y) const {
		 if (real(y) < min) return  complext<T>(min);
		 else return  complext<T>(real(y));
	 }
private:
	const T min;
};
// CLAMP functions
template<class  T> EXPORTGPUCORE
void Gadgetron::clamp_min(cuNDArray<T> *in_out, typename realType<T>::type min){
	thrust::transform(in_out->begin(),in_out->end(),in_out->begin(),min_clamp_functor<T>(min));
}

// Sum
//
template<class T> __global__ void
sum_kernel( T *in, T *out,
		  unsigned int stride, 
		  unsigned int number_of_batches, 
		  unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){

    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
 
    T val = in[in_idx];
 
    for( unsigned int i=1; i<number_of_batches; i++ ) 
      val += in[i*stride+in_idx];

    out[idx] = val; 
  }
}

// Sum
//
template<class T>  
boost::shared_ptr< cuNDArray<T> >
Gadgetron::sum( cuNDArray<T> *in, unsigned int dim,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
  
  // Some validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    BOOST_THROW_EXCEPTION(runtime_error("sum: underdimensioned."));
  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    BOOST_THROW_EXCEPTION(runtime_error( "sum: dimension out of range."));
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );
 
  // Find element stride
  unsigned int stride; std::vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(&dims);
  if( out.get() != 0x0 ) sum_kernel<T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,T,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// Expand
//
template<class T> __global__ void
expand_kernel( T *in, T *out,
		     unsigned int number_of_elements,
		     unsigned int new_dim_size )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = in[idx%(number_of_elements/new_dim_size)];
  }
}

// Expand
//
template<class T>  
boost::shared_ptr< cuNDArray<T> >
Gadgetron::expand( cuNDArray<T> *in, unsigned int new_dim_size,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
   
  unsigned int number_of_elements = in->get_number_of_elements()*new_dim_size;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );
 
  // Find element stride
  std::vector<unsigned int> dims = *in->get_dimensions();
  dims.push_back(new_dim_size);

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(&dims);
  if( out.get() != 0x0 ) expand_kernel<T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), number_of_elements, new_dim_size );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,T,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// SS
template<class REAL, class T> __inline__  __device__ REAL
_ss( unsigned int idx, T *in, unsigned int stride, unsigned int number_of_batches )
{
  unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
  REAL ss = REAL(0);

  for( unsigned int i=0; i<number_of_batches; i++ )
    ss += norm(in[i*stride+in_idx]);

  return ss;
}

// SS
template<class REAL, class T> __global__ void
ss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = _ss<REAL,T>(idx, in, stride, number_of_batches);
  }
}

// squaredNorm
template<class T>
boost::shared_ptr< cuNDArray<typename realType<T>::type> >
Gadgetron::squaredNorm( cuNDArray<T> *in, unsigned int dim,
	   cuNDA_device alloc_device, cuNDA_device compute_device )
{
	typedef typename realType<T>::type REAL;
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    BOOST_THROW_EXCEPTION(runtime_error( "ss: underdimensioned."));

  }

  if( dim > in->get_number_of_dimensions()-1 ){
    BOOST_THROW_EXCEPTION(runtime_error( "ss: dimension out of range."));

  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );

  // Find element stride
  unsigned int stride; std::vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if ( out.get() != 0x0 ) ss_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );

  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// RSS
template<class REAL, class T> __inline__  __device__ REAL
_rss( unsigned int idx, T *in, unsigned int stride, unsigned int number_of_batches )
{
  unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);
  REAL rss = REAL(0);
  
  for( unsigned int i=0; i<number_of_batches; i++ ) 
    rss += norm(in[i*stride+in_idx]);
  
  rss = sqrt(rss); 

  return rss;
}

// RSS
template<class REAL, class T> __global__ void
rss_kernel( T *in, REAL *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){
    out[idx] = _rss<REAL,T>(idx, in, stride, number_of_batches); 
  }
}
/*
// RSS
template<class REAL, class T>  
boost::shared_ptr< cuNDArray<REAL> >
_rss( cuNDArray<T> *in, unsigned int dim,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
    BOOST_THROW_EXCEPTION(runtime_error( "rss: underdimensioned."));

  }
 
  if( dim > in->get_number_of_dimensions()-1 ){
    BOOST_THROW_EXCEPTION(runtime_error( "rss: dimension out of range."));
  }

  unsigned int number_of_batches = in->get_size(dim);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );

  // Find element stride
  unsigned int stride; std::vector<unsigned int> dims;
  find_stride<T>( in, dim, &stride, &dims );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims); 
  if ( out.get() != 0x0 ) rss_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), stride, number_of_batches, number_of_elements );
  
  CHECK_FOR_CUDA_ERROR();
 
  // Restore
  restore<1,T,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}*/

// cRSS
template<class REAL, class T> __global__ void
crss_kernel( T *in, T *out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){

    REAL rss = _rss<REAL,T>(idx, in, stride, number_of_batches); 

    out[idx].vec[0] = rss;
    out[idx].vec[1] = REAL(0);
  }
}


// Build correlation matrix
template<class REAL, class T> __global__ void
correlation_kernel( T *in, T *corrm, unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int p = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int i = threadIdx.y;

  if( p < num_elements ){
    for( unsigned int j=0; j<i; j++){
      T tmp = in[i*num_elements+p]*conj(in[j*num_elements+p]);
      corrm[(j*num_batches+i)*num_elements+p] = tmp;
      corrm[(i*num_batches+j)*num_elements+p] = conj(tmp);
    }
    T tmp = in[i*num_elements+p];
    corrm[(i*num_batches+i)*num_elements+p] = tmp*conj(tmp);
  }
}

// Build correlation matrix
template<class T>
boost::shared_ptr< cuNDArray<T> >
Gadgetron::correlation( cuNDArray<T> *in,
		    cuNDA_device alloc_device, cuNDA_device compute_device )
{
	typedef typename realType<T>::type REAL;
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );

  // Validity checks
  if( !(in->get_number_of_dimensions()>1) ){
  	BOOST_THROW_EXCEPTION(runtime_error("correlation: underdimensioned."));
  }
 
  unsigned int number_of_batches = in->get_size(in->get_number_of_dimensions()-1);
  unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;

  int warp_size = cudaDeviceManager::Instance()->warp_size(old_device);
  int max_blockdim = cudaDeviceManager::Instance()->max_blockdim(old_device);
  dim3 blockDim(((max_blockdim/number_of_batches)/warp_size)*warp_size, number_of_batches);

  if( blockDim.x == 0 ){
  	BOOST_THROW_EXCEPTION(runtime_error("correlation: correlation dimension exceeds device capacity."));
  }
  
  dim3 gridDim((number_of_elements+blockDim.x-1)/blockDim.x);

  // Invoke kernel
  std::vector<unsigned int> dims = *in->get_dimensions(); dims.push_back(number_of_batches);
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(&dims);
  if( out.get() != 0x0 ) correlation_kernel<REAL,T><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,T,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// Real to complext
template<class REAL> __global__ void
real_to_complext_kernel( REAL *in, complext<REAL> *out, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < num_elements ){
  	complext<REAL> z;
    z.vec[0] = in[idx];
    z.vec[1] = REAL(0);
    out[idx] = z;
  }
}

// Convert real to complext
template<class REAL>  
boost::shared_ptr< cuNDArray<complext<REAL> > >
Gadgetron::real_to_complext( cuNDArray<REAL> *in,
			cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
 
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<complext<REAL> > > out = cuNDArray<complext<REAL> >::allocate(in->get_dimensions().get());
  if( out.get() != 0x0 ) real_to_complext_kernel<REAL><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements());
  
  CHECK_FOR_CUDA_ERROR();
  
  // Restore
  restore<1,REAL,complext<REAL>,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// complext to real by cropping the imaginary component 
template<class REAL> __global__ void
complext_to_real_kernel( complext<REAL> *in, REAL *out, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < num_elements ){
    out[idx] = in[idx].vec[0];
  }
}

// Convert complext to real by cropping the imaginary component 
template<class REAL>  
boost::shared_ptr< cuNDArray<REAL> >
Gadgetron::complext_to_real( cuNDArray<complext<REAL> > *in,
			cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<complext<REAL> > *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,complext<REAL>,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int ) ;
 
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(in->get_dimensions().get());  
  if( out.get() != 0x0 ) complext_to_real_kernel<REAL><<< gridDim, blockDim >>>( in_int->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements());
  
  CHECK_FOR_CUDA_ERROR();
  
  // Restore
  restore<1,complext<REAL>,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// Downsample
template<class REAL, unsigned int D> __global__ void
downsample_kernel( REAL *in, REAL *out,
			 vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
			 unsigned int num_elements, unsigned int num_batches )
{
	typedef vector_td<unsigned int,D> uintd;
  // We have started a thread for each output element
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
      if( weak_greater_equal( local_co, matrix_size_out ) ) continue; // To allow array dimensions of 1
      const unsigned int in_idx = co_to_idx<D>(co_in+local_co, matrix_size_in)+frame_offset*prod(matrix_size_in);
      actual_adds++;
      res += in[in_idx];
    }
    
    out[idx] = res/REAL(actual_adds);
  }
}

// Downsample
template<class REAL, unsigned int D>
boost::shared_ptr< cuNDArray<REAL> >
Gadgetron::downsample( cuNDArray<REAL> *in,
		  cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
     
  // A few sanity checks 
  if( in->get_number_of_dimensions() < D ){
    BOOST_THROW_EXCEPTION(runtime_error( "downsample: the number of array dimensions should be at least D"));

  }
  
  for( unsigned int d=0; d<D; d++ ){
    if( (in->get_size(d)%2) == 1 && in->get_size(d) != 1 ){
      BOOST_THROW_EXCEPTION(runtime_error( "downsample: uneven array dimensions larger than one not accepted"));
    }
  }
  
  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
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

  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, number_of_batches ) ;
  
  // Invoke kernel
  std::vector<unsigned int> dims = uintd_to_vector<D>(matrix_size_out);
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    dims.push_back(in->get_size(d));
  }
  
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if( out.get() != 0x0 ) 
    downsample_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in_int->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,REAL,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

// Nearest neighbor upsampling
template<class REAL, unsigned int D> __global__ void
upsample_nn_kernel( REAL *in, REAL *out,
		vector_td<unsigned int,D> matrix_size_in, vector_td<unsigned int,D> matrix_size_out,
		       unsigned int num_elements, unsigned int num_batches )
{
	typedef vector_td<unsigned int,D> uintd;
  // We have started a thread for each output element
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  
  if( idx < num_elements*num_batches ){    
    const unsigned int frame_idx = idx/num_elements;
    const uintd co_out = idx_to_co<D>( idx-frame_idx*num_elements, matrix_size_out );
    const uintd co_in = co_out >> 1;
    out[idx] = in[co_to_idx<D>(co_in, matrix_size_in)+frame_idx*prod(matrix_size_in)];
  }
}

// Nearest neighbor upsampling
template<class REAL, unsigned int D>
boost::shared_ptr< cuNDArray<REAL> >
Gadgetron::upsample_nn( cuNDArray<REAL> *in,
		   cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
     
  // A few sanity checks 
  if( in->get_number_of_dimensions() < D ){
    BOOST_THROW_EXCEPTION(runtime_error( "upsample: the number of array dimensions should be at least D" ));

  }
    
  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = matrix_size_in << 1;

  unsigned int number_of_elements = prod(matrix_size_out);
  unsigned int number_of_batches = 1;

  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;

  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, number_of_batches );
  
  // Invoke kernel
  std::vector<unsigned int> dims = uintd_to_vector<D>(matrix_size_out);
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    dims.push_back(in->get_size(d));
  }
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if( out.get() != 0x0 ) 
    upsample_nn_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in_int->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,REAL,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

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
	typedef vector_td<unsigned int,D> uintd;
  // We have started a thread for each output element
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
template<class REAL, unsigned int D>
boost::shared_ptr< cuNDArray<REAL> >
Gadgetron::upsample_lin( cuNDArray<REAL> *in,
		    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<REAL> *in_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,REAL,dummy,dummy>( compute_device, &cur_device, &old_device, in, &in_int );
     
  // A few sanity checks 
  if( in->get_number_of_dimensions() < D ){
    BOOST_THROW_EXCEPTION(runtime_error( "upsample: the number of array dimensions should be at least D"));
  }
    
  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
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

  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim, number_of_batches );
  
  // Invoke kernel
  std::vector<unsigned int> dims = uintd_to_vector<D>(matrix_size_out);
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    dims.push_back(in->get_size(d));
  }
  boost::shared_ptr< cuNDArray<REAL> > out = cuNDArray<REAL>::allocate(&dims);
  if( out.get() != 0x0 ) 
    upsample_lin_kernel<REAL,D><<< gridDim, blockDim >>>
      ( in_int->get_data_ptr(), out->get_data_ptr(), matrix_size_in, matrix_size_out, number_of_elements, number_of_batches );
  
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,REAL,REAL,dummy,dummy>( old_device, in, in_int, 0, alloc_device, out.get() );

  return out;
}

template<typename T>
struct sign_functor
{
	  __host__ __device__
	  T operator()(const T & x) const{
		  return sgn(x);
	  }
};

template <class T> void Gadgetron::inplace_sgn(cuNDArray<T>* x) {
		  thrust::device_ptr<T> dev_ptr(x->get_data_ptr());
		  thrust::transform(dev_ptr, dev_ptr + x->get_number_of_elements(),dev_ptr,sign_functor<T>());
}


// Reciprocal square root
template<class T> __global__ 
void reciprocal_sqrt_kernel( T *in_out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    in_out[idx] = rsqrt(in_out[idx]);
  }
}


// Normalized RSS
template<class REAL, class T> __global__ void
rss_normalize_kernel( T *in_out, unsigned int stride, unsigned int number_of_batches, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < number_of_elements ){

    REAL reciprocal_rss = 1/(_rss<REAL,T>(idx, in_out, stride, number_of_batches));
 
    unsigned int in_idx = (idx/stride)*stride*number_of_batches+(idx%stride);

    for( unsigned int i=0; i<number_of_batches; i++ ) {
      T out = in_out[i*stride+in_idx];
      out *= reciprocal_rss; // complex-scalar multiplication (element-wise operator)
      in_out[i*stride+in_idx] = out; 
    } 
  }
}

// Normalized RSS
template<class T>
void Gadgetron::rss_normalize( cuNDArray<T> *in_out, unsigned int dim,
			  cuNDA_device compute_device )
{
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int );

  // Validity checks
  if( !(in_out->get_number_of_dimensions()>1) ){
    BOOST_THROW_EXCEPTION(runtime_error( "rss_normalize: underdimensioned."));

  }
 
  if( dim > in_out->get_number_of_dimensions()-1 ){
  	BOOST_THROW_EXCEPTION(runtime_error("rss_normalize: dimension out of range."));

  }

  unsigned int number_of_batches = in_out->get_size(dim);
  unsigned int number_of_elements = in_out->get_number_of_elements()/number_of_batches;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, number_of_elements, &blockDim, &gridDim );
  // Find element stride
  unsigned int stride; std::vector<unsigned int> dims;
  find_stride<T>( in_out, dim, &stride, &dims );

  // Invoke kernel
  rss_normalize_kernel<typename realType<T>::type,T><<< gridDim, blockDim >>>( in_out_int->get_data_ptr(), stride, number_of_batches, number_of_elements );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 );


}





// Normalize (float)
template<class T> EXPORTGPUCORE
T Gadgetron::normalize( cuNDArray<T> *data, T new_max, cuNDA_device compute_device )
{

  unsigned int number_of_elements = data->get_number_of_elements();

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *data_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, data, &data_int );

  // Find the maximum value in the array
  int max_idx=amax(data_int);

  cudaThreadSynchronize();
  
  // Copy that value back to host memory
  T max_val;
  cudaMemcpy(&max_val, (data_int->get_data_ptr()+max_idx-1), sizeof(T), cudaMemcpyDeviceToHost);

  // Scale the array
  T scale = std::abs(new_max/max_val);
  *data_int *= scale;
  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, data, data_int, 1, compute_device );

  CHECK_FOR_CUDA_ERROR();
  return scale;
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
template<class T, unsigned int D> EXPORTGPUCORE
void Gadgetron::crop( typename uintd<D>::Type offset,
	    cuNDArray<T> *in, cuNDArray<T> *out,
	    cuNDA_device compute_device )
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

  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *out->get_dimensions() );
 
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }

  if( weak_greater(offset+matrix_size_out, matrix_size_in) ){
    BOOST_THROW_EXCEPTION(runtime_error( "crop: cropping size mismatch"));

  }

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int, *out_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in, &in_int, out, &out_int );
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );

  // Invoke kernel
  crop_kernel<T,D><<< gridDim, blockDim >>>
    ( offset, matrix_size_in, matrix_size_out, in_int->get_data_ptr(), out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, in, in_int, 2, compute_device, 0x0, out, out_int );

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
void Gadgetron::expand_with_zero_fill( cuNDArray<T> *in, cuNDArray<T> *out,
				  cuNDA_device compute_device )
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

  typename uintd<D>::Type matrix_size_in = vector_to_uintd<D>( *in->get_dimensions() );
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *out->get_dimensions() );
  
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in->get_number_of_dimensions(); d++ ){
    number_of_batches *= in->get_size(d);
  }

  if( weak_greater(matrix_size_in,matrix_size_out) ){
    std::runtime_error("expand: size mismatch, cannot expand");

  }
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int, *out_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in, &in_int, out, &out_int );

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim, number_of_batches );
 
  // Invoke kernel
  expand_with_zero_fill_kernel<T,D><<< gridDim, blockDim >>> ( matrix_size_in, matrix_size_out, in_int->get_data_ptr(), out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, in, in_int, 2, compute_device, 0x0, out, out_int );

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
void Gadgetron::zero_fill_border( typename uintd<D>::Type matrix_size_in, cuNDArray<T> *in_out,
			     cuNDA_device compute_device )
{ 
  typename uintd<D>::Type matrix_size_out = vector_to_uintd<D>( *in_out->get_dimensions() );
 
  if( weak_greater(matrix_size_in, matrix_size_out) ){
    BOOST_THROW_EXCEPTION(runtime_error("zero_fill: size mismatch, cannot zero fill"));

  }
 
  unsigned int number_of_batches = 1;
  for( unsigned int d=D; d<in_out->get_number_of_dimensions(); d++ ){
    number_of_batches *= in_out->get_size(d);
  }

 // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int ) ;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, prod(matrix_size_out), &blockDim, &gridDim );
 
  // Invoke kernel
  zero_fill_border_kernel<T,D><<< gridDim, blockDim >>>
  		( matrix_size_in, matrix_size_out, in_out_int->get_data_ptr(), number_of_batches, prod(matrix_size_out) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 );


}

// Zero fill border (circular)
template<class REAL, class T, unsigned int D> __global__ void
zero_fill_border_kernel( REAL radius, vector_td<int,D> dims, T *image,
			       unsigned int num_batches, unsigned int num_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
  const unsigned int frame_offset = idx/num_elements;
  
  if( idx < num_elements*num_batches ){
    const vector_td<int,D> co = idx_to_co<D>( idx-frame_offset*num_elements, dims ) - (dims>>1);
    if( REAL(norm_squared(co)) > radius*radius )
      image[idx] = T(0);
  }
}

// Zero fill border (circular, 2D)
template<class REAL, class T, unsigned int D> 
void Gadgetron::zero_fill_border( REAL radius, cuNDArray<T> *in_out,
			     cuNDA_device compute_device )
{

  
  unsigned int number_of_batches = 1;
  for( unsigned int d=2; d<in_out->get_number_of_dimensions(); d++ ){
    number_of_batches *= in_out->get_size(d);
  }

  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_out_int;

  // Perform device copy if array is not residing on the current device
  prepare<1,T,dummy,dummy>( compute_device, &cur_device, &old_device, in_out, &in_out_int );

  typename intd<D>::Type dims = vector_to_intd<D>(*in_out->get_dimensions());

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in_out->get_number_of_elements(), &blockDim, &gridDim, number_of_batches );
 
  // Invoke kernel
  zero_fill_border_kernel<REAL,T,D><<< gridDim, blockDim >>>
    ( radius, dims, in_out_int->get_data_ptr(), number_of_batches, prod(dims) );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<1,T,dummy,dummy,dummy>( old_device, in_out, in_out_int, 1 );


}

// Shrinkage
//

template<class REAL, class T> __global__ void 
shrink1_kernel( REAL gamma, T *in, T *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    T in_val = in[idx]; 
    REAL in_norm = abs(in_val);
    T _res;
    if( in_norm > REAL(0) )
      _res =  in_val/in_norm;
    else
      _res = T(0);
    REAL maximum = max( in_norm-gamma, REAL(0) );
    T res = maximum*_res;

    out[idx] = res;
  }
}

template<class REAL, class T> EXPORTGPUCORE
void Gadgetron::shrink1( REAL gamma, cuNDArray<T> *in, cuNDArray<T> *out )
{
  // TODO: multi-device handling

  if( !in || !out ){
    BOOST_THROW_EXCEPTION(runtime_error( "shrink1: 0x0 arrays not accepted" ));

  }

  if( in->get_number_of_elements() != out->get_number_of_elements() ){
    BOOST_THROW_EXCEPTION(runtime_error( "shrink1: i/o arrays must have an identical number of elements"));
  }
  
  // Get current Cuda device
  int cur_device;
  if( cudaGetDevice(&cur_device) != cudaSuccess ) {
    BOOST_THROW_EXCEPTION(runtime_error( "shrink1 : unable to get device no"));
  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );
  
  // Invoke kernel
  shrink1_kernel<REAL,T><<< gridDim, blockDim >>>( gamma, in->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();
  

}

template<class REAL, class T> __global__ void 
shrinkd_kernel( REAL gamma, REAL *s_k, T *in, T *out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    T in_val = in[idx]; 
    REAL s_k_val = s_k[idx];
    T _res;
    if( s_k_val > 0 )
      _res =  in_val/s_k_val;
    else
      _res = T(0);
    REAL maximum = max( s_k_val-gamma, REAL(0) );
    T res = maximum*_res;

    out[idx] = res;
  }
}

template<class REAL, class T> EXPORTGPUCORE
void Gadgetron::shrinkd( REAL gamma, cuNDArray<REAL> *s_k, cuNDArray<T> *in, cuNDArray<T> *out )
{
  // TODO: multi-device handling

  if( !in || !out || !s_k ){
    BOOST_THROW_EXCEPTION(runtime_error( "shrinkd: 0x0 arrays not accepted"));

  }

  if( in->get_number_of_elements() != out->get_number_of_elements() ){
    BOOST_THROW_EXCEPTION(runtime_error( "shrinkd: i/o arrays must have an identical number of elements"));

  }

  if( in->get_number_of_elements() != s_k->get_number_of_elements() ){
    BOOST_THROW_EXCEPTION(runtime_error( "shrinkd: i/o arrays must have an identical number of elements"));

  }
  
  // Get current Cuda device
  int cur_device;
  if( cudaGetDevice(&cur_device) != cudaSuccess ) {
    BOOST_THROW_EXCEPTION(runtime_error( "shrinkd : unable to get device no"));

  }

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in->get_number_of_elements(), &blockDim, &gridDim );
  
  // Invoke kernel
  shrinkd_kernel<REAL,T><<< gridDim, blockDim >>>( gamma, s_k->get_data_ptr(), in->get_data_ptr(), out->get_data_ptr(), in->get_number_of_elements() );
  
  CHECK_FOR_CUDA_ERROR();

}

// Mirror, but keep the origin unchanged
template<class T, unsigned int D> __global__ void
origin_mirror_kernel( vector_td<unsigned int,D> matrix_size, vector_td<unsigned int,D> origin, T *in, T *out, bool zero_fill )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;

  if( idx < prod(matrix_size) ){

  	vector_td<unsigned int,D> in_co = idx_to_co<D>( idx, matrix_size );
  	vector_td<unsigned int,D> out_co = matrix_size-in_co;
    
    bool wrap = false;
    for( unsigned int d=0; d<D; d++ ){
      if( out_co.vec[d] == matrix_size.vec[d] ){
	out_co.vec[d] = 0;
	wrap = true;
      }
    }
    
    const unsigned int in_idx = co_to_idx<D>(in_co, matrix_size);
    const unsigned int out_idx = co_to_idx<D>(out_co, matrix_size);

    if( wrap && zero_fill )
      out[out_idx] = T(0);
    else
      out[out_idx] = in[in_idx];
  }
}

// Mirror around the origin -- !! leaving the origin unchanged !!
// This creates empty space "on the left" that can be filled by zero (default) or the left-over entry.
template<class T, unsigned int D> EXPORTGPUCORE
void Gadgetron::origin_mirror( cuNDArray<T> *in, cuNDArray<T> *out, bool zero_fill, cuNDA_device compute_device )
{
  if( in == 0x0 || out == 0x0 ){
  	BOOST_THROW_EXCEPTION(runtime_error( "origin_mirror: 0x0 ndarray provided"));

  }

  if( !in->dimensions_equal(out) ){
  	BOOST_THROW_EXCEPTION(runtime_error("origin_mirror: image dimensions mismatch"));

  }
  
  if( in->get_number_of_dimensions() != D ){
  	std::stringstream ss;
    ss << "origin_mirror: number of image dimensions is not " << D;
    BOOST_THROW_EXCEPTION(runtime_error(ss.str()));
  }

  typename uintd<D>::Type matrix_size = vector_to_uintd<D>( *in->get_dimensions() );
 
  // Prepare internal array
  int cur_device, old_device;
  cuNDArray<T> *in_int, *out_int;

  // Perform device copy if array is not residing on the current device
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in, &in_int, out, &out_int );
  
  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, prod(matrix_size), &blockDim, &gridDim );

  // Invoke kernel
  origin_mirror_kernel<T,D><<< gridDim, blockDim >>> ( matrix_size, matrix_size>>1, in_int->get_data_ptr(), out_int->get_data_ptr(), zero_fill );
 
  CHECK_FOR_CUDA_ERROR();

  // Restore
  restore<2,T,dummy,T,dummy>( old_device, in, in_int, 2, compute_device, 0x0, out, out_int ) ;

}



// Minimum
//
template<class T> __global__ 
void minimum_kernel( T* in1,T* in2, T* out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    out[idx]=min(in1[idx],in2[idx]);
  }
} 


// Minimum
//
template<class T>  
boost::shared_ptr< cuNDArray<T> >
Gadgetron::minimum( cuNDArray<T> *in1,cuNDArray<T> *in2,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  int cur_device, old_device;
  cuNDArray<T> *in1_int;
  cuNDArray<T> *in2_int;


  if ( in1->get_number_of_elements() !=  in2->get_number_of_elements()){
    BOOST_THROW_EXCEPTION(runtime_error( "minimum: input arrays have different number of elements"));

  }
  // Prepare 
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in1, &in1_int, in2, &in2_int);

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in1->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(in1->get_dimensions().get());
  if( out.get() != 0x0 ) minimum_kernel<T><<< gridDim, blockDim >>>( in1_int->get_data_ptr(), in2_int->get_data_ptr(),out->get_data_ptr(), in1->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore 
  // Restore
  restore<2,T,T,T,dummy>( old_device, in1, in1_int, 0, compute_device, out.get(), in2, in2_int );

  return out;
}



// Maximum
//
template<class T> __global__ 
void maximum_kernel( T* in1,T* in2, T* out, unsigned int number_of_elements )
{
  const unsigned int idx = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x+threadIdx.x;
 
  if( idx<number_of_elements ){
    out[idx]=max(in1[idx],in2[idx]);
  }
} 


// Minimum
//
template<class T>  
boost::shared_ptr< cuNDArray<T> >
Gadgetron::maximum( cuNDArray<T> *in1,cuNDArray<T> *in2,
	    cuNDA_device alloc_device, cuNDA_device compute_device )
{
  int cur_device, old_device;
  cuNDArray<T> *in1_int;
  cuNDArray<T> *in2_int;


  if ( in1->get_number_of_elements() !=  in2->get_number_of_elements()){
    BOOST_THROW_EXCEPTION(runtime_error( "maximum: input arrays have different number of elements"));
  }
  // Prepare 
  prepare<2,T,T,dummy>( compute_device, &cur_device, &old_device, in1, &in1_int, in2, &in2_int  ) ;

  // Setup block/grid dimensions
  dim3 blockDim; dim3 gridDim;
  setup_grid( cur_device, in1->get_number_of_elements(), &blockDim, &gridDim );

  // Invoke kernel
  boost::shared_ptr< cuNDArray<T> > out = cuNDArray<T>::allocate(in1->get_dimensions().get());
  if( out.get() != 0x0 ) maximum_kernel<T><<< gridDim, blockDim >>>( in1_int->get_data_ptr(), in2_int->get_data_ptr(),out->get_data_ptr(), in1->get_number_of_elements() );

  CHECK_FOR_CUDA_ERROR();

  // Restore 
  // Restore
  restore<2,T,T,T,dummy>( old_device, in1, in1_int, 0, compute_device, out.get(), in2, in2_int );

  return out;
}


//
// Instantiation
//

// A few functions have integer support

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<int> >
Gadgetron::sum<int>( cuNDArray<int>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<1>::Type> >
Gadgetron::sum<intd<1>::Type >( cuNDArray<intd<1>::Type >*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<2>::Type> >
Gadgetron::sum<intd<2>::Type >( cuNDArray<intd<2>::Type >*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<3>::Type> >
Gadgetron::sum<intd<3>::Type >( cuNDArray<intd<3>::Type >*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<intd<4>::Type> >
Gadgetron::sum<intd<4>::Type >( cuNDArray<intd<4>::Type >*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<unsigned int> >
Gadgetron::sum<unsigned int>( cuNDArray<unsigned int>*, unsigned int, cuNDA_device, cuNDA_device);
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<1>::Type> >
Gadgetron::sum<uintd<1>::Type>( cuNDArray<uintd<1>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<2>::Type> >
Gadgetron::sum<uintd<2>::Type>( cuNDArray<uintd<2>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<3>::Type> >
Gadgetron::sum<uintd<3>::Type>( cuNDArray<uintd<3>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<uintd<4>::Type> >
Gadgetron::sum<uintd<4>::Type>( cuNDArray<uintd<4>::Type>*, unsigned int, cuNDA_device, cuNDA_device );


template EXPORTGPUCORE void
Gadgetron::crop<int,1>( uintd1, cuNDArray<int>*, cuNDArray<int>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,1>,1>( uintd1, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,2>,1>( uintd1, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,3>,1>( uintd1, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,4>,1>( uintd1, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,1>,2>( uintd2, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,2>,2>( uintd2, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,3>,2>( uintd2, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,4>,2>( uintd2, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,1>,3>( uintd3, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,2>,3>( uintd3, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,3>,3>( uintd3, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,4>,3>( uintd3, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,1>,4>( uintd4, cuNDArray<vector_td<int,1> >*, cuNDArray<vector_td<int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,2>,4>( uintd4, cuNDArray<vector_td<int,2> >*, cuNDArray<vector_td<int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,3>,4>( uintd4, cuNDArray<vector_td<int,3> >*, cuNDArray<vector_td<int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<int,4>,4>( uintd4, cuNDArray<vector_td<int,4> >*, cuNDArray<vector_td<int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<unsigned int,1>( uintd1, cuNDArray<unsigned int>*, cuNDArray<unsigned int>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,1>,1>( uintd1, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,2>,1>( uintd1, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,3>,1>( uintd1, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,4>,1>( uintd1, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,1>,2>( uintd2, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,2>,2>( uintd2, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,3>,2>( uintd2, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,4>,2>( uintd2, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,1>,3>( uintd3, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,2>,3>( uintd3, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,3>,3>( uintd3, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,4>,3>( uintd3, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,1>,4>( uintd4, cuNDArray<vector_td<unsigned int,1> >*, cuNDArray<vector_td<unsigned int,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,2>,4>( uintd4, cuNDArray<vector_td<unsigned int,2> >*, cuNDArray<vector_td<unsigned int,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,3>,4>( uintd4, cuNDArray<vector_td<unsigned int,3> >*, cuNDArray<vector_td<unsigned int,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<unsigned int,4>,4>( uintd4, cuNDArray<vector_td<unsigned int,4> >*, cuNDArray<vector_td<unsigned int,4> >*, cuNDA_device );

// Instanciation -- single precision

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::sum<float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<complext<float> > >
Gadgetron::sum<complext<float> >( cuNDArray<complext<float> >*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<1>::Type> >
Gadgetron::sum<floatd<1>::Type>( cuNDArray<floatd<1>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<2>::Type> >
Gadgetron::sum<floatd<2>::Type>( cuNDArray<floatd<2>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<3>::Type> >
Gadgetron::sum<floatd<3>::Type>( cuNDArray<floatd<3>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<floatd<4>::Type> >
Gadgetron::sum<floatd<4>::Type>( cuNDArray<floatd<4>::Type>*, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::expand<float>( cuNDArray<float>*, unsigned int, cuNDA_device, cuNDA_device);
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
Gadgetron::expand<float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device, cuNDA_device);


template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::abs<float_complext>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::abs<float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::minimum<float>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::maximum<float>( cuNDArray<float>*, cuNDArray<float>*,cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::correlation<float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
Gadgetron::correlation<float_complext>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );



template EXPORTGPUCORE void
Gadgetron::crop<float,1>( uintd1, cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<float,2>( uintd2, cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<complext<float>,1>( uintd1, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<complext<float>,2>( uintd2, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<complext<float>,3>( uintd3, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<complext<float>,4>( uintd4, cuNDArray<complext<float> >*, cuNDArray< complext<float> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,1>,1>( uintd1, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,2>,1>( uintd1, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,3>,1>( uintd1, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,4>,1>( uintd1, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,1>,2>( uintd2, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,2>,2>( uintd2, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,3>,2>( uintd2, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,4>,2>( uintd2, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,1>,3>( uintd3, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,2>,3>( uintd3, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,3>,3>( uintd3, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,4>,3>( uintd3, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,1>,4>( uintd4, cuNDArray<vector_td<float,1> >*, cuNDArray<vector_td<float,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,2>,4>( uintd4, cuNDArray<vector_td<float,2> >*, cuNDArray<vector_td<float,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,3>,4>( uintd4, cuNDArray<vector_td<float,3> >*, cuNDArray<vector_td<float,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<float,4>,4>( uintd4, cuNDArray<vector_td<float,4> >*, cuNDArray<vector_td<float,4> >*, cuNDA_device );


template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<float,1>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<float_complext,1>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<float,2>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<float_complext,2>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<float,3>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<float_complext,3>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<float,4>( cuNDArray<float>*, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<float_complext,4>( cuNDArray<float_complext>*, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float_complext> >
Gadgetron::real_to_complext<float>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::complext_to_real<float>( cuNDArray<float_complext>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::downsample<float,1>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::downsample<float,2>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::downsample<float,3>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::downsample<float,4>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::upsample_nn<float,1>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::upsample_nn<float,2>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::upsample_nn<float,3>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::upsample_nn<float,4>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::upsample_lin<float,1>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::upsample_lin<float,2>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::upsample_lin<float,3>( cuNDArray<float>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<float> >
Gadgetron::upsample_lin<float,4>( cuNDArray<float>*, cuNDA_device, cuNDA_device );


template EXPORTGPUCORE void Gadgetron::rss_normalize<float>( cuNDArray<float>*, unsigned int, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::rss_normalize<float_complext>( cuNDArray<float_complext>*, unsigned int, cuNDA_device );


template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,1>(uintd1, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<float_complext,1>(uintd1, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,2>(uintd2, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<float_complext,2>(uintd2, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,3>(uintd3, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<float_complext,3>(uintd3, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,4>(uintd4, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<float_complext,4>(uintd4, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,float,2>(float, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,float_complext,2>(float, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,float,3>(float, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,float_complext,3>(float, cuNDArray<float_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,float,4>(float, cuNDArray<float>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<float,float_complext,4>(float, cuNDArray<float_complext>*, cuNDA_device );



template EXPORTGPUCORE void Gadgetron::shrink1<float,float>( float, cuNDArray<float>*, cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::shrink1<float,float_complext>( float, cuNDArray<float_complext>*, cuNDArray<float_complext>* );

template EXPORTGPUCORE void Gadgetron::shrinkd<float,float>( float, cuNDArray<float>*, cuNDArray<float>*, cuNDArray<float>* );
template EXPORTGPUCORE void Gadgetron::shrinkd<float,float_complext>( float, cuNDArray<float>*, cuNDArray<float_complext>*, cuNDArray<float_complext>* );

template EXPORTGPUCORE
void Gadgetron::origin_mirror<float,1>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);
template EXPORTGPUCORE
void Gadgetron::origin_mirror<float,2>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);
template EXPORTGPUCORE
void Gadgetron::origin_mirror<float,3>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);
template EXPORTGPUCORE
void Gadgetron::origin_mirror<float,4>(cuNDArray<float>*, cuNDArray<float>*, bool, cuNDA_device);

template EXPORTGPUCORE void
Gadgetron::origin_mirror<float_complext,1>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
Gadgetron::origin_mirror<float_complext,2>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
Gadgetron::origin_mirror<float_complext,3>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
Gadgetron::origin_mirror<float_complext,4>(cuNDArray<float_complext>*, cuNDArray<float_complext>*, bool, cuNDA_device);


// Instanciation -- double precision

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::sum<double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext > >
Gadgetron::sum<complext<double> >( cuNDArray<double_complext >*, unsigned int, cuNDA_device, cuNDA_device);

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<1>::Type> >
Gadgetron::sum<doubled<1>::Type>( cuNDArray<doubled<1>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<2>::Type> >
Gadgetron::sum<doubled<2>::Type>( cuNDArray<doubled<2>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<3>::Type> >
Gadgetron::sum<doubled<3>::Type>( cuNDArray<doubled<3>::Type>*, unsigned int, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<doubled<4>::Type> >
Gadgetron::sum<doubled<4>::Type>( cuNDArray<doubled<4>::Type>*, unsigned int, cuNDA_device, cuNDA_device );



template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::expand<double>( cuNDArray<double>*, unsigned int, cuNDA_device, cuNDA_device);
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
Gadgetron::expand<double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device, cuNDA_device);


template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::abs<double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );
template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::abs<double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::correlation<double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
Gadgetron::correlation<double_complext>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );

//template EXPORTGPUCORE void axpy<double>( cuNDArray<double>*, cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<double,1>( uintd1, cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<complext<double> ,1>( uintd1, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<complext<double> ,2>( uintd2, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<complext<double> ,3>( uintd3, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<complext<double> ,4>( uintd4, cuNDArray<complext<double> >*, cuNDArray< complext<double> >*, cuNDA_device );



template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,1>,1>( uintd1, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,2>,1>( uintd1, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,3>,1>( uintd1, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,4>,1>( uintd1, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,1>,2>( uintd2, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,2>,2>( uintd2, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,3>,2>( uintd2, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,4>,2>( uintd2, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,1>,3>( uintd3, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,2>,3>( uintd3, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,3>,3>( uintd3, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,4>,3>( uintd3, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,1>,4>( uintd4, cuNDArray<vector_td<double,1> >*, cuNDArray<vector_td<double,1> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,2>,4>( uintd4, cuNDArray<vector_td<double,2> >*, cuNDArray<vector_td<double,2> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,3>,4>( uintd4, cuNDArray<vector_td<double,3> >*, cuNDArray<vector_td<double,3> >*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::crop<vector_td<double,4>,4>( uintd4, cuNDArray<vector_td<double,4> >*, cuNDArray<vector_td<double,4> >*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<double,1>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<double_complext,1>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<double,2>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<double_complext,2>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<double,3>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<double_complext,3>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<double,4>( cuNDArray<double>*, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void
Gadgetron::expand_with_zero_fill<double_complext,4>( cuNDArray<double_complext>*, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double_complext> >
Gadgetron::real_to_complext<double>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::complext_to_real<double>( cuNDArray<double_complext>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::downsample<double,1>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::downsample<double,2>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::downsample<double,3>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::downsample<double,4>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::upsample_nn<double,1>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::upsample_nn<double,2>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::upsample_nn<double,3>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::upsample_nn<double,4>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::upsample_lin<double,1>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::upsample_lin<double,2>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::upsample_lin<double,3>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE boost::shared_ptr< cuNDArray<double> >
Gadgetron::upsample_lin<double,4>( cuNDArray<double>*, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::rss_normalize<double>( cuNDArray<double>*, unsigned int, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::rss_normalize<double_complext>( cuNDArray<double_complext>*, unsigned int, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,1>(uintd1, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<double_complext,1>(uintd1, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,2>(uintd2, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<double_complext,2>(uintd2, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,3>(uintd3, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<double_complext,3>(uintd3, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,4>(uintd4, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<double_complext,4>(uintd4, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,double,2>(double, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,double_complext,2>(double, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,double,3>(double, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,double_complext,3>(double, cuNDArray<double_complext>*, cuNDA_device );

template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,double,4>(double, cuNDArray<double>*, cuNDA_device );
template EXPORTGPUCORE void Gadgetron::zero_fill_border<double,double_complext,4>(double, cuNDArray<double_complext>*, cuNDA_device );


template EXPORTGPUCORE void Gadgetron::shrink1<double,double>( double, cuNDArray<double>*, cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::shrink1<double,double_complext>( double, cuNDArray<double_complext>*, cuNDArray<double_complext>* );

template EXPORTGPUCORE void Gadgetron::shrinkd<double,double>( double, cuNDArray<double>*, cuNDArray<double>*, cuNDArray<double>* );
template EXPORTGPUCORE void Gadgetron::shrinkd<double,double_complext>( double, cuNDArray<double>*, cuNDArray<double_complext>*, cuNDArray<double_complext>* );

template EXPORTGPUCORE
void Gadgetron::origin_mirror<double,1>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);
template EXPORTGPUCORE
void Gadgetron::origin_mirror<double,2>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);
template EXPORTGPUCORE
void Gadgetron::origin_mirror<double,3>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);
template EXPORTGPUCORE
void Gadgetron::origin_mirror<double,4>(cuNDArray<double>*, cuNDArray<double>*, bool, cuNDA_device);

template EXPORTGPUCORE void
Gadgetron::origin_mirror<double_complext,1>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
Gadgetron::origin_mirror<double_complext,2>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
Gadgetron::origin_mirror<double_complext,3>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);
template EXPORTGPUCORE void
Gadgetron::origin_mirror<double_complext,4>(cuNDArray<double_complext>*, cuNDArray<double_complext>*, bool, cuNDA_device);

template EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> >
Gadgetron::squaredNorm( cuNDArray<float> *, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE
boost::shared_ptr< cuNDArray<float> >
Gadgetron::squaredNorm( cuNDArray<float_complext> *, unsigned int ,cuNDA_device , cuNDA_device);

template EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> >
Gadgetron::squaredNorm( cuNDArray<double> *, unsigned int, cuNDA_device, cuNDA_device );

template EXPORTGPUCORE
boost::shared_ptr< cuNDArray<double> >
Gadgetron::squaredNorm( cuNDArray<double_complext> *, unsigned int ,cuNDA_device , cuNDA_device);

template EXPORTGPUCORE
void Gadgetron::inplace_sgn(cuNDArray<float>* x);
template EXPORTGPUCORE
void Gadgetron::inplace_sgn(cuNDArray<double>* x);


template void Gadgetron::clamp<double>(cuNDArray<double>*,double,double);
template void Gadgetron::clamp<float>(cuNDArray<float>*,float,float);

template void Gadgetron::clamp_min<double>(cuNDArray<double>*,double);
template void Gadgetron::clamp_min<float>(cuNDArray<float>*,float);

template void Gadgetron::clamp_max<double>(cuNDArray<double>*,double);
template void Gadgetron::clamp_max<float>(cuNDArray<float>*,float);

template void Gadgetron::clamp_min<double_complext>(cuNDArray<double_complext>*,double);
template void Gadgetron::clamp_min<float_complext>(cuNDArray<float_complext>*,float);



template float Gadgetron::normalize<float>( cuNDArray<float> *, float, cuNDA_device );
template double Gadgetron::normalize<double>( cuNDArray<double> *, double, cuNDA_device );
