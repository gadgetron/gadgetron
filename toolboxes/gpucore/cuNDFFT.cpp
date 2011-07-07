#include "cuNDFFT.h"
#include "vector_td.h"
#include "ndarray_vector_td_utilities.h"

#include <cufft.h>
#include <cublas.h>
#include <cuComplex.h>

template<class T> cufftType_t get_transform_type();
template<> cufftType_t get_transform_type< cuFloatComplex  >() { return CUFFT_C2C; }
template<> cufftType_t get_transform_type< cuDoubleComplex >() { return CUFFT_Z2Z; }
template<> cufftType_t get_transform_type< float_complext::Type  >() { return CUFFT_C2C; }
template<> cufftType_t get_transform_type< double_complext::Type >() { return CUFFT_Z2Z; }

template<class T> cufftResult_t cuNDA_FFT_execute( cufftHandle plan, cuNDArray<T> *in_out, int direction );
template<> cufftResult_t cuNDA_FFT_execute<cuFloatComplex>( cufftHandle plan, cuNDArray<cuFloatComplex> *in_out, int direction ){
  return cufftExecC2C(plan, in_out->get_data_ptr(), in_out->get_data_ptr(), direction); }
template<> cufftResult_t cuNDA_FFT_execute<cuDoubleComplex>( cufftHandle plan, cuNDArray<cuDoubleComplex> *in_out, int direction ){
  return cufftExecZ2Z(plan, in_out->get_data_ptr(), in_out->get_data_ptr(), direction); }
template<> cufftResult_t cuNDA_FFT_execute<float_complext::Type>( cufftHandle plan, cuNDArray<float_complext::Type> *in_out, int direction ){
  return cufftExecC2C(plan, (cuFloatComplex*)in_out->get_data_ptr(), (cuFloatComplex*)in_out->get_data_ptr(), direction); }
template<> cufftResult_t cuNDA_FFT_execute<double_complext::Type>( cufftHandle plan, cuNDArray<double_complext::Type> *in_out, int direction ){
  return cufftExecZ2Z(plan, (cuDoubleComplex*)in_out->get_data_ptr(), (cuDoubleComplex*)in_out->get_data_ptr(), direction); }

template<class T> bool cuNDA_FFT_scal( T a, cuNDArray<T> *in_out );
template<> bool cuNDA_FFT_scal( float_complext::Type a, cuNDArray<float_complext::Type> *in_out ){
  return cuNDA_scal<float_complext::Type>( a, in_out ); }
template<> bool cuNDA_FFT_scal( double_complext::Type a, cuNDArray<double_complext::Type> *in_out ){
  return cuNDA_scal<double_complext::Type>( a, in_out ); }
template<> bool cuNDA_FFT_scal( cuFloatComplex a, cuNDArray<cuFloatComplex> *in_out ){
  return cuNDA_scal<float_complext::Type>( *((float_complext::Type*)&a), (cuNDArray<float_complext::Type>*)in_out ); }
template<> bool cuNDA_FFT_scal( cuDoubleComplex a, cuNDArray<cuDoubleComplex> *in_out ){
  return cuNDA_scal<double_complext::Type>( *((double_complext::Type*)&a), (cuNDArray<double_complext::Type>*)in_out ); }

template<class T> T operator* ( const T &z, const int &i ); 
template<> cuFloatComplex operator* ( const cuFloatComplex &z, const int &i ){ 
  cuFloatComplex res = z; res.x*=(float)i; res.y*=(float)i; return res; }
template<> cuDoubleComplex operator* ( const cuDoubleComplex &z, const int &i ){ 
  cuDoubleComplex res = z; res.x*=(double)i; res.y*=(double)i; return res; }
template<> float_complext::Type operator* ( const float_complext::Type &z, const int &i ){ 
  float_complext::Type res = z; res.vec[0]*=(float)i; res.vec[1]*=(float)i; return res; }
template<> double_complext::Type operator* ( const double_complext::Type &z, const int &i ){ 
  double_complext::Type res = z; res.vec[0]*=(double)i; res.vec[1]*=(double)i; return res; }

template<> cuFloatComplex get_one<cuFloatComplex>(){
  cuFloatComplex res; res.x = 1.0f; res.y = 0.0f; return res; }
template<> cuDoubleComplex get_one<cuDoubleComplex>(){
  cuDoubleComplex res; res.x = 1.0f; res.y = 0.0f; return res; }

template<> cuFloatComplex reciprocal<cuFloatComplex>(const cuFloatComplex z){
  return cuCdivf( make_cuFloatComplex(1.0f,0.0f), z );
}
template<> cuDoubleComplex reciprocal<cuDoubleComplex>(const cuDoubleComplex z){
  return cuCdiv( make_cuDoubleComplex(1.0,0.0), z );
}


template<class T> int 
cuNDFFT<T>::fft_int( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform, int direction, bool do_scale )
{
  std::vector<unsigned int> new_dim_order;
  std::vector<unsigned int> reverse_dim_order;
  std::vector<int> dims;
  std::vector<unsigned int> dim_count(input->get_number_of_dimensions(),0);
  
  unsigned int array_ndim = input->get_number_of_dimensions();
  boost::shared_ptr< std::vector<unsigned int> > array_dims = input->get_dimensions();
  
  dims = std::vector<int>(dims_to_transform->size(),0);
  for (unsigned int i = 0; i < dims_to_transform->size(); i++) {
    if ((*dims_to_transform)[i] >= array_ndim) {
      std::cerr << "cuNDFFT::fft Invalid dimensions specified for transform " << (*dims_to_transform)[i] << "max " << array_ndim << std::endl;
      return -1;
    }
    if (dim_count[(*dims_to_transform)[i]] > 0) {
      std::cerr << "cuNDFFT::fft Invalid dimensions (duplicates) specified for transform" << std::endl;
      return -1;
    }
    dim_count[(*dims_to_transform)[i]]++;
    dims[dims_to_transform->size()-1-i] = (*array_dims)[(*dims_to_transform)[i]];
  }
  
  new_dim_order = *dims_to_transform;
  for (unsigned int i = 0; i < array_ndim; i++) {
    if (!dim_count[i]) new_dim_order.push_back(i);
  }

  reverse_dim_order = std::vector<unsigned int>(array_ndim,0);
  for (unsigned int i = 0; i < array_ndim; i++) {
    reverse_dim_order[new_dim_order[i]] = i;
  }

  int ndim = dims.size();
  int batches = 0;
  int elements_in_ft = 1;
  for (unsigned int i = 0; i < dims.size(); i++) 
    elements_in_ft *= dims[i];
  batches = input->get_number_of_elements() / elements_in_ft;

  cufftHandle plan;
  cufftResult ftres;

  ftres = cufftPlanMany(&plan,ndim,&dims[0],&dims[0],1,elements_in_ft,&dims[0],1,elements_in_ft,get_transform_type<T>(),batches);
  if (ftres != CUFFT_SUCCESS) {
    std::cerr << "cuNDFFT FFT plan failed: " << ftres << std::endl;
    return -1;
  }

  //IFFTSHIFT
  if (input->permute(&new_dim_order,0,-1) < 0) {
    std::cerr << "cuNDFFT error permuting before FFT" << std::endl;
    return -1;
  }

  if( cuNDA_FFT_execute<T>( plan, input, direction ) != CUFFT_SUCCESS ) {
    std::cerr << "cuNDFFT FFT execute failed" << std::endl;
    return -1;
  }

  ftres = cufftDestroy( plan );
  if (ftres != CUFFT_SUCCESS) {
    std::cerr << "cuNDFFT FFT plan destroy failed: " << ftres << std::endl;
    return -1;
  }

  if (do_scale) {
    T scale = reciprocal<T>(get_one<T>()*elements_in_ft);
    if( !cuNDA_FFT_scal( scale, input ) ){
      std::cerr << "cuNDFFT rescaling failed " << std::endl;
      return -1;
    } 
  }
  
  //FFTSHIFT 
  if (input->permute(&reverse_dim_order,0,1) < 0) {
    std::cerr << "cuNDFFT error permuting after FFT" << std::endl;
    return -1;
  }
  
  return 0;
}

template<class T> int 
cuNDFFT<T>::fft( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform )
{
  return fft_int(input, dims_to_transform, CUFFT_FORWARD, false);
}

template<class T> int
cuNDFFT<T>::ifft( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform, bool do_scale )
{
  return fft_int(input, dims_to_transform, CUFFT_INVERSE, do_scale);
}


template<class T> int 
cuNDFFT<T>::fft( cuNDArray<T> *input, unsigned int dim_to_transform )
{
  std::vector<unsigned int> dims(1,dim_to_transform);
  return fft_int(input, &dims, CUFFT_FORWARD, false);
}
  
template<class T> int
cuNDFFT<T>::ifft( cuNDArray<T> *input, unsigned int dim_to_transform, bool do_scale )
{
  std::vector<unsigned int> dims(1,dim_to_transform);
  return fft_int(input, &dims, CUFFT_INVERSE, do_scale);
}

template<class T> int
cuNDFFT<T>::fft( cuNDArray<T> *input )
{
  std::vector<unsigned int> dims(input->get_number_of_dimensions(),0);
  for (unsigned int i = 0; i < dims.size(); i++) dims[i] = i;
  return fft_int(input, &dims, CUFFT_FORWARD, false);
}

template<class T> int
cuNDFFT<T>::ifft( cuNDArray<T> *input, bool do_scale )
{
  std::vector<unsigned int> dims(input->get_number_of_dimensions(),0);
  for (unsigned int i = 0; i < dims.size(); i++) dims[i] = i;
  return fft_int(input, &dims, CUFFT_INVERSE, do_scale);
}

// Instantiation
template class EXPORTGPUCORE cuNDFFT<cuFloatComplex>;
template class EXPORTGPUCORE cuNDFFT<cuDoubleComplex>;
template class EXPORTGPUCORE cuNDFFT<float_complext::Type>;
template class EXPORTGPUCORE cuNDFFT<double_complext::Type>;
