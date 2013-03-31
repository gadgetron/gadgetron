#include "cuFFT.h"
#include "vector_td.h"
#include "cuNDArray.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_operators.h"

#include <cufft.h>
#include <cublas.h>
#include <cuComplex.h>
#include <sstream>

namespace Gadgetron{
  
  template<class T> cufftType_t get_transform_type();
  template<> cufftType_t get_transform_type< cuFloatComplex  >() { return CUFFT_C2C; }
  template<> cufftType_t get_transform_type< cuDoubleComplex >() { return CUFFT_Z2Z; }
  template<> cufftType_t get_transform_type< float_complext  >() { return CUFFT_C2C; }
  template<> cufftType_t get_transform_type< double_complext >() { return CUFFT_Z2Z; }
  
  template<class T> cufftResult_t cuNDA_FFT_execute( cufftHandle plan, cuNDArray<T> *in_out, int direction );

  template<> cufftResult_t cuNDA_FFT_execute<cuFloatComplex>( cufftHandle plan, cuNDArray<cuFloatComplex> *in_out, int direction ){
    return cufftExecC2C(plan, in_out->get_data_ptr(), in_out->get_data_ptr(), direction); }

  template<> cufftResult_t cuNDA_FFT_execute<cuDoubleComplex>( cufftHandle plan, cuNDArray<cuDoubleComplex> *in_out, int direction ){
    return cufftExecZ2Z(plan, in_out->get_data_ptr(), in_out->get_data_ptr(), direction); }

  template<> cufftResult_t cuNDA_FFT_execute<float_complext>( cufftHandle plan, cuNDArray<float_complext> *in_out, int direction ){
    return cufftExecC2C(plan, (cuFloatComplex*)in_out->get_data_ptr(), (cuFloatComplex*)in_out->get_data_ptr(), direction); }

  template<> cufftResult_t cuNDA_FFT_execute<double_complext>( cufftHandle plan, cuNDArray<double_complext> *in_out, int direction ){
    return cufftExecZ2Z(plan, (cuDoubleComplex*)in_out->get_data_ptr(), (cuDoubleComplex*)in_out->get_data_ptr(), direction); }
  
  template<class T> void
  cuFFT<T>::fft_int( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform, int direction, bool do_scale )
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
    	std::stringstream ss;
    	ss << "cuFFT::fft Invalid dimensions specified for transform " << (*dims_to_transform)[i] << "max " << array_ndim;
	BOOST_THROW_EXCEPTION(runtime_error(ss.str()));	
      }
      if (dim_count[(*dims_to_transform)[i]] > 0) {
	BOOST_THROW_EXCEPTION(runtime_error("cuFFT::fft Invalid dimensions (duplicates) specified for transform"));	
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
      std::stringstream ss;
      ss << "cuFFT FFT plan failed: " << ftres;
      BOOST_THROW_EXCEPTION(runtime_error(ss.str()));      
    }
    
    //IFFTSHIFT
    *input = *permute(input,&new_dim_order,-1);
    
    if( cuNDA_FFT_execute<T>( plan, input, direction ) != CUFFT_SUCCESS ) {
      BOOST_THROW_EXCEPTION(runtime_error("cuFFT FFT execute failed"));      
    }
    
    ftres = cufftDestroy( plan );
    if (ftres != CUFFT_SUCCESS) {
      std::stringstream ss;
      ss << "cuFFT FFT plan destroy failed: " << ftres;
      BOOST_THROW_EXCEPTION(runtime_error(ss.str()));      
    }
    
    if (do_scale) {
      *input /= T(elements_in_ft);      
    }
    
    //FFTSHIFT 
    *input = *permute(input,&reverse_dim_order,1);        
  }
  
  template<class T> void
  cuFFT<T>::fft( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform )
  {
    fft_int(input, dims_to_transform, CUFFT_FORWARD, false);
  }
  
  template<class T> void
  cuFFT<T>::ifft( cuNDArray<T> *input, std::vector<unsigned int> *dims_to_transform, bool do_scale )
  {
    fft_int(input, dims_to_transform, CUFFT_INVERSE, do_scale);
  }
  
  template<class T> void
  cuFFT<T>::fft( cuNDArray<T> *input, unsigned int dim_to_transform )
  {
    std::vector<unsigned int> dims(1,dim_to_transform);
    fft_int(input, &dims, CUFFT_FORWARD, false);
  }
  
  template<class T> void
  cuFFT<T>::ifft( cuNDArray<T> *input, unsigned int dim_to_transform, bool do_scale )
  {
    std::vector<unsigned int> dims(1,dim_to_transform);
    fft_int(input, &dims, CUFFT_INVERSE, do_scale);
  }
  
  template<class T> void
  cuFFT<T>::fft( cuNDArray<T> *input )
  {
    std::vector<unsigned int> dims(input->get_number_of_dimensions(),0);
    for (unsigned int i = 0; i < dims.size(); i++) dims[i] = i;
    fft_int(input, &dims, CUFFT_FORWARD, false);
  }
  
  template<class T> void
  cuFFT<T>::ifft( cuNDArray<T> *input, bool do_scale )
  {
    std::vector<unsigned int> dims(input->get_number_of_dimensions(),0);
    for (unsigned int i = 0; i < dims.size(); i++) dims[i] = i;
    fft_int(input, &dims, CUFFT_INVERSE, do_scale);
  }
  
  // Instantiation
  //template class EXPORTGPUCORE cuFFT<cuFloatComplex>;
  //template class EXPORTGPUCORE cuFFT<cuDoubleComplex>;
  template class EXPORTGPUCORE cuFFT<float_complext>;
  template class EXPORTGPUCORE cuFFT<double_complext>;
}
