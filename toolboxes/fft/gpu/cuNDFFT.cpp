#include "cuNDFFT.h"
#include "vector_td.h"
#include "cuNDArray.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_operators.h"

#include <cufft.h>
#include <cuComplex.h>
#include <sstream>

namespace Gadgetron{

  template<class T> cuNDFFT<T>* cuNDFFT<T>::instance()
  {
    if (!__instance)
	  __instance = new cuNDFFT<T>;
	  return __instance;
  }

  template<class T> cuNDFFT<T>* cuNDFFT<T>::__instance = NULL;

  template<class T> cufftType_t get_transform_type();
  template<> cufftType_t get_transform_type<float>() { return CUFFT_C2C; }
  template<> cufftType_t get_transform_type<double>() { return CUFFT_Z2Z; }
  
  template<class T> cufftResult_t cuNDA_FFT_execute( cufftHandle plan, cuNDArray< complext<T> > *in_out, int direction );
  
  template<> cufftResult_t cuNDA_FFT_execute<float>( cufftHandle plan, cuNDArray<float_complext> *in_out, int direction ){
    return cufftExecC2C(plan, (cuFloatComplex*)in_out->get_data_ptr(), (cuFloatComplex*)in_out->get_data_ptr(), direction); }

  template<> cufftResult_t cuNDA_FFT_execute<double>( cufftHandle plan, cuNDArray<double_complext> *in_out, int direction ){
    return cufftExecZ2Z(plan, (cuDoubleComplex*)in_out->get_data_ptr(), (cuDoubleComplex*)in_out->get_data_ptr(), direction); }
  
  template<class T> void
  cuNDFFT<T>::fft_int( cuNDArray< complext<T> > *input, std::vector<size_t> *dims_to_transform, int direction, bool do_scale )
  {
    std::vector<size_t> new_dim_order;
    std::vector<size_t> reverse_dim_order;
    std::vector<size_t> dims;
    std::vector<size_t> dim_count(input->get_number_of_dimensions(),0);
    
    size_t array_ndim = input->get_number_of_dimensions();
    boost::shared_ptr< std::vector<size_t> > array_dims = input->get_dimensions();
    
    dims = std::vector<size_t>(dims_to_transform->size(),0);
    for (size_t i = 0; i < dims_to_transform->size(); i++) {
      if ((*dims_to_transform)[i] >= array_ndim) {
    	std::stringstream ss;
    	ss << "cuNDFFT::fft Invalid dimensions specified for transform " << (*dims_to_transform)[i] << "max " << array_ndim;
	throw std::runtime_error(ss.str());;
      }
      if (dim_count[(*dims_to_transform)[i]] > 0) {
	throw std::runtime_error("cuNDFFT::fft Invalid dimensions (duplicates) specified for transform");;
      }
      dim_count[(*dims_to_transform)[i]]++;
      dims[dims_to_transform->size()-1-i] = (*array_dims)[(*dims_to_transform)[i]];
    }
    
    new_dim_order = *dims_to_transform;
    for (size_t i = 0; i < array_ndim; i++) {
      if (!dim_count[i]) new_dim_order.push_back(i);
    }
    
    reverse_dim_order = std::vector<size_t>(array_ndim,0);
    for (size_t i = 0; i < array_ndim; i++) {
      reverse_dim_order[new_dim_order[i]] = i;
    }
    
    size_t ndim = dims.size();
    size_t batches = 0;
    size_t elements_in_ft = 1;
    for (size_t i = 0; i < dims.size(); i++)
      elements_in_ft *= dims[i];
    batches = input->get_number_of_elements() / elements_in_ft;
    
    cufftHandle plan;
    cufftResult ftres;
    
    std::vector<int> int_dims;
    for( unsigned int i=0; i<dims.size(); i++ )
      int_dims.push_back((int)dims[i]);

    ftres = cufftPlanMany(&plan,ndim,&int_dims[0], &int_dims[0], 1, elements_in_ft, &int_dims[0], 1, elements_in_ft, get_transform_type<T>(), batches);
    if (ftres != CUFFT_SUCCESS) {
      std::stringstream ss;
      ss << "cuNDFFT FFT plan failed: " << ftres;
      throw std::runtime_error(ss.str());;
    }
    
    bool must_permute = false;

    {
    	for (size_t i = 0; i < new_dim_order.size(); i++)
    		must_permute |= (i != new_dim_order[i]);
    }

    if (must_permute)
    	*input = *permute(input,&new_dim_order);

    if (direction == CUFFT_INVERSE)
    	for (size_t i =0; i < dims_to_transform->size(); i++)
    	timeswitch(input,dims_to_transform->at(i));
    
    if( cuNDA_FFT_execute<T>( plan, input, direction ) != CUFFT_SUCCESS ) {
      throw std::runtime_error("cuNDFFT FFT execute failed");;
    }
    
    ftres = cufftDestroy( plan );
    if (ftres != CUFFT_SUCCESS) {
      std::stringstream ss;
      ss << "cuNDFFT FFT plan destroy failed: " << ftres;
      throw std::runtime_error(ss.str());;
    }
    if (direction == CUFFT_FORWARD)
    	for (size_t i =0; i < dims_to_transform->size(); i++)
    		timeswitch(input,dims_to_transform->at(i));
    
    if (do_scale) {
      *input *= 1/std::sqrt(T(elements_in_ft));
    }
    
    if (must_permute)
    	*input = *permute(input,&reverse_dim_order);
  }
  
  template<class T> void
  cuNDFFT<T>::fft( cuNDArray< complext<T> > *input, std::vector<size_t> *dims_to_transform, bool do_scale )
  {
    fft_int(input, dims_to_transform, CUFFT_FORWARD, do_scale);
  }
  
  template<class T> void
  cuNDFFT<T>::ifft( cuNDArray< complext<T> > *input, std::vector<size_t> *dims_to_transform, bool do_scale )
  {
    fft_int(input, dims_to_transform, CUFFT_INVERSE, do_scale);
  }
  
  template<class T> void
  cuNDFFT<T>::fft( cuNDArray< complext<T> > *input, unsigned int dim_to_transform, bool do_scale )
  {
    std::vector<size_t> dims(1,dim_to_transform);
    fft_int(input, &dims, CUFFT_FORWARD, do_scale);
  }
  
  template<class T> void
  cuNDFFT<T>::ifft( cuNDArray< complext<T> > *input, unsigned int dim_to_transform, bool do_scale )
  {
    std::vector<size_t> dims(1,dim_to_transform);
    fft_int(input, &dims, CUFFT_INVERSE, do_scale);
  }
  
  template<class T> void
  cuNDFFT<T>::fft( cuNDArray< complext<T> > *input, bool do_scale )
  {
    std::vector<size_t> dims(input->get_number_of_dimensions(),0);
    for (size_t i = 0; i < dims.size(); i++) dims[i] = i;
    fft_int(input, &dims, CUFFT_FORWARD, do_scale);
  }
  
  template<class T> void
  cuNDFFT<T>::ifft( cuNDArray<complext<T> > *input, bool do_scale )
  {
    std::vector<size_t> dims(input->get_number_of_dimensions(),0);
    for (size_t i = 0; i < dims.size(); i++) dims[i] = i;
    fft_int(input, &dims, CUFFT_INVERSE, do_scale);
  }
  
  // Instantiation
  template class EXPORTGPUFFT cuNDFFT<float>;
  template class EXPORTGPUFFT cuNDFFT<double>;
}
