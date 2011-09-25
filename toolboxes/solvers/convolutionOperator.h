#pragma once

#include "matrixOperator.h"
#include "vector_td_utilities.h"

#include <boost/smart_ptr.hpp>
#include <vector>

template <class REAL, class COMPLEX_ARRAY_TYPE> class convolutionOperator : public matrixOperator<REAL, COMPLEX_ARRAY_TYPE>
{
  
public:
  
  convolutionOperator() : matrixOperator<REAL,COMPLEX_ARRAY_TYPE>() {}
  virtual ~convolutionOperator() {}
  
  virtual bool operator_fft( bool forwards_transform, COMPLEX_ARRAY_TYPE *image ) = 0;
  virtual bool operator_scale( bool conjugate_kernel, COMPLEX_ARRAY_TYPE *kernel, COMPLEX_ARRAY_TYPE *image ) = 0;
  virtual bool operator_xpy( COMPLEX_ARRAY_TYPE*, COMPLEX_ARRAY_TYPE* ) = 0;

  // Set the convolution kernel. Has to be of type complext even if COMPLEX_ARRAY_TYPE is real.
  virtual bool set_kernel( COMPLEX_ARRAY_TYPE *image_space_kernel )
  {     
    COMPLEX_ARRAY_TYPE *freq_kernel = new COMPLEX_ARRAY_TYPE(*image_space_kernel);
    if( !operator_fft( true, freq_kernel ) ){
      std::cout << std::endl << "convolutionOperator::set_kernel : fft failed" << std::endl;
      return false;
    }
    
    kernel_ = boost::shared_ptr<COMPLEX_ARRAY_TYPE>(freq_kernel);
    return true;
  }
  
  // Apply regularization image operator
  virtual int mult_MH_M( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out, bool accumulate = false )
  {    
    if( !kernel_.get() ){
      std::cout << std::endl << "convolutionOperator::mult_MH_M failed : kernel is not set" << std::endl;
      return -1;
    }
    
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() || in->get_number_of_elements() != kernel_->get_number_of_elements() ){
      std::cout << std::endl << "convolutionOperator::mult_MH_M failed : in/out image mismatch the kernel" << std::endl;
      return -1;
    }

    // TODO: optimize this with respect to memory bandwidth:
    // If !accumulate we can use 'out' over 'tmp'. Combine the two scaling operations to one.

    // Make a copy of the input
    COMPLEX_ARRAY_TYPE tmp(*in); 
    
    if( !tmp.get_data_ptr() ){
      std::cout << std::endl << "convolutionOperator::mult_MH_M failed : intermediate memory allocation failed" << std::endl;
      return -1;
    }
    
    // Forwards fft
    operator_fft( true, &tmp );

    // Multiply
    operator_scale( false, kernel_.get(), &tmp );
    operator_scale( true, kernel_.get(), &tmp );

    // Inverse fft
    operator_fft( false, &tmp );
    
    if( !accumulate ) 
      *out = tmp;
    else
      operator_xpy( &tmp, out );
    
    return 0;
  }
  
  
  virtual int mult_M( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out, bool accumulate = false )
  {
    if( !kernel_.get() ){
      std::cout << std::endl << "convolutionOperator::mult_M failed : kernel is not set" << std::endl;
      return -1;
    }
    
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() || in->get_number_of_elements() != kernel_->get_number_of_elements() ){
      std::cout << std::endl << "convolutionOperator::mult_M failed : in/out image mismatch the kernel" << std::endl;
      return -1;
    }

    // Make a copy of the input
    COMPLEX_ARRAY_TYPE tmp(*in); 
    
    if( !tmp.get_data_ptr() ){
      std::cout << std::endl << "convolutionOperator::mult_M failed : intermediate memory allocation failed" << std::endl;
      return -1;
    }
    
    // Forwards fft
    operator_fft( true, &tmp );

    // Multiply
    operator_scale( false, kernel_.get(), &tmp );

    // Inverse fft
    operator_fft( false, &tmp );
    
    if( !accumulate ) 
      *out = tmp;
    else
      operator_xpy( &tmp, out );
    
    return 0;
  }
  
  virtual int mult_MH( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out, bool accumulate = false )
  {
    if( !kernel_.get() ){
      std::cout << std::endl << "convolutionOperator::mult_MH failed : kernel is not set" << std::endl;
      return -1;
    }
    
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() || in->get_number_of_elements() != kernel_->get_number_of_elements() ){
      std::cout << std::endl << "convolutionOperator::mult_MH failed : in/out image mismatch the kernel" << std::endl;
      return -1;
    }

    // Make a copy of the input
    COMPLEX_ARRAY_TYPE tmp(*in); 
    
    if( !tmp.get_data_ptr() ){
      std::cout << std::endl << "convolutionOperator::mult_MH_M failed : intermediate memory allocation failed" << std::endl;
      return -1;
    }
    
    // Inverse fft
    operator_fft( false, &tmp );

    // Multiply
    operator_scale( true, kernel_.get(), &tmp );

   // Forwards fft
    operator_fft( true, &tmp );
    
    if( !accumulate ) 
      *out = tmp;
    else
      operator_xpy( &tmp, out );
    
    return 0;
  }
  
private:
  boost::shared_ptr<COMPLEX_ARRAY_TYPE> kernel_;
};
