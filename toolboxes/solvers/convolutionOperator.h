#pragma once

#include "linearOperator.h"
#include "vector_td_utilities.h"

#include <boost/smart_ptr.hpp>
#include <vector>

template <class REAL, class COMPLEX_ARRAY_TYPE, unsigned int D> class convolutionOperator 
	: public linearOperator<REAL, COMPLEX_ARRAY_TYPE>
{  
public:
  
  convolutionOperator() : linearOperator<REAL,COMPLEX_ARRAY_TYPE>() {}
  virtual ~convolutionOperator() {}

  // The functionality of these pure virtual functions must be provided in a derived class
  virtual bool operator_fft( bool forwards_transform, COMPLEX_ARRAY_TYPE *image ) = 0;
  virtual bool operator_scale( COMPLEX_ARRAY_TYPE *x, COMPLEX_ARRAY_TYPE *y, bool conjugate_kernel = false ) = 0; // Element-wise mult 
  virtual bool operator_xpy( COMPLEX_ARRAY_TYPE *x, COMPLEX_ARRAY_TYPE *y ) = 0; // x plus y
  virtual bool operator_mirror( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out ) = 0;
  virtual bool operator_expand( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out ) = 0;
  virtual bool operator_crop( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out ) = 0;

  // Set the convolution kernel
  virtual bool set_kernel( COMPLEX_ARRAY_TYPE *image_space_kernel )
  {     
    COMPLEX_ARRAY_TYPE *freq_kernel = new COMPLEX_ARRAY_TYPE(*image_space_kernel);
    if( !freq_kernel || !operator_fft( true, freq_kernel ) ){
      std::cout << std::endl << "convolutionOperator::set_kernel : fft failed" << std::endl;
      return false;
    }
    kernel_ = boost::shared_ptr<COMPLEX_ARRAY_TYPE>(freq_kernel);
 
    COMPLEX_ARRAY_TYPE *freq_kernel_adjoint = new COMPLEX_ARRAY_TYPE();
    if( !freq_kernel_adjoint->create( freq_kernel->get_dimensions().get() )){
      std::cout << std::endl << "convolutionOperator::set_kernel : adjoint kernel allocation failed" << std::endl;
      return false;
    }
    
    if( !operator_mirror( freq_kernel, freq_kernel_adjoint )){
      std::cout << std::endl << "convolutionOperator::set_kernel : mirroring adjoint kernel failed" << std::endl;
      return false;
    }
    
    adjoint_kernel_ = boost::shared_ptr<COMPLEX_ARRAY_TYPE>(freq_kernel_adjoint);

    return true;
  }
  
  // Apply image operators
  //

  virtual int mult_MH_M( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out, bool accumulate = false )
  {    
    if( !kernel_.get() ){
      std::cout << std::endl << "convolutionOperator::mult_MH_M failed : kernel is not set" << std::endl;
      return -1;
    }
    
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
      std::cout << std::endl << "convolutionOperator::mult_MH_M failed : in/out image dimensions mismatch" << std::endl;
      return -1;
    }

    bool use_oversampling;
    if( in->get_number_of_elements() == kernel_->get_number_of_elements() )
      use_oversampling = false;
    else if( (in->get_number_of_elements()<<D) == kernel_->get_number_of_elements() )
      use_oversampling = true;
    else{
      std::cout << std::endl << "convolutionOperator::mult_MH_M failed : in/out image dimensions mismatch the kernel" << std::endl;
      return -1;
    }
    
    // Intermediate variables
    COMPLEX_ARRAY_TYPE *tmp_out;

    if( use_oversampling ){

      boost::shared_ptr< std::vector<unsigned int> > osdims = kernel_->get_dimensions();

      tmp_out = new COMPLEX_ARRAY_TYPE();
      if( !tmp_out->create(osdims.get()) ){
	std::cout << std::endl << "convolutionOperator::mult_MH_M failed : memory allocation failed (1)" << std::endl;
	return -1;
      } 

      if( !operator_expand( in, tmp_out ) ){
	std::cout << std::endl << "convolutionOperator::mult_MH_M failed : cannot expand input image" << std::endl;
	return -1;
      }     
    }
    else if( accumulate ){
      tmp_out = new COMPLEX_ARRAY_TYPE(*in);
      if( !tmp_out->get_data_ptr() ){
	std::cout << std::endl << "convolutionOperator::mult_MH_M failed : memory allocation failed (2)" << std::endl;
	return -1;
      }
    }
    else{ 
      *out = *in;
      tmp_out = out;
    } 

    // Forwards fft
    operator_fft( true, tmp_out );

    // Multiply
    operator_scale( kernel_.get(), tmp_out );
    operator_scale( adjoint_kernel_.get(), tmp_out );

    // Inverse fft
    operator_fft( false, tmp_out );

    if( use_oversampling ) {
      if( !operator_crop( tmp_out, out )){
	std::cout << std::endl << "convolutionOperator::mult_MH_M failed : cannot crop oversampled image" << std::endl;
	return -1;
      }
      delete tmp_out;
    }    
    else if( accumulate ){
      operator_xpy( tmp_out, out );
      delete tmp_out;
    }
    
    return 0;
  }
  
  
  virtual int mult_M( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out, bool accumulate = false )
  {
    if( !kernel_.get() ){
      std::cout << std::endl << "convolutionOperator::mult_M failed : kernel is not set" << std::endl;
      return -1;
    }
    
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
      std::cout << std::endl << "convolutionOperator::mult_M failed : in/out image dimensions mismatch" << std::endl;
      return -1;
    }

    bool use_oversampling;
    if( in->get_number_of_elements() == kernel_->get_number_of_elements() )
      use_oversampling = false;
    else if( (in->get_number_of_elements()<<D) == kernel_->get_number_of_elements() )
      use_oversampling = true;
    else{
      std::cout << std::endl << "convolutionOperator::mult_M failed : in/out image dimensions mismatch the kernel" << std::endl;
      return -1;
    }
    
    // Intermediate variables
    COMPLEX_ARRAY_TYPE *tmp_out;

    if( use_oversampling ){

      boost::shared_ptr< std::vector<unsigned int> > osdims = kernel_->get_dimensions();

      tmp_out = new COMPLEX_ARRAY_TYPE();
      if( !tmp_out->create(osdims.get()) ){
	std::cout << std::endl << "convolutionOperator::mult_M failed : memory allocation failed (1)" << std::endl;
	return -1;
      } 

      if( !operator_expand( in, tmp_out ) ){
	std::cout << std::endl << "convolutionOperator::mult_M failed : cannot expand input image" << std::endl;
	return -1;
      }
    }
    else if( accumulate ){
      tmp_out = new COMPLEX_ARRAY_TYPE(*in);
      if( !tmp_out->get_data_ptr() ){
	std::cout << std::endl << "convolutionOperator::mult_M failed : memory allocation failed (2)" << std::endl;
	return -1;
      }
    }
    else{ 
      *out = *in;
      tmp_out = out;
    } 

    // Forwards fft
    operator_fft( true, tmp_out );

    // Multiply
    operator_scale( kernel_.get(), tmp_out );
 
    // Inverse fft
    operator_fft( false, tmp_out );

    if( use_oversampling ) {
      if( !operator_crop( tmp_out, out )){
	std::cout << std::endl << "convolutionOperator::mult_M failed : cannot crop oversampled image" << std::endl;
	return -1;
      }
      delete tmp_out;
    }    
    else if( accumulate ){
      operator_xpy( tmp_out, out );
      delete tmp_out;
    }
    
    return 0;
  }
  
  virtual int mult_MH( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out, bool accumulate = false )
  {
    if( !adjoint_kernel_.get() ){
      std::cout << std::endl << "convolutionOperator::mult_MH failed : kernel is not set" << std::endl;
      return -1;
    }
    
    if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
      std::cout << std::endl << "convolutionOperator::mult_MH failed : in/out image dimensions mismatch" << std::endl;
      return -1;
    }

    bool use_oversampling;
    if( in->get_number_of_elements() == adjoint_kernel_->get_number_of_elements() )
      use_oversampling = false;
    else if( (in->get_number_of_elements()<<D) == adjoint_kernel_->get_number_of_elements() )
      use_oversampling = true;
    else{
      std::cout << std::endl << "convolutionOperator::mult_MH failed : in/out image dimensions mismatch the kernel" << std::endl;
      return -1;
    }
    
    // Intermediate variables
    COMPLEX_ARRAY_TYPE *tmp_out;

    if( use_oversampling ){

      boost::shared_ptr< std::vector<unsigned int> > osdims = adjoint_kernel_->get_dimensions();

      tmp_out = new COMPLEX_ARRAY_TYPE();
      if( !tmp_out->create(osdims.get()) ){
	std::cout << std::endl << "convolutionOperator::mult_MH failed : memory allocation failed (1)" << std::endl;
	return -1;
      } 

      if( !operator_expand( in, tmp_out ) ){
	std::cout << std::endl << "convolutionOperator::mult_MH failed : cannot expand input image" << std::endl;
	return -1;
      }     
    }
    else if( accumulate ){
      tmp_out = new COMPLEX_ARRAY_TYPE(*in);
      if( !tmp_out->get_data_ptr() ){
	std::cout << std::endl << "convolutionOperator::mult_MH failed : memory allocation failed (2)" << std::endl;
	return -1;
      }
    }
    else{ 
      *out = *in;
      tmp_out = out;
    } 

    // Forwards fft
    operator_fft( true, tmp_out );

    // Multiply
    operator_scale( adjoint_kernel_.get(), tmp_out );

    // Inverse fft
    operator_fft( false, tmp_out );

    if( use_oversampling ) {
      if( !operator_crop( tmp_out, out )){
	std::cout << std::endl << "convolutionOperator::mult_MH failed : cannot crop oversampled image" << std::endl;
	return -1;
      }
      delete tmp_out;
    }    
    else if( accumulate ){
      operator_xpy( tmp_out, out );
      delete tmp_out;
    }
    
    return 0;
  }
  
private:
  boost::shared_ptr<COMPLEX_ARRAY_TYPE> kernel_;
  boost::shared_ptr<COMPLEX_ARRAY_TYPE> adjoint_kernel_;
};
