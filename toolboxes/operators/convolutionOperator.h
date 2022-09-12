/** \file convolutionOperator.h
    \brief Base class for all convolution operators.
*/

#pragma once

#include "linearOperator.h"
#include "vector_td_utilities.h"

#include <boost/shared_ptr.hpp>
#include <vector>

namespace Gadgetron{

  template <class COMPLEX_ARRAY_TYPE, unsigned int D> class convolutionOperator : public linearOperator<COMPLEX_ARRAY_TYPE>
  {  
  protected:
    typedef typename COMPLEX_ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;
    
  public:
    
    convolutionOperator() : linearOperator<COMPLEX_ARRAY_TYPE>() {}
    virtual ~convolutionOperator() {}
    
    // Set the convolution kernel
    virtual void set_kernel( COMPLEX_ARRAY_TYPE *image_space_kernel )
    {     
      if (!image_space_kernel) throw std::runtime_error("convolutionOperator: null pointer kernel provided");
      COMPLEX_ARRAY_TYPE *freq_kernel = new COMPLEX_ARRAY_TYPE(*image_space_kernel);
      operator_fft( true, freq_kernel );
      kernel_ = boost::shared_ptr<COMPLEX_ARRAY_TYPE>(freq_kernel);
      
      COMPLEX_ARRAY_TYPE *freq_kernel_adjoint = new COMPLEX_ARRAY_TYPE(freq_kernel->get_dimensions());      
      origin_mirror( freq_kernel, freq_kernel_adjoint );
      adjoint_kernel_ = boost::shared_ptr<COMPLEX_ARRAY_TYPE>(freq_kernel_adjoint);           
    }
    
    // Apply image operators
    //
    
    virtual void mult_MH_M( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out, bool accumulate = false )
    {    
      if( !kernel_.get() ){
	throw std::runtime_error( "convolutionOperator::mult_MH_M failed : kernel is not set");
      }
    
      if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
    	throw std::runtime_error( "convolutionOperator::mult_MH_M failed : in/out image dimensions mismatch");
      }
      
      bool use_oversampling;
      if( in->get_number_of_elements() == kernel_->get_number_of_elements() )
	use_oversampling = false;
      else if( (in->get_number_of_elements()<<D) == kernel_->get_number_of_elements() )
	use_oversampling = true;
      else{
	throw std::runtime_error( "convolutionOperator::mult_MH_M failed : in/out image dimensions mismatch the kernel");
      }
      
      // Intermediate variables
      COMPLEX_ARRAY_TYPE *tmp_out;

      if( use_oversampling ){
	boost::shared_ptr< std::vector<size_t> > osdims = kernel_->get_dimensions();
	tmp_out = new COMPLEX_ARRAY_TYPE(osdims);
	pad<ELEMENT_TYPE,D>( *in, *tmp_out );
      }
      else if( accumulate ){
	tmp_out = new COMPLEX_ARRAY_TYPE(*in);
      }
      else{ 
	*out = *in;
	tmp_out = out;
      } 

      // Forwards fft
      operator_fft( true, tmp_out );

      // Multiply
      *tmp_out *= *kernel_;
      *tmp_out *= *adjoint_kernel_;

      // Inverse fft
      operator_fft( false, tmp_out );

      if( use_oversampling ) {
	operator_crop( tmp_out, out );
	delete tmp_out;
      }    
      else if( accumulate ){
    	*out += *tmp_out;
	delete tmp_out;
      }    
    }
    
  
    virtual void mult_M( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out, bool accumulate = false )
    {
      if( !kernel_.get() ){
    	throw std::runtime_error("convolutionOperator::mult_M failed : kernel is not set");
      }
    
      if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
    	throw std::runtime_error( "convolutionOperator::mult_M failed : in/out image dimensions mismatch");
      }

      bool use_oversampling;
      if( in->get_number_of_elements() == kernel_->get_number_of_elements() )
	use_oversampling = false;
      else if( (in->get_number_of_elements()<<D) == kernel_->get_number_of_elements() )
	use_oversampling = true;
      else{
	throw std::runtime_error( "convolutionOperator::mult_M failed : in/out image dimensions mismatch the kernel");
      }
    
      // Intermediate variables
      COMPLEX_ARRAY_TYPE *tmp_out;

      if( use_oversampling ){
	boost::shared_ptr< std::vector<size_t> > osdims = kernel_->get_dimensions();
	tmp_out = new COMPLEX_ARRAY_TYPE(osdims);
	pad<ELEMENT_TYPE,D>( *in, *tmp_out );
      }
      else if( accumulate ){
	tmp_out = new COMPLEX_ARRAY_TYPE(*in);
      }
      else{ 
	*out = *in;
	tmp_out = out;
      } 

      // Forwards fft
      operator_fft( true, tmp_out );

      // Multiply
      *tmp_out *= *kernel_;
 
      // Inverse fft
      operator_fft( false, tmp_out );

      if( use_oversampling ) {
	operator_crop( tmp_out, out );
	delete tmp_out;
      }    
      else if( accumulate ){
    	*out += *tmp_out;
	delete tmp_out;
      }    
    }
  
    virtual void mult_MH( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out, bool accumulate = false )
    {
      if( !adjoint_kernel_.get() ){
	throw std::runtime_error("convolutionOperator::mult_MH failed : kernel is not set");
      }
      
      if( !in || !out || in->get_number_of_elements() != out->get_number_of_elements() ){
	throw std::runtime_error("convolutionOperator::mult_MH failed : in/out image dimensions mismatch");
      }

      bool use_oversampling;
      if( in->get_number_of_elements() == adjoint_kernel_->get_number_of_elements() )
	use_oversampling = false;
      else if( (in->get_number_of_elements()<<D) == adjoint_kernel_->get_number_of_elements() )
	use_oversampling = true;
      else{
    	throw std::runtime_error( "convolutionOperator::mult_MH failed : in/out image dimensions mismatch the kernel");
      }
      
      // Intermediate variables
      COMPLEX_ARRAY_TYPE *tmp_out;

      if( use_oversampling ){
	boost::shared_ptr< std::vector<size_t> > osdims = adjoint_kernel_->get_dimensions();
	tmp_out = new COMPLEX_ARRAY_TYPE(osdims);
	pad<ELEMENT_TYPE,D>( *in, *tmp_out );
      }
      else if( accumulate ){
	tmp_out = new COMPLEX_ARRAY_TYPE(*in);
      }
      else{ 
	*out = *in;
	tmp_out = out;
      } 
      
      // Forwards fft
      operator_fft( true, tmp_out );

      // Multiply
      *tmp_out *= *adjoint_kernel_;

      // Inverse fft
      operator_fft( false, tmp_out );

      if( use_oversampling ) {
	operator_crop( tmp_out, out );
	delete tmp_out;
      }    
      else if( accumulate ){
    	*out += *tmp_out;
	delete tmp_out;
      }
    }

  protected:
  
    virtual void operator_fft( bool forwards_transform, COMPLEX_ARRAY_TYPE *image ) = 0;    
    virtual void origin_mirror( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out ) = 0;

    virtual void operator_crop( COMPLEX_ARRAY_TYPE *in, COMPLEX_ARRAY_TYPE *out ){
      typename uint64d<D>::Type offset = from_std_vector<size_t,D>(*(in->get_dimensions().get()))>>2;
       typename uint64d<D>::Type size = from_std_vector<size_t,D>(*(out->get_dimensions().get()))>>2;
      crop<ELEMENT_TYPE,D>( offset,size, *in, *out );
    }
    
  private:
    boost::shared_ptr<COMPLEX_ARRAY_TYPE> kernel_;
    boost::shared_ptr<COMPLEX_ARRAY_TYPE> adjoint_kernel_;
  };
}
