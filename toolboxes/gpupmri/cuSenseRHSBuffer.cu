#include "cuSenseRHSBuffer.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template<class REAL, unsigned int D> 
int cuSenseRHSBuffer<REAL,D>::set_sense_operator( boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D> > op )
{
  if( op.get() == 0x0 ){
	std::cerr << "cuSenseRHSBuffer::set_sense_operator: invalid array pointer" << std::endl;
	return -1;
  }

  if( !op->is_setup() ){
	std::cerr << "cuSenseRHSBuffer::set_sense_operator: operator has not been set up" << std::endl;
	return -1;
  }

  sense_op_ = op;

  _uintd image_dims = op->get_plan()->get_matrix_size_os();
  unsigned int num_coils = op->get_number_of_coils();

  std::vector<unsigned int> dims = uintd_to_vector<D>( image_dims );
  dims.push_back(num_coils);

  if( acc_buffer_.create(&dims) < 0 ){
	std::cerr << "cuSenseRHSBuffer::set_sense_operator: failed to create accumulation buffer" << std::endl;
	return -1;
  }
	
  dims.push_back(cycle_length_); // Number of frames in the cyclic buffer

  if( cyc_buffer_.create(&dims) < 0 ){
	std::cerr << "cuSenseRHSBuffer::set_sense_operator: failed to create cyclic buffer" << std::endl;
	return -1;
  }

  if( !cuNDA_clear(&acc_buffer_, get_zero<_complext>()) || !cuNDA_clear(&cyc_buffer_, get_zero<_complext>() )){
	std::cerr << "cuSenseRHSBuffer::set_sense_operator: failed to clear buffers" << std::endl;
	return -1;
  }

  return 0;
}

template<class REAL, unsigned int D> int 
cuSenseRHSBuffer<REAL,D>::add_frame_data( cuNDArray<_complext> *samples, cuNDArray<_reald> *trajectory )
{
	if( !samples ){
		std::cerr << "cuSenseRHSBuffer::add_frame_data: illegal input pointer" << std::endl;
		return -1;
	}

	// Make array containing the k-th buffer from the cyclic buffer
	cuNDArray<_complext> cur_buffer;
	if( cur_buffer.create( acc_buffer_.get_dimensions().get(), cyc_buffer_.get_data_ptr()+cur_idx_*acc_buffer_.get_number_of_elements() ) < 0 ){
		std::cerr << "cuSenseRHSBuffer::add_frame_data: failed creating sub-array" << std::endl;
		return -1;
	}

	// Subtract the k-th buffer from the accumulation buffer and clear current buffer
	if( cur_sub_idx_ == 0 ){
  	  if( !cuNDA_axpy( get_zero<_complext>()-get_one<_complext>(), &cur_buffer, &acc_buffer_ )){
		  std::cerr << "cuSenseRHSBuffer::add_frame_data: failed subtracting buffers" << std::endl;
		  return -1;
	  }
	  if( !cuNDA_clear( &cur_buffer ) ){
		  std::cerr << "cuSenseRHSBuffer::add_frame_data: failed to clear current buffer" << std::endl;
		  return -1;
	  }
	}

	// Preprocess frame
	if( sense_op_->get_plan()->preprocess( trajectory, NFFT_plan<REAL,D>::NFFT_PREP_BACKWARDS ) < 0 ){
		std::cerr << "cuSenseRHSBuffer::add_frame_data: NFFT preprocessing failed" << std::endl;
		return -1;	
	}

	// Convolve to acquire k-space frame
	if( !sense_op_->get_plan()->convolve( samples, &cur_buffer, sense_op_->get_dcw().get(), NFFT_plan<REAL,D>::NFFT_BACKWARDS, true ) ){
		std::cerr << "cuSenseRHSBuffer::add_frame_data: NFFT convolution failed" << std::endl;
		return -1;	
	}

	// Accumulate k-space for CSM/regularization estimation
	if( cur_sub_idx_ == sub_cycle_length_-1 ){
		if( !cuNDA_axpy( get_one<_complext>(), &cur_buffer, &acc_buffer_ )){
		  std::cerr << "cuSenseRHSBuffer::add_frame_data: accumulation failed" << std::endl;
		  return -1;
	  }
	  acc_buffer_empty_ = false;
  	  cur_idx_++; 
	  if( cur_idx_ == cycle_length_ ) cur_idx_ = 0;
	}

	cur_sub_idx_++;
	if( cur_sub_idx_ == sub_cycle_length_ ) cur_sub_idx_ = 0;
	
	return 0;
}

template<class REAL, unsigned int D>
boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> > cuSenseRHSBuffer<REAL,D>::get_acc_coil_images()
{
  // Prepare return image
  _uintd image_dims = sense_op_->get_plan()->get_matrix_size();
  _uintd image_dims_os = sense_op_->get_plan()->get_matrix_size_os();
  unsigned int num_coils = sense_op_->get_number_of_coils();

  std::vector<unsigned int> dims = uintd_to_vector<D>( image_dims );
  dims.push_back(num_coils);

  cuNDArray<_complext> *image = new cuNDArray<_complext>(); 
  if( image->create(&dims) == 0x0 ){
	  std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: device memory allocation failed" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();
  }

  // Check if we are ready to reconstruct. If not return an image of ones...
  if( acc_buffer_empty_ ){
	  cuNDA_clear(image, get_one<_complext>());
	  return boost::shared_ptr< cuNDArray<_complext> >(image);
  }

  // Complete gridding of k-space CSM image
  //

  // Copy accumulation buffer before in-place FFT
  cuNDArray<_complext> acc_copy = acc_buffer_;

  // FFT
  if( !sense_op_->get_plan()->fft( &acc_copy, NFFT_plan<REAL,D>::NFFT_BACKWARDS )){
	  std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: fft failed" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();
  }

  // Deapodize
  if( !sense_op_->get_plan()->deapodize( &acc_copy )){
	  std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: deapoization failed" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();
  }

  // Remove oversampling
  if( !cuNDA_crop<_complext,D>( (image_dims_os-image_dims)>>1, &acc_copy, image )){
	  std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: image crop failed" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();
  }

  return boost::shared_ptr< cuNDArray<_complext> >(image);
}

template<class REAL, unsigned int D>
boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> > cuSenseRHSBuffer<REAL,D>::get_cur_rhs()
{
  if( sense_op_->get_csm().get() == 0x0 ){
    std::cerr << "cuSenseRHSBuffer::get_cur_rhs: csm not set" << std::endl;
	return boost::shared_ptr< cuNDArray<_complext> >();
  }

  // Copy current buffer before in-place FFT
  cuNDArray<_complext> cur_buffer;
  if( cur_buffer.create( acc_buffer_.get_dimensions().get(), cyc_buffer_.get_data_ptr()+((cur_idx_+cycle_length_-1)%cycle_length_)*acc_buffer_.get_number_of_elements() ) < 0 ){
	  std::cerr << "cuSenseRHSBuffer::add_frame_data: failed creating sub-array" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();;
  }
  cuNDArray<_complext> cur_copy = cur_buffer;

  // FFT
  if( !sense_op_->get_plan()->fft( &cur_copy, NFFT_plan<REAL,D>::NFFT_BACKWARDS )){
	  std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: fft failed" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();
  }

  // Deapodize
  if( !sense_op_->get_plan()->deapodize( &cur_copy )){
	  std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: deapoization failed" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();
  }

  // Remove oversampling
  _uintd image_dims = sense_op_->get_plan()->get_matrix_size();
  _uintd image_dims_os = sense_op_->get_plan()->get_matrix_size_os();
  unsigned int num_coils = sense_op_->get_number_of_coils();

  std::vector<unsigned int> dims = uintd_to_vector<D>( image_dims );
  dims.push_back(num_coils);

  cuNDArray<_complext> image;
  if( image.create(&dims) == 0x0 ){
	  std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: device memory allocation failed" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();
  }

  if( !cuNDA_crop<_complext,D>( (image_dims_os-image_dims)>>1, &cur_copy, &image )){
	  std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: image crop failed" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();
  }

  dims.pop_back();
  cuNDArray<_complext> *combined_image = new cuNDArray<_complext>();
  if( image.create(&dims) == 0x0 ){
	  std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: device memory allocation failed" << std::endl;
	  return boost::shared_ptr< cuNDArray<_complext> >();
  }

  if( sense_op_->mult_csm_conj_sum( &image, combined_image ) < 0 ) {
    std::cerr << "cuSenseRHSBuffer::get_cur_rhs: Unable to multiply with conjugate of sensitivity maps and sum" << std::endl;
    return boost::shared_ptr< cuNDArray<_complext> >(); 
  }
  
  return boost::shared_ptr< cuNDArray<_complext> >(combined_image); 
}

//
// Instantiations
//

template class EXPORTGPUPMRI cuSenseRHSBuffer<float,2>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<float,3>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<float,4>;

template class EXPORTGPUPMRI cuSenseRHSBuffer<double,2>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<double,3>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<double,4>;
