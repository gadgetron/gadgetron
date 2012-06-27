#include "cuSenseRHSBuffer.h"
#include "vector_td_utilities.h"
#include "ndarray_vector_td_utilities.h"

template<class REAL, unsigned int D, bool ATOMICS>
 int cuSenseRHSBuffer<REAL,D,ATOMICS>::clear()
{
  if( sense_op_.get() == 0x0 ){
    std::cerr << "cuSenseRHSBuffer::clear: nothing to clear; sense operator not yet set." << std::endl;
    return -1;
  }

  if( !cuNDA_clear(&acc_buffer_, _complext(0)) || !cuNDA_clear(&cyc_buffer_, _complext(0) )){
    std::cerr << "cuSenseRHSBuffer::clear: failed to clear buffers" << std::endl;
    return -1;
  }

  cur_idx_ = cur_sub_idx_ = 0;
  acc_buffer_empty_ = true;

  return 0;
}

template<class REAL, unsigned int D, bool ATOMICS>
int cuSenseRHSBuffer<REAL,D,ATOMICS>::set_sense_operator( boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D,ATOMICS> > op )
{
  if( op.get() == 0x0 ){
    std::cerr << "cuSenseRHSBuffer::set_sense_operator: invalid array pointer" << std::endl;
    return -1;
  }

  if( !op->is_setup() ){
    std::cerr << "cuSenseRHSBuffer::set_sense_operator: operator has not been set up" << std::endl;
    return -1;
  }

  if( num_coils_ == 0 ){
    std::cerr << "cuSenseRHSBuffer::set_sense_operator: number of coils not set" << std::endl;
    return -1;
  }

  sense_op_ = op;

  _uintd image_dims = op->get_plan()->get_matrix_size_os();

  std::vector<unsigned int> dims = uintd_to_vector<D>( image_dims );
  dims.push_back(num_coils_);

  if( acc_buffer_.create(&dims) < 0 ){
    std::cerr << "cuSenseRHSBuffer::set_sense_operator: failed to create accumulation buffer" << std::endl;
    return -1;
  }

  dims.push_back(cycle_length_); // Number of frames in the cyclic buffer

  if( cyc_buffer_.create(&dims) < 0 ){
    std::cerr << "cuSenseRHSBuffer::set_sense_operator: failed to create cyclic buffer" << std::endl;
    return -1;
  }

  if( clear() < 0 ){
    std::cerr << "cuSenseRHSBuffer::set_sense_operator: failed to clear buffers" << std::endl;
    return -1;
  }

  return 0;
}

template<class REAL, unsigned int D, bool ATOMICS> 
int cuSenseRHSBuffer<REAL,D,ATOMICS>::add_frame_data( cuNDArray<_complext> *samples, cuNDArray<_reald> *trajectory )
{
  if( sense_op_.get() == 0x0 ){
    std::cerr << "cuSenseRHSBuffer::add_frame_data: sense_operator not set" << std::endl;
    return -1;
  }

  if( !samples || !trajectory ){
    std::cerr << "cuSenseRHSBuffer::add_frame_data: illegal input pointer" << std::endl;
    return -1;
  }

  // Make array containing the "current" buffer from the cyclic buffer
  //
  cuNDArray<_complext> cur_buffer;
  if( cur_buffer.create( acc_buffer_.get_dimensions().get(), 
			 cyc_buffer_.get_data_ptr()+cur_idx_*acc_buffer_.get_number_of_elements() ) < 0 ){
    std::cerr << "cuSenseRHSBuffer::add_frame_data: failed creating sub-array" << std::endl;
    return -1;
  }

  // Preprocess frame
  //
  if( sense_op_->get_plan()->preprocess( trajectory, NFFT_plan<REAL,D,ATOMICS>::NFFT_PREP_NC2C ) < 0 ){
    std::cerr << "cuSenseRHSBuffer::add_frame_data: NFFT preprocessing failed" << std::endl;
    return -1;	
  }

  // Convolve to form k-space frame (accumulation mode)
  //
  if( !sense_op_->get_plan()->convolve( samples, &cur_buffer, sense_op_->get_dcw().get(), NFFT_plan<REAL,D,ATOMICS>::NFFT_CONV_NC2C, true ) ){
    std::cerr << "cuSenseRHSBuffer::add_frame_data: NFFT convolution failed" << std::endl;
    return -1;	
  }

  // Update the accumulation buffer (if it is time...)
  //
  if( cur_sub_idx_ == sub_cycle_length_-1 ){

    // Buffer complete, add to accumulation buffer
    //
    if( !cuNDA_axpy( _complext(1), &cur_buffer, &acc_buffer_ )){
      std::cerr << "cuSenseRHSBuffer::add_frame_data: accumulation failed" << std::endl;
      return -1;
    }

    acc_buffer_empty_ = false;

    // Start filling the next buffer in the cycle ...
    //
    cur_idx_++; 
    if( cur_idx_ == cycle_length_ ) cur_idx_ = 0;

    // ... but first subtract this next buffer from the accumulation buffer
    //
    if( cur_buffer.create( acc_buffer_.get_dimensions().get(), cyc_buffer_.get_data_ptr()+cur_idx_*acc_buffer_.get_number_of_elements() ) < 0 ){
      std::cerr << "cuSenseRHSBuffer::add_frame_data: failed creating sub-array" << std::endl;
      return -1;
    }

    if( !cuNDA_axpy( _complext(0)-_complext(1), &cur_buffer, &acc_buffer_ )){
      std::cerr << "cuSenseRHSBuffer::add_frame_data: failed subtracting buffers" << std::endl;
      return -1;
    }

    // Clear new buffer before refilling
    //
    if( !cuNDA_clear( &cur_buffer ) ){
      std::cerr << "cuSenseRHSBuffer::add_frame_data: failed to clear current buffer" << std::endl;
      return -1;
    }
  }

  cur_sub_idx_++;
  if( cur_sub_idx_ == sub_cycle_length_ ) cur_sub_idx_ = 0;

  return 0;
}

template<class REAL, unsigned int D, bool ATOMICS>
boost::shared_ptr< cuNDArray<complext<REAL> > > cuSenseRHSBuffer<REAL,D,ATOMICS>::get_acc_coil_images( bool normalize )
{
  if( sense_op_.get() == 0x0 ){
    std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: sense_operator not set" << std::endl;
    return boost::shared_ptr< cuNDArray<complext<REAL> > >();
  }

  // Prepare return image
  _uintd image_dims = sense_op_->get_plan()->get_matrix_size();
  _uintd image_dims_os = sense_op_->get_plan()->get_matrix_size_os();

  std::vector<unsigned int> dims = uintd_to_vector<D>( image_dims );
  dims.push_back(num_coils_);

  cuNDArray<_complext> *image = new cuNDArray<_complext>(); 
  if( image->create(&dims) == 0x0 ){
    std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: device memory allocation failed" << std::endl;
    return boost::shared_ptr< cuNDArray<_complext> >();
  }

  // Check if we are ready to reconstruct. If not return an image of ones...
  if( acc_buffer_empty_ ){
    cuNDA_clear(image, _complext(1));
    return boost::shared_ptr< cuNDArray<_complext> >(image);
  }

  // Complete gridding of k-space CSM image
  //

  // Copy accumulation buffer before in-place FFT
  cuNDArray<_complext> acc_copy = acc_buffer_;

  // FFT
  if( !sense_op_->get_plan()->fft( &acc_copy, NFFT_plan<REAL,D,ATOMICS>::NFFT_BACKWARDS )){
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

  if( normalize ){
    REAL scale = REAL(1)/(((REAL)cycle_length_-REAL(1))*(REAL)sub_cycle_length_);
    if( !cuNDA_scal( scale, image ) ){
      std::cerr << "cuSenseRHSBuffer::get_acc_coil_images: normalization failed" << std::endl;
      return boost::shared_ptr< cuNDArray<_complext> >();
    }
  }

  return boost::shared_ptr< cuNDArray<_complext> >(image);
}


//
// Instantiations
//

template class EXPORTGPUPMRI cuSenseRHSBuffer<float,2,true>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<float,2,false>;

template class EXPORTGPUPMRI cuSenseRHSBuffer<float,3,true>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<float,3,false>;

template class EXPORTGPUPMRI cuSenseRHSBuffer<float,4,true>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<float,4,false>;

template class EXPORTGPUPMRI cuSenseRHSBuffer<double,2,true>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<double,2,false>;

template class EXPORTGPUPMRI cuSenseRHSBuffer<double,3,true>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<double,3,false>;

template class EXPORTGPUPMRI cuSenseRHSBuffer<double,4,true>;
template class EXPORTGPUPMRI cuSenseRHSBuffer<double,4,false>;
