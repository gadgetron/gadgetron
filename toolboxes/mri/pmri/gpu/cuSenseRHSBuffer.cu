#include "cuSenseRHSBuffer.h"
#include "vector_td_utilities.h"
#include "NFFT_utils.h"

namespace Gadgetron{

  template<class REAL, unsigned int D, bool ATOMICS>
  void cuSenseRHSBuffer<REAL,D,ATOMICS>::clear()
  {
    if( sense_op_.get() == 0x0 ){
      throw std::runtime_error("cuSenseRHSBuffer::clear: nothing to clear; sense operator not yet set." );
    }

    Gadgetron::clear(&acc_buffer_);
    Gadgetron::clear(&cyc_buffer_);

    cur_idx_ = cur_sub_idx_ = 0;
    acc_buffer_empty_ = true;
  }

  template<class REAL, unsigned int D, bool ATOMICS>
  void cuSenseRHSBuffer<REAL,D,ATOMICS>::set_sense_operator( boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D,ATOMICS> > op )
  {
    if( op.get() == 0x0 ){
      throw std::runtime_error("cuSenseRHSBuffer::set_sense_operator: invalid array pointer");
    }

    if( !op->is_setup() ){
      throw std::runtime_error("cuSenseRHSBuffer::set_sense_operator: operator has not been set up");
    }

    if( num_coils_ == 0 ){
      throw std::runtime_error("cuSenseRHSBuffer::set_sense_operator: number of coils not set");
    }

    sense_op_ = op;

    _uintd image_dims = op->get_plan()->get_matrix_size_os();

    std::vector<unsigned int> dims = to_std_vector( image_dims );
    dims.push_back(num_coils_);
    acc_buffer_.create(&dims);

    dims.push_back(cycle_length_); // Number of frames in the cyclic buffer
    cyc_buffer_.create(&dims);
    clear();
  }

  template<class REAL, unsigned int D, bool ATOMICS> 
  void cuSenseRHSBuffer<REAL,D,ATOMICS>::add_frame_data( cuNDArray<_complext> *samples, cuNDArray<_reald> *trajectory )
  {
    if( sense_op_.get() == 0x0 ){
      throw std::runtime_error("cuSenseRHSBuffer::add_frame_data: sense_operator not set");
    }

    if( !samples || !trajectory ){
      throw std::runtime_error("cuSenseRHSBuffer::add_frame_data: illegal input pointer");
    }

    // Make array containing the "current" buffer from the cyclic buffer
    //
    cuNDArray<_complext> cur_buffer(acc_buffer_.get_dimensions().get(),
				    cyc_buffer_.get_data_ptr()+cur_idx_*acc_buffer_.get_number_of_elements());

    // Preprocess frame
    //
    sense_op_->get_plan()->preprocess( trajectory, cuNFFT_plan<REAL,D,ATOMICS>::NFFT_PREP_NC2C );

    // Convolve to form k-space frame (accumulation mode)
    //
    sense_op_->get_plan()->convolve( samples, &cur_buffer, sense_op_->get_dcw().get(), cuNFFT_plan<REAL,D,ATOMICS>::NFFT_CONV_NC2C, true );

    // Update the accumulation buffer (if it is time...)
    //
    if( cur_sub_idx_ == sub_cycle_length_-1 ){

      // Buffer complete, add to accumulation buffer
      //
      acc_buffer_ += cur_buffer;
      acc_buffer_empty_ = false;

      // Start filling the next buffer in the cycle ...
      //
      cur_idx_++; 
      if( cur_idx_ == cycle_length_ ) cur_idx_ = 0;

      // ... but first subtract this next buffer from the accumulation buffer
      //
      cur_buffer.create( acc_buffer_.get_dimensions().get(), cyc_buffer_.get_data_ptr()+cur_idx_*acc_buffer_.get_number_of_elements() );
      acc_buffer_ -= cur_buffer;

      // Clear new buffer before refilling
      //
      Gadgetron::clear(&cur_buffer);
    }

    cur_sub_idx_++;
    if( cur_sub_idx_ == sub_cycle_length_ ) cur_sub_idx_ = 0;
  }

  template<class REAL, unsigned int D, bool ATOMICS>
  boost::shared_ptr< cuNDArray<complext<REAL> > > cuSenseRHSBuffer<REAL,D,ATOMICS>::get_acc_coil_images( bool normalize )
  {
    if( sense_op_.get() == 0x0 ){
      throw std::runtime_error("cuSenseRHSBuffer::get_acc_coil_images: sense_operator not set");
    }

    // Prepare return image
    _uintd image_dims = sense_op_->get_plan()->get_matrix_size();
    _uintd image_dims_os = sense_op_->get_plan()->get_matrix_size_os();

    std::vector<unsigned int> dims = to_std_vector( image_dims );
    dims.push_back(num_coils_);

    cuNDArray<_complext> *image = new cuNDArray<_complext>(&dims);
    // Check if we are ready to reconstruct. If not return an image of ones...
    if( acc_buffer_empty_ ){
      fill(image,_complext(1));
      return boost::shared_ptr< cuNDArray<_complext> >(image);
    }

    // Complete gridding of k-space CSM image
    //

    // Copy accumulation buffer before in-place FFT
    cuNDArray<_complext> acc_copy = acc_buffer_;

    // FFT
    sense_op_->get_plan()->fft( &acc_copy, cuNFFT_plan<REAL,D,ATOMICS>::NFFT_BACKWARDS );

    // Deapodize
    sense_op_->get_plan()->deapodize( &acc_copy );

    // Remove oversampling
    crop<_complext,D>( (image_dims_os-image_dims)>>1, &acc_copy, image );

    if( normalize ){
      REAL scale = REAL(1)/(((REAL)cycle_length_-REAL(1))*(REAL)sub_cycle_length_);
      *image *= scale;
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


  template class EXPORTGPUPMRI cuSenseRHSBuffer<double,2,false>;
  template class EXPORTGPUPMRI cuSenseRHSBuffer<double,3,false>;
  template class EXPORTGPUPMRI cuSenseRHSBuffer<double,4,false>;
}
