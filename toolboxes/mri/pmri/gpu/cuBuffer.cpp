#include "cuBuffer.h"
#include "vector_td_utilities.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_utils.h"

namespace Gadgetron{

  template<class REAL, unsigned int D>
  cuBuffer<REAL,D>::cuBuffer()
  {
    acc_buffer_ = boost::shared_ptr< cuNDArray<_complext> >(new cuNDArray<_complext>);
    cyc_buffer_ = boost::shared_ptr< cuNDArray<_complext> >(new cuNDArray<_complext>);
    num_coils_ = 0;
    cur_idx_ = cur_sub_idx_ = 0;
    cycle_length_ = 0; sub_cycle_length_ = 0;
    acc_buffer_empty_ = true;
    Gadgetron::clear(matrix_size_);
    Gadgetron::clear(matrix_size_os_);
    W_ = REAL(0);
  }
  
  template<class REAL, unsigned int D>
  void cuBuffer<REAL,D>::clear()
  {
    Gadgetron::clear(acc_buffer_.get());
    Gadgetron::clear(cyc_buffer_.get());
    cur_idx_ = cur_sub_idx_ = 0;
    acc_buffer_empty_ = true;
  }

  template<class REAL, unsigned int D>
  void cuBuffer<REAL,D>
  ::setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W, 
           unsigned int num_coils, unsigned int num_cycles, unsigned int num_sub_cycles )
  {      
    bool matrix_size_changed = (matrix_size_ == matrix_size);
    bool matrix_size_os_changed = (matrix_size_os_ == matrix_size_os);
    bool kernel_changed = (W_ == W);
    bool num_coils_changed = (num_coils_ == num_coils );
    bool num_cycles_changed = (cycle_length_ == num_cycles+1);

    matrix_size_ = matrix_size;
    matrix_size_os_ = matrix_size_os;
    W_ = W;
    num_coils_ = num_coils;
    cycle_length_ = num_cycles+1; // +1 as we need a "working buffer" in a addition to 'cycle_length' full ones
    sub_cycle_length_ = num_sub_cycles;

    if( !nfft_plan_ || matrix_size_changed || matrix_size_os_changed || kernel_changed ){
      nfft_plan_ = NFFT<cuNDArray,REAL,D>::make_plan( matrix_size_, matrix_size_os_, W );
    }
    
    std::vector<size_t> dims = to_std_vector(matrix_size_os_);    
    dims.push_back(num_coils_);

    if( acc_buffer_->get_number_of_elements() == 0 || matrix_size_os_changed || num_coils_changed ){
      acc_buffer_->create(dims);
      Gadgetron::clear( acc_buffer_.get() );
    }

    dims.push_back(cycle_length_);
    if( cyc_buffer_->get_number_of_elements() == 0 || matrix_size_os_changed || num_coils_changed ){
      cyc_buffer_->create(dims);      
      Gadgetron::clear( cyc_buffer_.get() );
    }
    else if( num_cycles_changed ){
      // Reuse the old buffer content in this case...
      // This happens automatically (in all cases?) with the current design?
    }
  }
  
  template<class REAL, unsigned int D>
  bool cuBuffer<REAL,D>::add_frame_data( cuNDArray<_complext> *samples, cuNDArray<_reald> *trajectory )
  {
    if( !samples || !trajectory ){
      throw std::runtime_error("cuBuffer::add_frame_data: illegal input pointer");
    }

    if( num_coils_ != samples->get_size(samples->get_number_of_dimensions()-1) ){
      throw std::runtime_error("cuBuffer::add_frame_data: unexpected number of coils according to setup");
    }

    //if( dcw_.get() == 0x0 ){
    //throw std::runtime_error("cuBuffer::density compensation weights not set");
    //}
    
    // Make array containing the "current" buffer from the cyclic buffer
    //

    cuNDArray<_complext> cur_buffer(*acc_buffer_->get_dimensions(),
				    cyc_buffer_->get_data_ptr()+cur_idx_*acc_buffer_->get_number_of_elements());

    // Preprocess frame
    //

    nfft_plan_->preprocess( *trajectory, NFFT_prep_mode::NC2C );
    
    // Convolve to form k-space frame (accumulation mode)
    //

    {
      if (dcw_) {
        auto samples_rescaled = *samples;
        samples_rescaled *= *dcw_;
        nfft_plan_->convolve(samples_rescaled, cur_buffer, NFFT_conv_mode::NC2C, true);
      } else {
        nfft_plan_->convolve(*samples, cur_buffer, NFFT_conv_mode::NC2C, true);
      }

    }
    // Update the accumulation buffer (if it is time...)
    //

    bool cycle_completed = false;

    if( cur_sub_idx_ == sub_cycle_length_-1 ){

      cycle_completed = true;
      
      // Buffer complete, add to accumulation buffer
      //

      *acc_buffer_ += cur_buffer;
      acc_buffer_empty_ = false;

      // Start filling the next buffer in the cycle ...
      //

      cur_idx_++; 
      if( cur_idx_ == cycle_length_ ) cur_idx_ = 0;

      // ... but first subtract this next buffer from the accumulation buffer
      //

      cur_buffer.create( *acc_buffer_->get_dimensions(), cyc_buffer_->get_data_ptr()+cur_idx_*acc_buffer_->get_number_of_elements() );
      *acc_buffer_ -= cur_buffer;

      // Clear new buffer before refilling
      //

      Gadgetron::clear(&cur_buffer);
    }

    cur_sub_idx_++;
    if( cur_sub_idx_ == sub_cycle_length_ ) cur_sub_idx_ = 0;

    return cycle_completed;
  }

  template<class REAL, unsigned int D>
  boost::shared_ptr< cuNDArray<complext<REAL> > > cuBuffer<REAL,D>::get_accumulated_coil_images()
  {
    std::vector<size_t> dims = to_std_vector(matrix_size_);
    dims.push_back(num_coils_);

    acc_image_ = boost::shared_ptr< cuNDArray<_complext> >( new cuNDArray<_complext>(dims) );
				    
    // Check if we are ready to reconstruct. If not return an image of ones...
    if( acc_buffer_empty_ ){
      fill(acc_image_.get(),_complext(1));
      return acc_image_;
    }

    // Finalize gridding of k-space CSM image (convolution has been done already)
    //

    // Copy accumulation buffer before in-place FFT
    cuNDArray<_complext> acc_copy = *acc_buffer_;

    // FFT
    nfft_plan_->fft( acc_copy, NFFT_fft_mode::BACKWARDS );
    
    // Deapodize
    nfft_plan_->deapodize( acc_copy );
    
    // Remove oversampling
    crop<_complext,D>( (matrix_size_os_-matrix_size_)>>1, matrix_size_, acc_copy, *acc_image_ );
    
    //if( normalize ){
    //REAL scale = REAL(1)/(((REAL)cycle_length_-REAL(1))*(REAL)sub_cycle_length_);
    //*acc_image_ *= scale;
    //}
    
    return acc_image_;
  }

  //
  // Instantiations
  //
  
  template class EXPORTGPUPMRI cuBuffer<float,2>;
  template class EXPORTGPUPMRI cuBuffer<float,3>;
  template class EXPORTGPUPMRI cuBuffer<float,4>;

  template class EXPORTGPUPMRI cuBuffer<double,2>;
  template class EXPORTGPUPMRI cuBuffer<double,3>;
  template class EXPORTGPUPMRI cuBuffer<double,4>;
}
