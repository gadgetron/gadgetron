#include "GrappaCalibrationBuffer.h"

#include "Gadgetron.h"

GrappaCalibrationBuffer::GrappaCalibrationBuffer(std::vector<unsigned int> dimensions,
						 GrappaWeights<float>* w,
						 GrappaWeightsCalculator<float>* weights_calculator)
  : weights_(w)
  , weights_calculator_(weights_calculator)
  , buffer_counter_(dimensions[1])
  , biggest_gap_current_(0)
  , acceleration_factor_(0)
  , last_line_(0)
  , weights_invalid_(true)
{
  dimensions_ = dimensions;
  if (!buffer_.create(&dimensions_)) {
    GADGET_DEBUG1("Unable to allocate memory for GRAPPA buffer");
  }
  
}

int GrappaCalibrationBuffer::add_data(GadgetMessageAcquisition* m1, hoNDArray< std::complex<float> >* m2)
{
  if (!buffer_.get_data_ptr()) {
    GADGET_DEBUG1("Buffer not allocated, cannot add data");
    return GADGET_FAIL;
  }
  
  unsigned int samples =  m1->samples;
  unsigned int line = m1->idx.line;
  unsigned int partition = m1->idx.partition;
  //unsigned int slice = m1->idx.slice; //We should probably check this

  if (samples != dimensions_[0]) {
    GADGET_DEBUG1("Wrong number of samples received\n");
    return GADGET_FAIL;    
  }

  std::complex<float>* b = buffer_.get_data_ptr();
  std::complex<float>* d = m2->get_data_ptr();

  size_t offset= 0;
  //Copy the data for all the channels
  for (int c = 0; c < m1->channels; c++) {
    offset = 
      c*dimensions_[0]*dimensions_[1]*dimensions_[2] +
      partition*dimensions_[0]*dimensions_[1] +
      line*dimensions_[0];
    memcpy(b+offset,d+c*samples,sizeof(std::complex<float>)*samples);
  }

  int buf_update  = buffer_counter_.update_line(line,m1->position,m1->quarternion);
  if ( buf_update < 0) {
    GADGET_DEBUG2("Unable to update buffer counter for line %d\n", line);
    return GADGET_FAIL;
  }

  //Let's figure out if we should start a weight calculation job
  
  //This means that the orientation changed
  if (buf_update == 1) {
    weights_invalid_ = true;
  }

  bool is_first_scan_in_slice =
    (m1->flags & GADGET_FLAG_FIRST_ACQ_IN_SLICE);

  if (is_first_scan_in_slice) {
    biggest_gap_current_ = 0;
  } else {
    unsigned int gap = abs(static_cast<int>(last_line_) - static_cast<int>(line));
    if (gap > biggest_gap_current_) biggest_gap_current_ = gap;
  }
  last_line_ = line;

  bool is_last_scan_in_slice =
    (m1->flags & GADGET_FLAG_LAST_ACQ_IN_SLICE);

  
  if (is_last_scan_in_slice) { 
    unsigned int min_ky, max_ky;

    if (biggest_gap_current_ != acceleration_factor_) {
      acceleration_factor_ = biggest_gap_current_;
      weights_invalid_ = true;
    }
 
    if (buffer_counter_.get_region_of_support(min_ky, max_ky) < 0) {
      GADGET_DEBUG1("Unable to query min_ky, max_ky\n");
      return GADGET_FAIL;
    }
    
    if (weights_invalid_ && ((max_ky-min_ky) > acceleration_factor_)) {
      std::vector< std::pair<unsigned int, unsigned int> > sampled_region;
      sampled_region.push_back(std::pair<unsigned int, unsigned int>(0, samples-1));
      sampled_region.push_back(std::pair<unsigned int, unsigned int>(min_ky, max_ky));

      std::vector<unsigned int> uncombined_channel_weights;

      if (!weights_calculator_) {
	GADGET_DEBUG1("Weights calculator not defined\n");
	return GADGET_FAIL;
      }

      weights_calculator_->add_job( &buffer_,
				   sampled_region,
				   acceleration_factor_,
				   weights_,
				   uncombined_channel_weights,
				   true);

      weights_invalid_ = false;
    }
  }


  return GADGET_OK;
}
