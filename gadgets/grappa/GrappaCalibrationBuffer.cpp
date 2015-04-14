#include "GrappaCalibrationBuffer.h"
#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron{

  GrappaCalibrationBuffer::GrappaCalibrationBuffer(std::vector<size_t> dimensions,
                                                   boost::shared_ptr<GrappaWeights<float> > w,
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
    try {
      buffer_.create(&dimensions_);
      buffer_.fill(std::complex<float>(0.0,0.0));
    } catch (std::runtime_error & err){
      GEXCEPTION(err,"Unable to allocate memory for GRAPPA buffer");
    }
  
  }

  int GrappaCalibrationBuffer::add_data(ISMRMRD::AcquisitionHeader* m1, hoNDArray< std::complex<float> >* m2,
					 unsigned short line_offset, unsigned short partition_offset)
  {
    if (!buffer_.get_data_ptr()) {
      GDEBUG("Buffer not allocated, cannot add data");
      return GADGET_FAIL;
    }
  
    unsigned int samples =  m1->number_of_samples;
    unsigned int line = m1->idx.kspace_encode_step_1 + line_offset;
    unsigned int partition = m1->idx.kspace_encode_step_2 + partition_offset;
    unsigned int slice = m1->idx.slice; //We should probably check this

    if (samples != dimensions_[0]) {
      GDEBUG("Wrong number of samples received\n");
      return GADGET_FAIL;    
    }

    std::complex<float>* b = buffer_.get_data_ptr();
    std::complex<float>* d = m2->get_data_ptr();

    size_t offset= 0;
    //Copy the data for all the channels
    for (int c = 0; c < m1->active_channels; c++) {
      offset = 
        c*dimensions_[0]*dimensions_[1]*dimensions_[2] +
        partition*dimensions_[0]*dimensions_[1] +
        line*dimensions_[0];
      memcpy(b+offset,d+c*samples,sizeof(std::complex<float>)*samples);
    }

    int buf_update  = buffer_counter_.update_line(line, m1->position,
                                                  m1->read_dir, m1->phase_dir, m1->slice_dir);

    if ( buf_update < 0) {
      GDEBUG("Unable to update buffer counter for line %d\n", line);
      return GADGET_FAIL;
    }

    //Let's figure out if we should start a weight calculation job
  
    //This means that the orientation changed
    if (buf_update == 1) {
      weights_invalid_ = true;
    }

    bool is_first_scan_in_slice = m1->isFlagSet(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);


    //Depending on the sequence used, we could get into trouble if the sequence switches slice acquisition scheme before finishing a slice.
    bool acquiring_sequentially = line > last_line_;

    if (is_first_scan_in_slice) {
      biggest_gap_current_ = 0;
    } else if (acquiring_sequentially){
      unsigned int gap = std::abs(static_cast<int>(last_line_) - static_cast<int>(line));
      if (gap != biggest_gap_current_) biggest_gap_current_ = gap;
    } else {
      biggest_gap_current_ = 0;
    }
    last_line_ = line;


    bool is_last_scan_in_slice = m1->isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);

    if (is_last_scan_in_slice && acquiring_sequentially) {
      unsigned int min_ky, max_ky;

      if (biggest_gap_current_ != acceleration_factor_) {
        acceleration_factor_ = biggest_gap_current_;
        weights_invalid_ = true;
      }
 
      if (buffer_counter_.get_region_of_support(min_ky, max_ky) < 0) {
        GDEBUG("Unable to query min_ky, max_ky\n");
        return GADGET_FAIL;
      }
    
      //If there is nothing on the queue, we might as well recalculate
      if (weights_calculator_->msg_queue()->message_count() < 1) {
        //GDEBUG("Queue is empty, invalidating weights\n");
        weights_invalid_ = true;
      } else {
        //GDEBUG("Queue is NOT EMPTY, calculation not triggered\n");
      }

      if (weights_invalid_ && ((max_ky-min_ky) > acceleration_factor_)) {
        std::vector< std::pair<unsigned int, unsigned int> > sampled_region;
        sampled_region.push_back(std::pair<unsigned int, unsigned int>(0, samples-1));
        sampled_region.push_back(std::pair<unsigned int, unsigned int>(min_ky, max_ky));

        std::vector<unsigned int> uncombined_channel_weights;

        //GDEBUG_STREAM("==========================================================================");
        //GDEBUG_STREAM("compute weights on scan : " << m1->scan_counter);
        //GDEBUG("sampled_region[0] = %d,%d\n", sampled_region[0].first, sampled_region[0].second);
        //GDEBUG("sampled_region[1] = %d,%d\n", sampled_region[1].first, sampled_region[1].second);

        if (!weights_calculator_) {
          GDEBUG("Weights calculator not defined\n");
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
}
