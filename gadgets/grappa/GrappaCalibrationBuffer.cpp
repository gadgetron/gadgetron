#include "GrappaCalibrationBuffer.h"

#include "Gadgetron.h"

GrappaCalibrationBuffer::GrappaCalibrationBuffer(std::vector<unsigned int> dimensions,
						 GrappaWeights<float>* w,
						 GrappaWeightsCalculator<float>* weights_calculator)
  : weights_(w)
  , weights_calculator_(weights_calculator)
  , buffer_counter_(dimensions[1])
{
  if (!buffer_.create(dimensions)) {
    GADGET_DEBUG1("Unable to allocate memory for GRAPPA buffer");
  }
  dimensions_ = dimensions;
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

  if (buffer_counter_.update_line(line,m1->position,m1->quarternion) < 0) {
    GADGET_DEBUG1("Unable to update buffer counter\n");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}
