#include <ace/OS_NS_stdlib.h>

#include "Gadgetron.h"
#include "GrappaGadget.h"
#include "GadgetXml.h"
#include "FFT.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

GrappaGadget::GrappaGadget()
{

}
 
GrappaGadget::~GrappaGadget()
{
  for (unsigned int i = 0; i < buffers_.size(); i++) {
    if (buffers_[i]) delete buffers_[i];
    buffers_[i] = 0;

    if (weights_[i]) delete weights_[i];
    weights_[i] = 0;
  }

}

int GrappaGadget::process_config(ACE_Message_Block* mb)
{
  TiXmlDocument doc;
  doc.Parse(mb->rd_ptr());

  GadgetXMLNode n = GadgetXMLNode(&doc).get<GadgetXMLNode>(std::string("gadgetron"))[0];

  std::vector<long> dims = n.get<long>(std::string("encoding.kspace.matrix_size.value"));

  if (dims.size() < 3) {
    GADGET_DEBUG2("Matrix dimensions have the wrong length: %d\n", dims.size());
    return GADGET_FAIL;
  }

  if (dims[2] == 0) {
    dims[2] = 1;
  }

  dimensions_.push_back(n.get<long>(std::string("encoding.kspace.readout_length.value"))[0]);
  //  dimensions_.push_back(dims[0]);
  dimensions_.push_back(dims[1]);
  dimensions_.push_back(dims[2]);
  dimensions_.push_back(n.get<long>(std::string("encoding.channels.value"))[0]);
  dimensions_.push_back(n.get<long>(std::string("encoding.slices.value"))[0]);

  image_dimensions_.push_back(dimensions_[0] / 2); //TODO: fix this in general
  image_dimensions_.push_back(dimensions_[1]);
  image_dimensions_.push_back(dimensions_[2]);
  image_dimensions_.push_back(dimensions_[3]);


  weights_ = std::vector< GrappaWeights<float>* >(dimensions_[4],0);
  buffers_ = std::vector<GrappaCalibrationBuffer* >(dimensions_[4],0);
  time_stamps_ = std::vector<ACE_UINT32>(dimensions_[4],0);

 
 //Let's figure out if we have channels that are supposed to be uncombined
  boost::shared_ptr<std::string> uncomb_str = this->get_string_value("uncombined_channels");
  std::vector<std::string> uncomb;
  boost::split(uncomb, *uncomb_str, boost::is_any_of(","));
  for (unsigned int i = 0; i < uncomb.size(); i++) {
    std::string ch = boost::algorithm::trim_copy(uncomb[i]);
    if (ch.size() > 0) {
      unsigned int channel_id = static_cast<unsigned int>(ACE_OS::atoi(ch.c_str()));
      weights_calculator_.add_uncombined_channel(channel_id);
    }
  }

  for (unsigned int i = 0; i < buffers_.size(); i++) {
    weights_[i] = new GrappaWeights<float>();

    //Let's set some default GRAPPA weights, so that we have something to work with the first couple of frames.
    std::vector<unsigned int> wdims = image_dimensions_;
    if (weights_calculator_.get_number_of_uncombined_channels()) {
      wdims.push_back(weights_calculator_.get_number_of_uncombined_channels()+1);
    }

    hoNDArray< std::complex<float> > tmp_w;
    if (!tmp_w.create(&wdims)) {
      GADGET_DEBUG1("Unable to create temporary array with dimensions\n");
      return GADGET_FAIL;
    }
    tmp_w.clear(std::complex<float>(1.0,0));
    weights_[i]->update(&tmp_w);
    

    buffers_[i] = new GrappaCalibrationBuffer(image_dimensions_,
					      weights_[i],
					      &weights_calculator_);
  }

 
  if (weights_calculator_.open() < 0) {
    GADGET_DEBUG1("Failed to open GrappaWeightsCalculator\n");
    return GADGET_FAIL;
  }

  image_data_ = std::vector< GadgetContainerMessage< hoNDArray< std::complex<float> > >* >(dimensions_[4],0);
  for (unsigned int i = 0; i < image_data_.size(); i++) {
    if (create_image_buffer(i) != GADGET_OK) {
      GADGET_DEBUG1("Unable to create image buffers");
      return GADGET_FAIL;
    }
  }

  return GADGET_OK;
}

int GrappaGadget::
process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  GadgetMessageAcquisition* acq_head = m1->getObjectPtr();
  
  unsigned int samples =  acq_head->samples;
  unsigned int line = acq_head->idx.line;
  unsigned int partition = acq_head->idx.partition;
  unsigned int slice = acq_head->idx.slice;

  if (samples != image_dimensions_[0]) {
    GADGET_DEBUG1("GrappaGadget: wrong number of samples received\n");
    return GADGET_FAIL;    
  }

  if (slice >= image_data_.size()) {
    GADGET_DEBUG1("Invalid slice number received\n");
    return GADGET_FAIL;
  }

  if (!image_data_[0]) {
    if (create_image_buffer(slice) != GADGET_OK) {
      GADGET_DEBUG1("Failed to allocate new slice buffer\n");
      return GADGET_FAIL;
    }
  }

  std::complex<float>* b = image_data_[slice]->getObjectPtr()->get_data_ptr();
  std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();

  size_t offset= 0;
  //Copy the data for all the channels
  for (int c = 0; c < m1->getObjectPtr()->channels; c++) {
    offset = 
      c*image_dimensions_[0]*image_dimensions_[1]*image_dimensions_[2] +
      partition*image_dimensions_[0]*image_dimensions_[1] +
      line*image_dimensions_[0];
    
    memcpy(b+offset,d+c*samples,sizeof(std::complex<float>)*samples);
  }


  bool is_last_scan_in_slice =
    (m1->getObjectPtr()->flags & GADGET_FLAG_LAST_ACQ_IN_SLICE);
  
  bool is_first_scan_in_slice =
    (m1->getObjectPtr()->flags & GADGET_FLAG_FIRST_ACQ_IN_SLICE);

  if (is_first_scan_in_slice) {
    time_stamps_[slice] = m1->getObjectPtr()->time_stamp;
  }

  if (is_last_scan_in_slice) {

    GadgetContainerMessage<GadgetMessageImage>* cm1 = 
      new GadgetContainerMessage<GadgetMessageImage>();
    
    GadgetContainerMessage< hoNDArray<std::complex<float> > >* cm2 = 
      new GadgetContainerMessage< hoNDArray<std::complex<float> > >();
    
    
    std::vector<unsigned int> combined_dims(3,0);
    combined_dims[0] = image_dimensions_[0];
    combined_dims[1] = image_dimensions_[1];
    combined_dims[2] = image_dimensions_[2];

    if (weights_calculator_.get_number_of_uncombined_channels()) {
      combined_dims.push_back(weights_calculator_.get_number_of_uncombined_channels()+1); 
    }
 
    if (!cm2->getObjectPtr()->create(&combined_dims)) {
      GADGET_DEBUG1("Unable to create combined image array\n");
      return GADGET_FAIL;
    }

    cm1->cont(cm2);
    

    cm1->getObjectPtr()->matrix_size[0] = image_dimensions_[0];
    cm1->getObjectPtr()->matrix_size[1] = image_dimensions_[1];
    cm1->getObjectPtr()->matrix_size[2] = image_dimensions_[2];
    cm1->getObjectPtr()->channels       = 1+weights_calculator_.get_number_of_uncombined_channels();
    cm1->getObjectPtr()->data_idx_min       = m1->getObjectPtr()->min_idx;
    cm1->getObjectPtr()->data_idx_max       = m1->getObjectPtr()->max_idx;
    cm1->getObjectPtr()->data_idx_current   = m1->getObjectPtr()->idx;	
    cm1->getObjectPtr()->time_stamp         = time_stamps_[slice];

    memcpy(cm1->getObjectPtr()->position,m1->getObjectPtr()->position,
	   sizeof(float)*3);

    memcpy(cm1->getObjectPtr()->quarternion,m1->getObjectPtr()->quarternion,
	   sizeof(float)*4);


    FFT<float>::instance()->ifft(image_data_[slice]->getObjectPtr(),0);
    FFT<float>::instance()->ifft(image_data_[slice]->getObjectPtr(),1);
    FFT<float>::instance()->ifft(image_data_[slice]->getObjectPtr(),2);

    //apply weights
    float scale_factor = dimensions_[0]*dimensions_[1];

    int appl_result = weights_[slice]->apply(image_data_[slice]->getObjectPtr(), cm2->getObjectPtr(), scale_factor);
    if (appl_result < 0) {
      GADGET_DEBUG2("Failed to apply GRAPPA weights: error code %d\n", appl_result);
      return GADGET_FAIL;
    }

    if (this->next()->putq(cm1) < 0) {
      GADGET_DEBUG1("Failed to pass image on to next Gadget in chain\n");
      return GADGET_FAIL;
    }

    image_data_[slice]->getObjectPtr()->clear(std::complex<float>(0.0f,0.0f));
  }

  if (buffers_[slice]->add_data(m1->getObjectPtr(),m2->getObjectPtr()) < 0) {
    GADGET_DEBUG1("Failed to add incoming data to grappa calibration buffer\n");
    return GADGET_FAIL;
  }

  m1->release();
  return GADGET_OK;
}


int GrappaGadget::create_image_buffer(unsigned int slice)
{
  if (slice >= image_data_.size()) {
    return GADGET_FAIL;
  }

  if (image_data_[slice] != 0) {
    image_data_[slice]->release();
    image_data_[slice] = 0;
  }

  image_data_[slice] = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
  if (!image_data_[slice]->getObjectPtr()->create(&image_dimensions_)) {
    GADGET_DEBUG1("Unable to create image buffers");
    return GADGET_FAIL;
  } 

  std::fill(image_data_[slice]->getObjectPtr()->get_data_ptr(),
	    image_data_[slice]->getObjectPtr()->get_data_ptr()+image_data_[slice]->getObjectPtr()->get_number_of_elements(),
	    std::complex<float>(0.0f,0.0f));

  return GADGET_OK;

}

GADGET_FACTORY_DECLARE(GrappaGadget)
