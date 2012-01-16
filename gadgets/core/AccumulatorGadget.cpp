#include "AccumulatorGadget.h"
#include "GadgetXml.h"

AccumulatorGadget::AccumulatorGadget()
  :buffer_(0)
  , image_counter_(0)
  , image_series_(0)
{

}
 
AccumulatorGadget::~AccumulatorGadget()
{
  if (buffer_) delete buffer_;
}

int AccumulatorGadget::process_config(ACE_Message_Block* mb)
{

  TiXmlDocument doc;
  doc.Parse(mb->rd_ptr());

  GadgetXMLNode n =
		  GadgetXMLNode(&doc).get<GadgetXMLNode>(std::string("gadgetron"))[0];

  std::vector<long> dims =
		  n.get<long>(std::string("encoding.kspace.matrix_size.value"));

  if (dims.size() < 3) {
	dims.push_back(0);
  }

  if (dims[2] == 0) {
    dims[2] = 1;
  }

  dimensions_.push_back(
		  n.get<long>(std::string("encoding.kspace.readout_length.value"))[0]);
  dimensions_.push_back(dims[1]);
  dimensions_.push_back(dims[2]);
  dimensions_.push_back(
		  n.get<long>(std::string("encoding.channels.value"))[0]);
  dimensions_.push_back(
		  n.get<long>(std::string("encoding.slices.value"))[0]);


  if (!(buffer_ = new hoNDArray< std::complex<float> >())) {
	  GADGET_DEBUG1("Failed create buffer\n");
	  return GADGET_FAIL;
  }

  if (!buffer_->create(&dimensions_)) {
	  GADGET_DEBUG1("Failed allocate buffer array\n");
	  return GADGET_FAIL;
  }

  image_series_ = this->get_int_value("image_series");

  return GADGET_OK;
}

int AccumulatorGadget::
process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  if (!buffer_) {
	  GADGET_DEBUG1("Buffer not allocated\n");
	  return GADGET_FAIL;
  }


  std::complex<float>* b =
		  buffer_->get_data_ptr();

  std::complex<float>* d =
		  m2->getObjectPtr()->get_data_ptr();

  int samples =  m1->getObjectPtr()->samples;
  int line = m1->getObjectPtr()->idx.line;
  int partition = m1->getObjectPtr()->idx.partition;
  int slice = m1->getObjectPtr()->idx.slice;

  if (samples != static_cast<int>(dimensions_[0])) {
	  GADGET_DEBUG1("Wrong number of samples received\n");
	  return GADGET_FAIL;
  }

  size_t offset= 0;
  //Copy the data for all the channels
  for (int c = 0; c < m1->getObjectPtr()->channels; c++) {
    offset = 
      slice*dimensions_[0]*dimensions_[1]*dimensions_[2]*dimensions_[3] +
      c*dimensions_[0]*dimensions_[1]*dimensions_[2] +
      partition*dimensions_[0]*dimensions_[1] +
      line*dimensions_[0];
    
    memcpy(b+offset,
    	d+c*samples,
    	sizeof(std::complex<float>)*samples);
  }
  
  bool is_last_scan_in_slice =
    (m1->getObjectPtr()->flags & GADGET_FLAG_LAST_ACQ_IN_SLICE);
  
  if (is_last_scan_in_slice) {
    GadgetContainerMessage<GadgetMessageImage>* cm1 = 
      new GadgetContainerMessage<GadgetMessageImage>();
    
    cm1->getObjectPtr()->flags = 0;

    GadgetContainerMessage< hoNDArray< std::complex<float> > >* cm2 = 
      new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
    
    cm1->cont(cm2);
    
    std::vector<unsigned int> img_dims(4);
    img_dims[0] = dimensions_[0];
    img_dims[1] = dimensions_[1];
    img_dims[2] = dimensions_[2];
    img_dims[3] = dimensions_[3];
    
    if (!cm2->getObjectPtr()->create(&img_dims)) {
      GADGET_DEBUG1("Unable to allocate new image array\n");
      cm1->release();
      return -1;
    }
    
    size_t data_length = dimensions_[0]*dimensions_[1]*
    		dimensions_[2]*dimensions_[3];
    
    offset = slice*data_length;
    
    memcpy(cm2->getObjectPtr()->get_data_ptr(),b+offset,
	   sizeof(std::complex<float>)*data_length);
    
    cm1->getObjectPtr()->matrix_size[0]     = img_dims[0];
    cm1->getObjectPtr()->matrix_size[1]     = img_dims[1];
    cm1->getObjectPtr()->matrix_size[2]     = img_dims[2];
    cm1->getObjectPtr()->channels           = img_dims[3];
    cm1->getObjectPtr()->data_idx_min       = m1->getObjectPtr()->min_idx;
    cm1->getObjectPtr()->data_idx_max       = m1->getObjectPtr()->max_idx;
    cm1->getObjectPtr()->data_idx_current   = m1->getObjectPtr()->idx;	

    memcpy(cm1->getObjectPtr()->position,
    		m1->getObjectPtr()->position,
	   sizeof(float)*3);

    memcpy(cm1->getObjectPtr()->quarternion,
    		m1->getObjectPtr()->quarternion,
	   sizeof(float)*4);
 
    cm1->getObjectPtr()->table_position =
    		m1->getObjectPtr()->table_position;

    cm1->getObjectPtr()->image_format = GADGET_IMAGE_COMPLEX_FLOAT;
    cm1->getObjectPtr()->image_index = ++image_counter_;
    cm1->getObjectPtr()->image_series_index = image_series_;

    if (this->next()->putq(cm1) < 0) {
    	return GADGET_FAIL;
    }
  } 

  m1->release();
  return GADGET_OK;
}

GADGET_FACTORY_DECLARE(AccumulatorGadget)
