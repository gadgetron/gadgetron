#include "AccumulatorGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

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

/**
 *   Expects ISMRMRD XML configuration
 *
 */
int AccumulatorGadget::process_config(ACE_Message_Block* mb)
{

	boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

	ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
	if (e_seq.size() != 1) {
		GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
		GADGET_DEBUG1("This simple AccumulatorGadget only supports one encoding space\n");
		return GADGET_FAIL;
	}

	ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
	ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
	ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

	GADGET_DEBUG2("Matrix size: %d, %d, %d\n", e_space.matrixSize().x(), e_space.matrixSize().y(), e_space.matrixSize().z());
	dimensions_.push_back(e_space.matrixSize().x());
	dimensions_.push_back(e_space.matrixSize().y());
	dimensions_.push_back(e_space.matrixSize().z());

	slices_ = e_limits.slice().present() ? e_limits.slice().get().maximum()+1 : 1;

  return GADGET_OK;
}

int AccumulatorGadget::
process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
	GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  if (!buffer_) {
	  dimensions_.push_back(m1->getObjectPtr()->active_channels);
	  dimensions_.push_back(slices_);

	  if (!(buffer_ = new hoNDArray< std::complex<float> >())) {
		  GADGET_DEBUG1("Failed create buffer\n");
		  return GADGET_FAIL;
	  }

	  if (!buffer_->create(&dimensions_)) {
		  GADGET_DEBUG1("Failed allocate buffer array\n");
		  return GADGET_FAIL;
	  }

	  image_series_ = this->get_int_value("image_series");

  }


  std::complex<float>* b =
		  buffer_->get_data_ptr();

  std::complex<float>* d =
		  m2->getObjectPtr()->get_data_ptr();

  int samples =  m1->getObjectPtr()->number_of_samples;
  int line = m1->getObjectPtr()->idx.kspace_encode_step_1;
  int partition = m1->getObjectPtr()->idx.kspace_encode_step_2;
  int slice = m1->getObjectPtr()->idx.slice;

  if (samples > static_cast<int>(dimensions_[0])) {
	  GADGET_DEBUG1("Wrong number of samples received\n");
	  return GADGET_FAIL;
  }

  size_t offset= 0;
  //Copy the data for all the channels
  for (int c = 0; c < m1->getObjectPtr()->active_channels; c++) {
    offset = 
      slice*dimensions_[0]*dimensions_[1]*dimensions_[2]*dimensions_[3] +
      c*dimensions_[0]*dimensions_[1]*dimensions_[2] +
      partition*dimensions_[0]*dimensions_[1] +
      line*dimensions_[0] + (dimensions_[0]>>1)-m1->getObjectPtr()->center_sample;
    
    memcpy(b+offset,
    	d+c*samples,
    	sizeof(std::complex<float>)*samples);
  }
  
  bool is_last_scan_in_slice = ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags);
  
  if (is_last_scan_in_slice) {
    GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = 
      new GadgetContainerMessage<ISMRMRD::ImageHeader>();
    
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
    cm1->getObjectPtr()->slice   = m1->getObjectPtr()->idx.slice;

    memcpy(cm1->getObjectPtr()->position,
    		m1->getObjectPtr()->position,
	   sizeof(float)*3);

    memcpy(cm1->getObjectPtr()->quaternion,
    		m1->getObjectPtr()->quaternion,
	   sizeof(float)*4);
 
    memcpy(cm1->getObjectPtr()->patient_table_position,
    		m1->getObjectPtr()->patient_table_position, sizeof(float)*3);

    cm1->getObjectPtr()->image_data_type = ISMRMRD::DATA_COMPLEX_FLOAT;
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
