#include "AccumulatorGadget.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{
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
  ISMRMRD::IsmrmrdHeader h;
  ISMRMRD::deserialize(mb->rd_ptr(),h);

  if (h.encoding.size() != 1) {
    GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
    GDEBUG("This simple AccumulatorGadget only supports one encoding space\n");
    return GADGET_FAIL;
  }


  ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
  ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
  ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;
  
  GDEBUG("Matrix size: %d, %d, %d\n", r_space.matrixSize.x, e_space.matrixSize.y, e_space.matrixSize.z);
  dimensions_.push_back(r_space.matrixSize.x);
  dimensions_.push_back(e_space.matrixSize.y);
  dimensions_.push_back(e_space.matrixSize.z);
  
  field_of_view_.push_back(r_space.fieldOfView_mm.x);
  field_of_view_.push_back(e_space.fieldOfView_mm.y);
  field_of_view_.push_back(e_space.fieldOfView_mm.z);
  GDEBUG("FOV: %f, %f, %f\n", r_space.fieldOfView_mm.x, e_space.fieldOfView_mm.y, e_space.fieldOfView_mm.z);
  
  slices_ = e_limits.slice? e_limits.slice->maximum+1 : 1;

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
		  GDEBUG("Failed create buffer\n");
		  return GADGET_FAIL;
	  }

	  try {buffer_->create(dimensions_);}
	  catch (std::runtime_error &err){
		  GEXCEPTION(err,"Failed allocate buffer array\n");
		  return GADGET_FAIL;
	  }

	  image_series_ = image_series.value();

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
	  GDEBUG("Wrong number of samples received\n");
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
  
  bool is_last_scan_in_slice = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
  
  if (is_last_scan_in_slice) {
    GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 = 
      new GadgetContainerMessage<ISMRMRD::ImageHeader>();

    // On some platforms, it is necessary to initialize the image header
    memset(cm1->getObjectPtr(),0,sizeof(ISMRMRD::ImageHeader));
    
    cm1->getObjectPtr()->clearAllFlags();

    GadgetContainerMessage< hoNDArray< std::complex<float> > >* cm2 = 
      new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
    
    cm1->cont(cm2);
    
    std::vector<size_t> img_dims(4);
    img_dims[0] = dimensions_[0];
    img_dims[1] = dimensions_[1];
    img_dims[2] = dimensions_[2];
    img_dims[3] = dimensions_[3];
    
    try{cm2->getObjectPtr()->create(img_dims);}
    catch (std::runtime_error &err){
      GEXCEPTION(err,"Unable to allocate new image array\n");
      cm1->release();
      return -1;
    }
    
    size_t data_length = dimensions_[0]*dimensions_[1]*
    		dimensions_[2]*dimensions_[3];
    
    offset = slice*data_length;
    
    memcpy(cm2->getObjectPtr()->get_data_ptr(),b+offset,
	   sizeof(std::complex<float>)*data_length);
    
    cm1->getObjectPtr()->matrix_size[0]     = (uint16_t)img_dims[0];
    cm1->getObjectPtr()->matrix_size[1]     = (uint16_t)img_dims[1];
    cm1->getObjectPtr()->matrix_size[2]     = (uint16_t)img_dims[2];
    cm1->getObjectPtr()->field_of_view[0]   = field_of_view_[0];
    cm1->getObjectPtr()->field_of_view[1]   = field_of_view_[1];
    cm1->getObjectPtr()->field_of_view[2]   = field_of_view_[2];
    cm1->getObjectPtr()->channels           = (uint16_t)img_dims[3];
    cm1->getObjectPtr()->slice   = m1->getObjectPtr()->idx.slice;

    memcpy(cm1->getObjectPtr()->position,
    		m1->getObjectPtr()->position,
	   sizeof(float)*3);

    memcpy(cm1->getObjectPtr()->read_dir,
                m1->getObjectPtr()->read_dir,
           sizeof(float)*3);

    memcpy(cm1->getObjectPtr()->phase_dir,
                m1->getObjectPtr()->phase_dir,
           sizeof(float)*3);

    memcpy(cm1->getObjectPtr()->slice_dir,
                m1->getObjectPtr()->slice_dir,
           sizeof(float)*3);

    memcpy(cm1->getObjectPtr()->patient_table_position,
    		m1->getObjectPtr()->patient_table_position, sizeof(float)*3);

    cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;
    cm1->getObjectPtr()->image_index = (uint16_t)(++image_counter_);
    cm1->getObjectPtr()->image_series_index = (uint16_t)image_series_;

    if (this->next()->putq(cm1) < 0) {
    	return GADGET_FAIL;
    }
  } 

  m1->release();
  return GADGET_OK;
}

GADGET_FACTORY_DECLARE(AccumulatorGadget)
}
