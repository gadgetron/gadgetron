#include "AccumulatorGadget.h"
#include "GadgetXml.h"

AccumulatorGadget::AccumulatorGadget()
  :buffer_(0)
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

  dimensions_.push_back(GetIntParameterValueFromXML(&doc, "encoding", "readout_length"));
  dimensions_.push_back(GetIntParameterValueFromXML(&doc, "encoding", "matrix_y"));
  dimensions_.push_back(GetIntParameterValueFromXML(&doc, "encoding", "matrix_z"));
  dimensions_.push_back(GetIntParameterValueFromXML(&doc, "encoding", "channels"));
  dimensions_.push_back(GetIntParameterValueFromXML(&doc, "encoding", "slices"));

  if (!(buffer_ = new NDArray< std::complex<float> >())) {
    ACE_DEBUG( (LM_ERROR, ACE_TEXT("Failed to allocate buffer array")) );
    return -1;
  }

  if (!buffer_->create(dimensions_)) {
    ACE_DEBUG( (LM_ERROR, ACE_TEXT("Failed to create buffer array")) );
    return -1;    
  }

  return 0;
}

int AccumulatorGadget::
process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
	GadgetContainerMessage< NDArray< std::complex<float> > >* m2)
{
  ACE_TRACE(( ACE_TEXT("AccumulatorGadget::process") ));

  if (!buffer_) {
    ACE_DEBUG( (LM_ERROR, ACE_TEXT("Accumulator: Buffer array not allocated")) );
    m1->release();
    return -1;    
  }


  std::complex<float>* b = buffer_->get_data_ptr();
  std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();

  
  int samples =  m1->getObjectPtr()->samples;
  int line = m1->getObjectPtr()->idx.line;
  int partition = m1->getObjectPtr()->idx.partition;
  int slice = m1->getObjectPtr()->idx.slice;

  if (samples != dimensions_[0]) {
    ACE_DEBUG( (LM_ERROR, ACE_TEXT("Accumulator: wrong number of samples received")) );
    m1->release();
    return -1;    
  }

  size_t offset= 0;
  //Copy the data for all the channels
  for (int c = 0; c < m1->getObjectPtr()->channels; c++) {
    offset = 
      slice*dimensions_[0]*dimensions_[1]*dimensions_[2]*dimensions_[3] +
      c*dimensions_[0]*dimensions_[1]*dimensions_[2] +
      partition*dimensions_[0]*dimensions_[1] +
      line*dimensions_[0];
    
    memcpy(b+offset,d+c*samples,sizeof(std::complex<float>)*samples);
  }
  
  bool is_last_scan_in_slice =
    (m1->getObjectPtr()->flags & GADGET_FLAG_LAST_ACQ_IN_SLICE);
  
  if (is_last_scan_in_slice) {
    GadgetContainerMessage<GadgetMessageImage>* cm1 = 
      new GadgetContainerMessage<GadgetMessageImage>();
    
    GadgetContainerMessage< NDArray< std::complex<float> > >* cm2 = 
      new GadgetContainerMessage<NDArray< std::complex<float> > >();
    
    cm1->cont(cm2);
    
    std::vector<int> img_dims(4);
    img_dims[0] = dimensions_[0];
    img_dims[1] = dimensions_[1];
    img_dims[2] = dimensions_[2];
    img_dims[3] = dimensions_[3];
    
    if (!cm2->getObjectPtr()->create(img_dims)) {
      ACE_DEBUG( (LM_ERROR, ACE_TEXT("Unable to allocate new image array")) );
      m1->release();
      cm1->release();
      return -1;
    }
    
    size_t data_length = dimensions_[0]*dimensions_[1]*dimensions_[2]*dimensions_[3];
    
    offset = slice*data_length;
    
    memcpy(cm2->getObjectPtr()->get_data_ptr(),b+offset,
	   sizeof(std::complex<float>)*data_length);
    

    cm1->getObjectPtr()->matrix_size[0] = img_dims[0];
    cm1->getObjectPtr()->matrix_size[1] = img_dims[1];
    cm1->getObjectPtr()->matrix_size[2] = img_dims[2];
    cm1->getObjectPtr()->channels       = img_dims[3];
    cm1->getObjectPtr()->data_idx_min       = m1->getObjectPtr()->min_idx;
    cm1->getObjectPtr()->data_idx_max       = m1->getObjectPtr()->max_idx;
    cm1->getObjectPtr()->data_idx_current   = m1->getObjectPtr()->idx;	

    memcpy(cm1->getObjectPtr()->position,m1->getObjectPtr()->position,
	   sizeof(float)*3);

    memcpy(cm1->getObjectPtr()->quarternion,m1->getObjectPtr()->quarternion,
	   sizeof(float)*4);
 
    m1->release();
    
    return this->next()->putq(cm1);
  } 


  m1->release();
  
  return 0;
}

GADGET_FACTORY_DECLARE(AccumulatorGadget)
