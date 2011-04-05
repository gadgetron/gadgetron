#ifndef IMAGEWRITER_H
#define IMAGEWRITER_H

#include <fstream>

#include "GadgetSocketReceiver.h"

class ImageWriter : public GadgetImageMessageReader
{

 public:
  ImageWriter() 
    : number_of_calls_(0)
    {}


  virtual ACE_Message_Block* read(ACE_SOCK_Stream* socket) {
    ACE_Message_Block* mb = GadgetImageMessageReader::read(socket);

    if (!mb) {
      GADGET_DEBUG1("Read failed in parent\n");
      return 0;
    }
    
    GadgetContainerMessage<GadgetMessageImage> * img_head_mb = 
      dynamic_cast<GadgetContainerMessage<GadgetMessageImage> *>(mb);

    if (!img_head_mb) {
      GADGET_DEBUG1("Failed in dynamic cast\n");
      mb->release();
      return 0;
    }


    GadgetContainerMessage<hoNDArray< std::complex<float> > > * img_data_mb = 
      dynamic_cast<GadgetContainerMessage<hoNDArray< std::complex<float> > > *>(img_head_mb->cont());

    if (!img_data_mb) {
      GADGET_DEBUG1("Failed in dynamic cast\n");
      mb->release();
      return 0;
    }
    
    if (this->process_image(img_head_mb->getObjectPtr(), img_data_mb->getObjectPtr()) < 0) {
      GADGET_DEBUG1("Failed to process image\n");
      mb->release();
      return 0;
    }

    return mb;
  }
  
  virtual int process_image(GadgetMessageImage* img_head, 
			    hoNDArray< std::complex<float> >* data)
  {
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("Image Writer writing image\n")) );

    char filename[1024];
    sprintf(filename, "out_%05d.cplx", (int)number_of_calls_);

    std::ofstream outfile;    
    outfile.open (filename, std::ios::out|std::ios::binary);

    int ndim = 3;
    int dims[3];
    size_t elements = 1;
    dims[0] = img_head->matrix_size[0]; elements*=dims[0];
    dims[1] = img_head->matrix_size[1]; elements*=dims[1];
    dims[2] = img_head->matrix_size[2]; elements*=dims[2];

    outfile.write((char*)&ndim,sizeof(int));
    outfile.write((char*)dims,sizeof(int)*3);
    outfile.write((char*)data->get_data_ptr(),sizeof(float)*elements*2);
    outfile.close();

    number_of_calls_++;

    return 0;
  }

 protected:
  size_t number_of_calls_;

};

#endif //IMAGE_WRITER


