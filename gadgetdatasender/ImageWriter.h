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

  virtual int process_image(GadgetMessageImage* img_head, 
			    NDArray< std::complex<float> >* data)
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


