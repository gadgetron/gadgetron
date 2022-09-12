#include "CropAndCombineGadget.h"

namespace Gadgetron{
int CropAndCombineGadget::
process( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	 GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{


  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 = 
    new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

  std::vector<size_t> new_dimensions(3);
  new_dimensions[0] = m2->getObjectPtr()->get_size(0)>>1;
  new_dimensions[1] = m2->getObjectPtr()->get_size(1);
  new_dimensions[2] = m2->getObjectPtr()->get_size(2);

  try{m3->getObjectPtr()->create(new_dimensions);}
  catch (std::runtime_error &err){
  	GEXCEPTION(err,"CropAndCombineGadget, failed to allocate new array\n");
    return -1;
  }

  size_t dimx     = m3->getObjectPtr()->get_size(0);
  size_t dimx_old = m2->getObjectPtr()->get_size(0);

  size_t dimy = m3->getObjectPtr()->get_size(1);
  size_t dimz = m3->getObjectPtr()->get_size(2);

  size_t channels = m2->getObjectPtr()->get_size(3);

  std::complex<float>* d1 = m2->getObjectPtr()->get_data_ptr();
  std::complex<float>* d2 = m3->getObjectPtr()->get_data_ptr();

  size_t img_block_old = dimx_old*dimy*dimz;

  for (size_t z = 0; z < dimz; z++) {
    for (size_t y = 0; y < dimy; y++) {
      for (size_t x = 0; x < dimx; x++) {
	float mag = 0;
	float phase = 0;
	size_t offset_1 = z*dimy*dimx_old+y*dimx_old+x+((dimx_old-dimx)>>1);
	size_t offset_2 = z*dimy*dimx+y*dimx+x;
	for (size_t c = 0; c < channels; c++) {
	  float mag_tmp = norm(d1[offset_1 + c*img_block_old]);
	  phase += mag_tmp*arg(d1[offset_1 + c*img_block_old]);
	  mag += mag_tmp;
	}

	d2[offset_2] = std::polar(std::sqrt(mag),phase);
      }
    }
  }

  //Now add the new array to the outgoing message
  m1->cont(m3);
  m2->release();

  //Modify header to match
  m1->getObjectPtr()->matrix_size[0] = m1->getObjectPtr()->matrix_size[0]>>1;
  m1->getObjectPtr()->channels = 1;

  m1->getObjectPtr()->field_of_view[0] = m1->getObjectPtr()->field_of_view[0]/2;

  return this->next()->putq(m1);
}

GADGET_FACTORY_DECLARE(CropAndCombineGadget)
}
