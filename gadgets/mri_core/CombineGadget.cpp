#include "CombineGadget.h"

namespace Gadgetron{

  CombineGadget::CombineGadget() {}
  CombineGadget::~CombineGadget() {}

int CombineGadget::
process( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	 GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  // Get the dimensions
  size_t nx = m2->getObjectPtr()->get_size(0);
  size_t ny = m2->getObjectPtr()->get_size(1);
  size_t nz = m2->getObjectPtr()->get_size(2);
  size_t nc = m2->getObjectPtr()->get_size(3);

  // Create a new message with an hoNDArray for the combined image
  GadgetContainerMessage< hoNDArray<std::complex<float> > >* m3 = 
    new GadgetContainerMessage< hoNDArray<std::complex<float> > >();

  std::vector<size_t> dimensions(3);
  dimensions[0] = nx;
  dimensions[1] = ny; 
  dimensions[2] = nz;

  try{m3->getObjectPtr()->create(dimensions);}
  catch (std::runtime_error &err){
  	GEXCEPTION(err,"CombineGadget, failed to allocate new array\n");
    return -1;
  }

  std::complex<float>* d1 = m2->getObjectPtr()->get_data_ptr();
  std::complex<float>* d2 = m3->getObjectPtr()->get_data_ptr();

  size_t img_block = nx*ny*nz;

  for (size_t z = 0; z < nz; z++) {
    for (size_t y = 0; y < ny; y++) {
      for (size_t x = 0; x < nx; x++) {
	float mag = 0;
	float phase = 0;
	size_t offset = z*ny*nx+y*nx+x;
	for (size_t c = 0; c < nc; c++) {
	  float mag_tmp = norm(d1[offset + c*img_block]);
	  phase += mag_tmp*arg(d1[offset + c*img_block]);
	  mag += mag_tmp;
	}
	d2[offset] = std::polar(std::sqrt(mag),phase/mag);
      }
    }
  }

  // Modify header to match the size and change the type to real
  m1->getObjectPtr()->channels = 1;

  // Now add the new array to the outgoing message
  m1->cont(m3);

  // Release the old data
  m2->release();

  return this->next()->putq(m1);
}

GADGET_FACTORY_DECLARE(CombineGadget)
}
