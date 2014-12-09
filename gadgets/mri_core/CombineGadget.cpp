#include "GadgetIsmrmrdReadWrite.h"
#include "CombineGadget.h"
#include "mri_core.h"
#include "hoNDArray_math.h"

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

  try{m3->getObjectPtr()->create(&dimensions);}
  catch (std::runtime_error &err){
  	GADGET_DEBUG_EXCEPTION(err,"CombineGadget, failed to allocate new array\n");
    return -1;
  }

  // Create a temporary real array to store the magnitude (ssos)
  hoNDArray< float > temp(nx,ny,nz);

  // sqrt of sum of squares
  coilmap_norm(*m2->getObjectPtr(), temp);

  // assign to the real part of the complex result
  real_to_complex(temp, *m3->getObjectPtr());
  
  // Modify header to match the size
  m1->getObjectPtr()->channels = 1;

  // Now add the new array to the outgoing message
  m1->cont(m3);

  // Release the old data
  m2->release();

  return this->next()->putq(m1);
}

GADGET_FACTORY_DECLARE(CombineGadget)
}
