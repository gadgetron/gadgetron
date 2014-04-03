#include "GadgetIsmrmrdReadWrite.h"
#include "CombineGadget.h"

namespace Gadgetron{

  CombineGadget::CombineGadget() {}
  CombineGadget::~CombineGadget() {}

int CombineGadget::
process( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	 GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  // Get the dimensions
  int nx = m2->getObjectPtr()->get_size(0);
  int ny = m2->getObjectPtr()->get_size(1);
  int nz = m2->getObjectPtr()->get_size(2);
  int nc = m2->getObjectPtr()->get_size(3);

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

  // Warp the input and output in armadillo column vectors
  // Reshape the input into a 2D array (coil is the outermost dimension)
  arma::cx_fmat adata_in  = arma::cx_fmat( m2->getObjectPtr()->get_data_ptr(), nx*ny*nz, nc, false, true );
  arma::cx_fmat adata_out = arma::cx_fmat( m3->getObjectPtr()->get_data_ptr(), nx*ny*nz, 1,  false, true );

  // Initialize the output to zero
  adata_out.zeros();

  // Square root of the sum of squares over the coil dimension
  // divide by the number of channels to keep noise variance scaled properly
  // TODO: check the scale
  adata_out.set_real((1.0/((float)nc)) * arma::sqrt(arma::sum(arma::square(arma::abs(adata_in)),1)));

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
