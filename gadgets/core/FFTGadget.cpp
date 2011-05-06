#include "FFTGadget.h"
#include "FFT.h"

int FFTGadget::process( GadgetContainerMessage< GadgetMessageImage>* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{
  FFT<float>::instance()->ifft(m2->getObjectPtr(),0);
  FFT<float>::instance()->ifft(m2->getObjectPtr(),1);
  FFT<float>::instance()->ifft(m2->getObjectPtr(),2);

  std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();
  for (unsigned int i = 0; 
       i <  m2->getObjectPtr()->get_number_of_elements(); 
       i++) 
    {
      d[i] *= m2->getObjectPtr()->get_number_of_elements();
    } 

  return this->next()->putq(m1);
}

GADGET_FACTORY_DECLARE(FFTGadget)
