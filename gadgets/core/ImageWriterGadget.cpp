#include "ImageWriterGadget.h"

#include <fstream>

int ImageWriterGadget ::
process( GadgetContainerMessage< GadgetMessageImage>* m1,
	 GadgetContainerMessage< NDArray< std::complex<float> > >* m2)
{
    ACE_DEBUG( (LM_DEBUG, ACE_TEXT("ImageWriterGadget writing image\n")) );

    char filename[1024];
    sprintf(filename, "out_%05d.cplx", (int)calls_);

    std::ofstream outfile;    
    outfile.open (filename, std::ios::out|std::ios::binary);

    int ndim = m2->getObjectPtr()->get_number_of_dimensions();
    int* dims = new int[ndim];
    size_t elements = 1;
    for (int d = 0; d < ndim; d++) {
      dims[d] = m2->getObjectPtr()->get_size(d);
      elements *= dims[d];
    }
    outfile.write((char*)&ndim,sizeof(int));
    outfile.write((char*)dims,sizeof(int)*ndim);
    outfile.write((char*)m2->getObjectPtr()->get_data_ptr(),sizeof(float)*elements*2);
    outfile.close();
    delete [] dims;

    calls_++;
    return this->next()->putq(m1);
}
