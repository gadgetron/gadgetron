#include "readers/GadgetIsmrmrdReader.h"
#include "ImageWriterGadget.h"

#include <fstream>
namespace Gadgetron{
template<typename T>
int ImageWriterGadget<T> ::
process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1,
	 GadgetContainerMessage< hoNDArray< T > >* m2)
{
    GDEBUG("Writing image\n");

    char filename[1024];
    switch (sizeof(T)) {
     case (8): //Complex float
     	snprintf(filename, 1024, "out_%05d.cplx", (int)this->calls_);
     	break;
     case (4): //Real floats
 		snprintf(filename, 1024, "out_%05d.real", (int)this->calls_);
 		break;
     case (2): //Unsigned short
 		snprintf(filename, 1024, "out_%05d.short", (int)this->calls_);
 		break;
     default:
     	snprintf(filename, 1024, "out_%05d.cplx", (int)this->calls_);
     	break;
     }

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
    outfile.write((char*)m2->getObjectPtr()->get_data_ptr(),sizeof(T)*elements);
    outfile.close();
    delete [] dims;

    this->calls_++;
    return this->next()->putq(m1);
}

GADGET_FACTORY_DECLARE(ImageWriterGadgetUSHORT)
GADGET_FACTORY_DECLARE(ImageWriterGadgetFLOAT)
GADGET_FACTORY_DECLARE(ImageWriterGadgetCPLX)
}
