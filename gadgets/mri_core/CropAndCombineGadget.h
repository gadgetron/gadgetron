#ifndef CROPANDCOMBINEGADGET_H
#define CROPANDCOMBINEGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{

class EXPORTGADGETSMRICORE CropAndCombineGadget :
public Gadget2<ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
{
public:
	CropAndCombineGadget();

protected:

	virtual int process_config( ACE_Message_Block* mb );

	virtual int process( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

	bool no_cropping_;
};


}

#endif //CROPANDCOMBINEGADGET_H
