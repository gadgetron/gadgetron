#ifndef COILREDUCTIONGADGET_H_
#define COILREDUCTIONGADGET_H_

#include "gadgetron_core_export.h"
#include "Gadget.h"
#include "ismrmrd.h"
#include "hoNDArray.h"

#include <complex>

namespace Gadgetron{

class EXPORTGADGETSCORE CoilReductionGadget :
public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
{
public:
	GADGET_DECLARE(CoilReductionGadget);

	CoilReductionGadget();
	virtual ~CoilReductionGadget();

	virtual int process_config(ACE_Message_Block* mb);
	virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);

protected:
	std::vector<unsigned short> coil_mask_;
	unsigned int coils_in_;
	unsigned int coils_out_;

};
}
#endif /* COILREDUCTIONGADGET_H_ */
