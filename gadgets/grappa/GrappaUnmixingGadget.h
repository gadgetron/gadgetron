#ifndef GRAPPAUNMIXINGGADGET_H_
#define GRAPPAUNMIXINGGADGET_H_

#include "gadgetron_grappa_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd.h"
#include "GrappaWeights.h"

#include <complex>

namespace Gadgetron{

struct EXPORTGADGETSGRAPPA GrappaUnmixingJob
{
	boost::shared_ptr< GrappaWeights<float> > weights_;
};

class EXPORTGADGETSGRAPPA GrappaUnmixingGadget: public Gadget3<GrappaUnmixingJob, ISMRMRD::ImageHeader, hoNDArray<std::complex<float> > > {
public:
	GADGET_DECLARE(GrappaUnmixingGadget);

	GrappaUnmixingGadget();
	virtual ~GrappaUnmixingGadget();
protected:
	virtual int process(GadgetContainerMessage<GrappaUnmixingJob>* m1,
			GadgetContainerMessage<ISMRMRD::ImageHeader>* m2, GadgetContainerMessage<hoNDArray<std::complex<float> > >* m3);

};
}
#endif /* GRAPPAUNMIXINGGADGET_H_ */
