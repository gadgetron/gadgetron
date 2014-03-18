#ifndef IMAGEFINISHGADGET_H
#define IMAGEFINISHGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "GadgetMRIHeaders.h"
#include "GadgetStreamController.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{

	template <typename T> class EXPORTGADGETSMRICORE ImageFinishGadget : 
	public Gadget2<ISMRMRD::ImageHeader,hoNDArray< T > >
	{
	protected:
		virtual int process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, 
			GadgetContainerMessage< hoNDArray< T > >* m2);
	};

	class EXPORTGADGETSMRICORE ImageFinishGadgetUSHORT :
		public ImageFinishGadget<ACE_UINT16> {};

	class EXPORTGADGETSMRICORE ImageFinishGadgetFLOAT :
		public ImageFinishGadget<float> {};

	class EXPORTGADGETSMRICORE ImageFinishGadgetCPLX :
		public ImageFinishGadget< std::complex<float> > {};
}

#endif //IMAGEFINISHGADGET_H
