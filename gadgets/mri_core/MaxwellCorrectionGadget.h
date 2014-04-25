#ifndef MaxwellCorrectionGadget_H
#define MaxwellCorrectionGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{  

    class EXPORTGADGETSMRICORE MaxwellCorrectionGadget :
        public Gadget2< ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {

    public:
        GADGET_DECLARE(MaxwellCorrectionGadget);
        MaxwellCorrectionGadget();
        virtual ~MaxwellCorrectionGadget();


    protected:
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);
	
    private:
	std::vector<double> maxwell_coefficients_;
	bool maxwell_coefficients_present_;
    };
}

#endif //MaxwellCorrectionGadget_H
