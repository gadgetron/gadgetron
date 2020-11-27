#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <ctime>
#include "GadgetMRIHeaders.h"
#include "ismrmrd/meta.h"

namespace Gadgetron
{
    class EXPORTGADGETSMRICORE NoiseSummaryGadget : public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
        GADGET_DECLARE(NoiseSummaryGadget);
        
        using BaseClass = Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >;
        
        NoiseSummaryGadget();
        virtual ~NoiseSummaryGadget();

        virtual int close(unsigned long flags);

    protected:
	GADGET_PROPERTY(noise_file, std::string, "Name of noise file", "");

        std::string noise_dependency_folder_;
        bool processed_in_close_;

        virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
                            GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
        {
            m1->release();
            return GADGET_OK;
        }
    };
}
