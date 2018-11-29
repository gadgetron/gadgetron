#ifndef FlowPhaseSubtractionGadget_H
#define FlowPhaseSubtractionGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace Gadgetron{
  
    class EXPORTGADGETSMRICORE FlowPhaseSubtractionGadget :
        public Gadget2< ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {

    public:
        GADGET_DECLARE(FlowPhaseSubtractionGadget);

        FlowPhaseSubtractionGadget();
        virtual ~FlowPhaseSubtractionGadget();

    protected:
        virtual int process_config(ACE_Message_Block* mb);

        virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);

    private:
        unsigned int sets_;

        std::vector<std::vector<float> buffer_;


    };
}

#endif //FlowPhaseSubtractionGadget_H
