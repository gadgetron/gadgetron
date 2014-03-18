#ifndef FlowPhaseSubtractionGadget_H
#define FlowPhaseSubtractionGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace Gadgetron{
  
    class EXPORTGADGETSMRICORE FlowPhaseSubtractionGadget :
        public Gadget2< ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {

    public:

        FlowPhaseSubtractionGadget();
        virtual ~FlowPhaseSubtractionGadget();

    protected:
        virtual int process_config(ACE_Message_Block* mb);

        virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);

    private:
        unsigned int sets_;
	boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > buffer_;
    };
}

#endif //FlowPhaseSubtractionGadget_H
