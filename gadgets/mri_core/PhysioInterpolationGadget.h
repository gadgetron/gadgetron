#ifndef PhysioInterpolationGadget_H
#define PhysioInterpolationGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{  
    class EXPORTGADGETSMRICORE  PhysioInterpolationGadget :
        public Gadget2< ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
    {

    public:
        GADGET_DECLARE(PhysioInterpolationGadget);

        PhysioInterpolationGadget();
        virtual ~PhysioInterpolationGadget();

    protected:
        virtual int process_config(ACE_Message_Block* mb);

        virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1,
            GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);
	
	virtual int close(unsigned long flags); //All the work is done here in this Gadget

	unsigned short phys_time_index_;
	unsigned short phases_to_reconstruct_;
	unsigned short mode_; //0=seperate series for each complete RR,
	                      //1=First complete RR interval only
	

    private:
	ACE_Message_Queue<ACE_MT_SYNCH> buffer_;
	std::vector<float> time_stamps_;
    };
}

#endif //PhysioInterpolationGadget_H
