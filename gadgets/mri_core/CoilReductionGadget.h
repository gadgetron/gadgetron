#ifndef COILREDUCTIONGADGET_H_
#define COILREDUCTIONGADGET_H_

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{

class EXPORTGADGETSMRICORE CoilReductionGadget :
  public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
    {
    public:
      
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
