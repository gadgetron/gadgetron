#ifndef RegistrationAveragingGadget_H
#define RegistrationAveragingGadget_H

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "complext.h"
#include "gadgetron_moco_export.h"

#include <ismrmrd.h>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace Gadgetron{  

  /**
     This is an abstract gadget class and consequently should not be included in any xml configuration file.
     "Instatiate" instead the cpuRegistrationAveragingGadget or gpuRegistrationAveragingGadget.
  */
  template<class ARRAY_TYPE, unsigned int D> class EXPORTGADGETS_MOCO RegistrationAveragingGadget :
    public Gadget2<ISMRMRD::ImageHeader, ARRAY_TYPE>
  {
    
  public:
    GADGET_DECLARE(RegistrationAveragingGadget);

    RegistrationAveragingGadget();
    virtual ~RegistrationAveragingGadget();

  protected:
    virtual int process_config(ACE_Message_Block *mb);

    virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader > *m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2);
	
    virtual int close(unsigned long flags); // All the work is done here in this gadget

    virtual int setup_solver() = 0;

  private:
    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> phase_images_;
    boost::shared_ptr< opticalFlowSolver<ARRAY_TYPE,D> > of_solver_;
    ARRAY_TYPE::element_type alpha_;
    ARRAY_TYPE::element_type beta_;
  };
}

#endif //RegistrationAveragingGadget_H
