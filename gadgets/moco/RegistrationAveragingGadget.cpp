#include "RegistrationAveragingGadget.h"
#include "Gadgetron.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "GadgetronTimer.h"
#include <numeric>
#include "Spline.h"


#ifdef USE_OMP
#include <omp.h>
#endif 

namespace Gadgetron{

  template<class ARRAY_TYPE, unsigned int D> 
  RegistrationAveragingGadget<ARRAY_TYPE,D>::RegistrationAveragingGadget() 
  {
    set_parameter(std::string("alpha").c_str(), "0.05");
    set_parameter(std::string("beta").c_str(), "1.0");
  }

  template<class ARRAY_TYPE, unsigned int D> 
  RegistrationAveragingGadget<ARRAY_TYPE,D>::~RegistrationAveragingGadget() {}
  
  template<class ARRAY_TYPE, unsigned int D> int 
  RegistrationAveragingGadget<ARRAY_TYPE,D>::process_config( ACE_Message_Block *mb )
  {
    alpha_ = ARRAY_TYPE::element_type(get_double_value("alpha"));
    beta_  = ARRAY_TYPE::element_type(get_double_value("beta"));
    return GADGET_OK;
  }

  template<class ARRAY_TYPE, unsigned int D> int 
  RegistrationAveragingGadget<ARRAY_TYPE,D>::process( GadgetContainerMessage<ISMRMRD::ImageHeader> *m1,
						      GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2 )
  {
    
    //
    // Put the incoming images on the appropriate queue 
    // - based on the phase index
    //

    // At first pass allocate the buffer array
    //

    if( phase_images_.get() == 0x0 ){
      unsigned int 
    }

    
    unsigned int phase 
    
    if (buffer_.enqueue_tail(m3) < 0) {
      GADGET_DEBUG1("Failed to add image to buffer\n");
      m3->release();
      return GADGET_FAIL;
    }

    time_stamps_.push_back(m1->getObjectPtr()->physiology_time_stamp[phys_time_index_]);

    if (this->next()->putq(m1) < 0) {
      GADGET_DEBUG1("Unable to put data on next Gadgets Q\n");
      return GADGET_FAIL;
    }

    return GADGET_OK;
  }

  int RegistrationAveragingGadget::close( unsigned long flags ) 
  {
    GADGET_DEBUG1("RegistrationAveragingGadget::close...\n");

    int ret = Gadget::close(flags);
    return ret;
  }

  GADGET_FACTORY_DECLARE(RegistrationAveragingGadget)
}
