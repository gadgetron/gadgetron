#pragma once

#include "gadgetroncore_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd.h"

#include <complex>
namespace Gadgetron{
class EXPORTGADGETSCORE NoiseAdjustGadget : 
public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
{
 public:
  GADGET_DECLARE(NoiseAdjustGadget);
  
  NoiseAdjustGadget();

 protected:
  bool noise_decorrelation_calculated_;
  hoNDArray< std::complex<double> > noise_covariance_matrix_;
  unsigned long int number_of_noise_samples_;
  float noise_dwell_time_us_;
  float acquisition_dwell_time_us_;
  float noise_bw_scale_factor_;
  float receiver_noise_bandwidth_;
  bool is_configured_;

  virtual int process_config(ACE_Message_Block* mb);
  virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
		      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

  
};

}
