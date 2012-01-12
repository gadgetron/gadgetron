#ifndef SPIRALGADGET_H
#define SPIRALGADGET_H

#include <complex>
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "spiral_gadgets_export.h"
#include "vector_td.h"
#include "NFFT.h"

#include <boost/shared_ptr.hpp>


class EXPORTGADGETSSPIRAL SpiralGadget : 
public Gadget2< GadgetMessageAcquisition, hoNDArray< std::complex<float> > >
{
  
 public:
  GADGET_DECLARE(SpiralGadget);

  SpiralGadget();
  virtual ~SpiralGadget();

 protected:
  virtual int process_config(ACE_Message_Block* mb);
  virtual int process(GadgetContainerMessage< GadgetMessageAcquisition >* m1,
		      GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);

 private:
  int samples_to_skip_start_;
  int samples_to_skip_end_;
  int samples_per_interleave_;
  boost::shared_ptr< hoNDArray<floatd2::Type> > host_traj_;
  boost::shared_ptr< hoNDArray<float> > host_weights_;
  cuNDArray<float> gpu_weights_;

  hoNDArray<float_complext::Type>* host_data_buffer_;
  std::vector<unsigned int> image_dimensions_;
  NFFT_plan<float, 2> plan_;

};

#endif //SPIRALGADGET_H
