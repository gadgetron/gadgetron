#ifndef GRAPPAGADGET_H
#define GRAPPAGADGET_H

#include <complex>

#include "gadgetrongrappa_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "GrappaCalibrationBuffer.h"
#include "ismrmrd.h"
namespace Gadgetron{
struct EXPORTGADGETSGRAPPA GrappaBufferInfo
{
  float           position[3];
  float           quaternion[4];
  unsigned int    acceleration_factor;
};

class EXPORTGADGETSGRAPPA GrappaGadget : 
public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
{
  
 public:
  GADGET_DECLARE(GrappaGadget);

  GrappaGadget();
  virtual ~GrappaGadget();

 protected:
  virtual int process_config(ACE_Message_Block* mb);
  virtual int process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1,
		  GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2 );

  virtual int create_image_buffer(unsigned int slice);

  //We have to overwrite close in this gadget to make sure we wait for the weights calculator.
  virtual int close(unsigned long flags);

  virtual int initial_setup();

  bool first_call_;
 private:
  std::vector< GrappaCalibrationBuffer* > buffers_;
  std::vector<unsigned int> dimensions_;
  std::vector<unsigned int> image_dimensions_;
  std::vector< GadgetContainerMessage<  hoNDArray< std::complex<float> > >* > image_data_;
  std::vector< boost::shared_ptr<GrappaWeights<float> > > weights_;
  GrappaWeightsCalculator<float> weights_calculator_;
  std::vector<ACE_UINT32> time_stamps_;
  int image_counter_;
  int image_series_;
  int target_coils_;
  float phase_encoding_resolution_;
  unsigned int line_offset_;
};
}
#endif //GRAPPAGADGET_H
