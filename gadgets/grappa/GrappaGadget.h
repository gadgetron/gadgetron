#ifndef GRAPPAGADGET_H
#define GRAPPAGADGET_H

#include <complex>

#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "GrappaCalibrationBuffer.h"

struct GrappaBufferInfo
{
  float           position[3];
  float           quarternion[4];   
  unsigned int    acceleration_factor;
};

class GrappaGadget : 
public Gadget2< GadgetMessageAcquisition, hoNDArray< std::complex<float> > >
{
  
 public:
  GADGET_DECLARE(AccumulatorGadget);

  GrappaGadget();
  ~GrappaGadget();

 protected:
  virtual int process_config(ACE_Message_Block* mb);
  virtual int process(GadgetContainerMessage< GadgetMessageAcquisition >* m1,
		      GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);

  virtual int create_image_buffer(unsigned int slice);

  std::vector< GrappaCalibrationBuffer* > buffers_;
  std::vector<unsigned int> dimensions_;
  std::vector<unsigned int> image_dimensions_;
  std::vector< GadgetContainerMessage<  hoNDArray< std::complex<float> > >* > image_data_;
  std::vector< GrappaWeights<float>* > weights_;
  GrappaWeightsCalculator<float> weights_calculator_;
};

#endif //GRAPPAGADGET_H
