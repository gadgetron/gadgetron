/*
 * gpuSenseGadget.h
 *
 *  Created on: Nov 17, 2014
 *      Author: dch
 */

#ifndef GPUSENSEGADGET_H_
#define GPUSENSEGADGET_H_

#include "Gadget.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include "vector_td.h"
#include "GenericReconJob.h"
#include "cuNDArray.h"
namespace Gadgetron {

class gpuSenseGadget: public Gadget2<ISMRMRD::ImageHeader, GenericReconJob>{
public:
  gpuSenseGadget();
  virtual ~gpuSenseGadget();
  virtual int process_config(ACE_Message_Block* mb);
  
protected:
  GADGET_PROPERTY(deviceno,int,"GPU device number", 0);
  GADGET_PROPERTY(setno,int,"Set number to process", 0);
  GADGET_PROPERTY(sliceno,int,"Slice number to process",0);
  GADGET_PROPERTY(oversampling_factor, float, "Oversampling factor for NFFT", 1.5);
  GADGET_PROPERTY(kernel_width, float, "Kernel width for NFFT", 5.5);
  GADGET_PROPERTY(save_individual_frames, bool, "Save individual frames", true);
  GADGET_PROPERTY(output_convergence, bool, "Output convergence information", false);
  GADGET_PROPERTY(rotations_to_discard, int, "Number of rotations to dump", 0);
  GADGET_PROPERTY(output_timing, bool, "Output timing information", false);

  virtual int put_frames_on_que(int frames,int rotations, GenericReconJob* j, cuNDArray<float_complext>* cgresult, int channels = 1);
  int channels_;
  int device_number_;
  int set_number_;
  int slice_number_;
  
  uint64d2 matrix_size_;
  uint64d2 matrix_size_os_;
  uint64d2 matrix_size_seq_;
  double oversampling_factor_;
  double kernel_width_;
  unsigned int rotations_to_discard_;
  
  bool output_convergence_;
  bool output_timing_;
  bool save_individual_frames_;
  
  int frame_counter_;
};

} /* namespace Gadgetron */
#endif /* GPUSENSEGADGET_H_ */
