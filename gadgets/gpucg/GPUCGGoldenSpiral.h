#ifndef GPUCGGOLDENSPIRALGADGET_H
#define GPUCGGOLDENSPIRALGADGET_H

#include "GPUCGGadget.h"

class GPUCGGoldenSpiralGadget : public GPUCGGadget
{

 public:
  GPUCGGoldenSpiralGadget(bool pass_on_data = false, int slice = 0);
  ~GPUCGGoldenSpiralGadget();

 protected:
  virtual int set_base_parameters(ConfigParser* cp);
  virtual int process_config(ACE_Message_Block* mb);

  virtual int copy_samples_for_profile(float* host_base_ptr,
				       std::complex<float>* data_base_ptr,
				       int profile_no,
				       int channel_no);

  virtual int calculate_trajectory();
  virtual int calculate_density_compensation();


  int Interleaves_;
  int ADCsPerInterleave_;
  int SamplesPerADC_;
  int SamplesToSkipStart_;
  int SamplesToSkipEnd_; 
  int SamplingTime_ns_;
  int Reordering_;
  double MaxGradient_Gcm_;
  double MaxSlewRate_Gcms_;
  double krmax_cm_;
  double FOVCoeff_1_;

  int     host_allocated_traj_samples_;
  float2* host_trajectory_ptr_;
  float*  host_density_weight_ptr_;
  int allocated_dev_traj_samples_;

};

#endif //GPUCGGOLDENSPIRALGADGET_H
