#ifndef GPUCGGADGET_H
#define GPUCGGADGET_H
#pragma once

#include <complex>

#include "gadgetron_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "vector_td.h"
#include "hoNDArray.h"
#include "cuNDArray.h"
#include "NFFT.h"
#include "cuCG.h"
#include "cgOperatorNonCartesianSense.h"
#include "cgOperatorSenseRHSBuffer.h"
#include "cuCGImageOperator.h"

class EXPORTGADGETSCGSENSE GPUCGGadget : 
public Gadget2<GadgetMessageAcquisition, hoNDArray< std::complex<float> > >
{

public:
  GADGET_DECLARE(GPUCGGadget);

  GPUCGGadget();
  virtual ~GPUCGGadget();

 protected:

  virtual int process_config(ACE_Message_Block* mb);
  virtual int process(GadgetContainerMessage< GadgetMessageAcquisition >* m1,
                      GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);

  virtual int copy_samples_for_profile(float* host_base_ptr,
				       std::complex<float>* data_base_ptr,
				       int profile_no,
				       int channel_no);

  virtual int calculate_trajectory() = 0;
  virtual int calculate_density_compensation() = 0;
  virtual int upload_samples();
  virtual int allocate_csm_buffer();

  ACE_Message_Queue<ACE_MT_SYNCH> buffer_;
  int slice_no_;
  int profiles_per_frame_;
  int shared_profiles_;
  int channels_;
  int samples_per_profile_;
  int device_number_;
  uintd2::Type matrix_size_;
  uintd2::Type matrix_size_os_;
  unsigned int number_of_iterations_;
  float oversampling_;
  float kernel_width_;
  float kappa_;

  bool is_configured_;

  int current_profile_offset_;
  int current_frame_number_;

  int allocated_samples_;

  unsigned int kspace_acc_coil_images_os_size_;
  unsigned int csm_buffer_length_;

  // Thes one is a legacy name from the previous implementation
  // TODO: implement CSM buffer in new design
  cuNDArray<float_complext::Type> csm_acc_coil_image_os_dev_ptr_; 
  
  NFFT_plan< float, 2 > plan_; // used only for csm and regularization estimation, the CG solver holds its own.

  // Define conjugate gradient solver
  cuCG<float, float_complext::Type> cg_;

  // Define non-Cartesian Sense Encofing operator
  boost::shared_ptr< cgOperatorNonCartesianSense<float,2> > E_;

  // Define preconditioner
  boost::shared_ptr< cuCGPrecondWeight<float_complext::Type> > D_;

  // Define regularization image operator
  boost::shared_ptr< cuCGImageOperator<float,float_complext::Type> > R_;

  // Define rhs operator (for regularization)
  boost::shared_ptr< cgOperatorSenseRHSBuffer<float,2> > rhs_buffer_;
    
  // CSM
  boost::shared_ptr< cuNDArray<float_complext::Type> > csm_;

  // Trajectory
  boost::shared_ptr< cuNDArray<floatd2::Type> > traj_;

  // Density compensation weights
  boost::shared_ptr< cuNDArray<float> > dcw_;	

  // Host data array TODO: do we nee dto store this?
  boost::shared_ptr< hoNDArray<float_complext::Type> > host_samples_;
	
  // Device data array
  boost::shared_ptr< cuNDArray<float_complext::Type> > device_sample_s;

};

#endif //GPUCGGADGET
