#ifndef GPUCGGADGET_H
#define GPUCGGADGET_H

#include "ace/Message_Queue.h"

#include <complex>

#include "Gadgetron.h"
#include "NDArray.h"
#include "GadgetMRIHeaders.h"
#include "GadgetXml.h"

// Cuda includes
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cutil.h>
#include <cufft.h>
#include <cublas.h>
#include <math_constants.h>

// Preprocessing headers
#include "preprocess_sense.hcu"
#include "preprocess_radial.hcu"
#include "preprocess.hcu"
#include "preprocess_private.hcu"

// Other "once to be lib" headers
#include "NSense.hcu"
#include "image_utilities.hcu"
#include "uint_util.hcu"
#include "FFT.hcu"
#include "NFFT.hcu"

class GPUCGGadget : 
public Gadget2<GadgetMessageAcquisition, NDArray< std::complex<float> > >
{
 public:
  GPUCGGadget();
  virtual ~GPUCGGadget();


 protected:
  virtual 
    int process(GadgetContainerMessage<GadgetMessageAcquisition>* m1,
		GadgetContainerMessage< NDArray< std::complex<float> > >* m2);

  virtual int process_config(ACE_Message_Block* mb);

  virtual int set_base_parameters(TiXmlNode* xmlnode);

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
  uint2 matrix_size_;
  uint2 matrix_size_os_;
  unsigned int number_of_iterations_;
  float oversampling_;
  float kernel_width_;
  float kappa_;
  uint2 domain_size_grid_;
  unsigned int domain_size_samples_;
  unsigned int domain_size_coils_;
  uint2 fixed_dims_;
  float gc_factor_;

  bool is_configured_;

  int current_profile_offset_;
  int current_frame_number_;

  int             allocated_samples_;
  float*          data_host_ptr_;
  cuFloatComplex* data_dev_ptr_;
  float2*         trajectory_dev_ptr_;
  float*          dcw_dev_ptr_;

  //cuFloatComplex* coil_images_dev_ptr;
  cuFloatComplex* csm_buffer_dev_ptr_;
  cuFloatComplex* csm_acc_coil_image_os_dev_ptr_; 
  unsigned int kspace_acc_coil_images_os_size_;
  cuFloatComplex* image_dev_ptr_;
  unsigned int csm_buffer_length_;
  
  mr_recon::NFFT_iteration_plan< uint2, float2, 0 > *plan_generic_;

};

#endif //GPUCGGADGET
