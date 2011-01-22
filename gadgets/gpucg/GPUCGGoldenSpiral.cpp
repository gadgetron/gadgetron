#include "GPUCGGoldenSpiral.h"
#include "ConfigParser.h"
#include <iostream>
#include <fstream>

void calc_traj(double* xgrad, double* ygrad, 
	       int ngrad, int Nints, double Tgsamp, 
	       double krmax, double** x_trajectory, double** y_trajectory,
	       double** weights);


void calc_vds(double slewmax,double gradmax,double Tgsample,
	      double Tdsample,int Ninterleaves, double* fov, 
	      int numfov,double krmax, int ngmax, 
	      double** xgrad,double** ygrad,int* numgrad);

#ifndef PI
   #define PI 3.14159265
#endif

GPUCGGoldenSpiralGadget::GPUCGGoldenSpiralGadget()
  : Interleaves_(0)
  , ADCsPerInterleave_(0)
  , SamplesPerADC_(0)
  , SamplesToSkipStart_(0)
  , SamplesToSkipEnd_(0)
  , SamplingTime_ns_(0)
  , Reordering_(0)
  , MaxGradient_Gcm_(0.0f)
  , MaxSlewRate_Gcms_(0.0f)
  , krmax_cm_(0.0f)
  , FOVCoeff_1_(0.0f)
  , host_allocated_traj_samples_(0)
  , host_trajectory_ptr_(0)
  , host_density_weight_ptr_(0)
  , allocated_dev_traj_samples_(0)
{
 
  profiles_per_frame_ = 16;
  shared_profiles_ = 0;
}


GPUCGGoldenSpiralGadget::~GPUCGGoldenSpiralGadget()
{
  if (host_trajectory_ptr_) delete [] host_trajectory_ptr_;
  if (host_density_weight_ptr_) delete [] host_density_weight_ptr_;
}

int GPUCGGoldenSpiralGadget::set_base_parameters(ConfigParser* cp)
{
  
  profiles_per_frame_ = Interleaves_;

  channels_ = cp->getIntVal("encoding","channels");

  if (matrix_size_.x == 0 && matrix_size_.y == 0) {
    matrix_size_ = make_uint2(cp->getIntVal("encoding","matrix_x"), 
			      cp->getIntVal("encoding","matrix_y"));
  }

  return GADGET_OK;
}

int GPUCGGoldenSpiralGadget::process_config(ACE_Message_Block* mb)
{
  ConfigParser cp;
  cp.parse(mb->rd_ptr());


  if (!is_configured_) {
    if (!cp.findSection(std::string("spiral"))) {
      GADGET_DEBUG1("Unable to locate spiral section of configuration.\n");
      return GADGET_FAIL;
    }
    
    Interleaves_        = cp.getIntVal(std::string("spiral"), 
				       std::string("Interleaves"));
    
    ADCsPerInterleave_  = cp.getIntVal(std::string("spiral"),
				       std::string("ADCsPerInterleave"));
    
    SamplesPerADC_      = cp.getIntVal(std::string("spiral"),
				       std::string("SamplesPerADC"));
    
    SamplesToSkipStart_ =  cp.getIntVal(std::string("spiral"),
					std::string("SamplesToSkipStart"));
    
    SamplesToSkipEnd_   =  cp.getIntVal(std::string("spiral"),
					std::string("SamplesToSkipEnd"));
    
    SamplingTime_ns_    =  cp.getIntVal(std::string("spiral"),
					std::string("SamplingTime_ns"));
    
    Reordering_         =  cp.getIntVal(std::string("spiral"),
					std::string("Reordering"));
    
    MaxGradient_Gcm_    =  cp.getFloatVal(std::string("spiral"),
					  std::string("MaxGradient_Gcm"));
    
    MaxSlewRate_Gcms_   =  cp.getFloatVal(std::string("spiral"),
					  std::string("MaxSlewRate_Gcms"));
    
    
    krmax_cm_           =  cp.getFloatVal(std::string("spiral"),
					  std::string("krmax_cm"));
    
    
    FOVCoeff_1_         =  cp.getFloatVal(std::string("spiral"),
					  std::string("FOVCoeff_1"));


    //Calculate trajectory
    {

      double  Tdsamp = SamplingTime_ns_/(1.0e9); /*	Data Sample period (s) */
      //double  Tgsamp = 1e-5;                         /*	Data Sample period (s) */

      double  *xgrad; 	/* 	X-component of gradient.	*/
      double  *ygrad;     /*	Y-component of gradient.	*/
      double  *x_trajectory;
      double  *y_trajectory;
      double  *weighting;
      int     ngrad;	

      //Calc gradient waveform
      calc_vds(MaxSlewRate_Gcms_,MaxGradient_Gcm_,
	       Tdsamp,Tdsamp,
	       Interleaves_,&FOVCoeff_1_,1,
	       krmax_cm_,1e5,
	       &xgrad,&ygrad,&ngrad);
      

      GADGET_DEBUG2("Number of gradient samples: %d\n", ngrad);
      if ((ADCsPerInterleave_*SamplesPerADC_-ngrad) != SamplesToSkipEnd_) {

	ngrad = (ADCsPerInterleave_*SamplesPerADC_-SamplesToSkipEnd_) < ngrad ? 
	  ADCsPerInterleave_*SamplesPerADC_-SamplesToSkipEnd_ : ngrad;

	SamplesToSkipEnd_ = ADCsPerInterleave_*SamplesPerADC_-ngrad;
      }
      GADGET_DEBUG2("Number of gradient samples: %d\n", ngrad);

      /* Calcualte the trajectory and weights*/
      calc_traj(xgrad, ygrad, ngrad, 
		Interleaves_, Tdsamp, krmax_cm_, &x_trajectory, &y_trajectory, &weighting);

      samples_per_profile_ = ngrad;
      host_trajectory_ptr_ = new float2[ngrad*Interleaves_];
      host_density_weight_ptr_ = new float[ngrad*Interleaves_];
      host_allocated_traj_samples_ = ngrad*Interleaves_;

      if (!host_trajectory_ptr_ || ! host_density_weight_ptr_) {
	GADGET_DEBUG1("Unable to allocate memory for trajectory and weights\n");
	return GADGET_FAIL;
      }

      float2 maxk = make_float2(0.0f, 0.0f);
      float2 mink = make_float2(0.0f, 0.0f);

      for (int i = 0; i < (ngrad*Interleaves_); i++) {
	host_trajectory_ptr_[i].x = -x_trajectory[i]/2;
	host_trajectory_ptr_[i].y = -y_trajectory[i]/2;
	host_density_weight_ptr_[i] = weighting[i];
	if (host_trajectory_ptr_[i].x > maxk.x) 
	  maxk.x = host_trajectory_ptr_[i].x;
	if (host_trajectory_ptr_[i].y > maxk.y) 
	  maxk.y = host_trajectory_ptr_[i].y;
	if (host_trajectory_ptr_[i].x < mink.x) 
	  mink.x = host_trajectory_ptr_[i].x;
	if (host_trajectory_ptr_[i].y < mink.y) 
	  mink.y = host_trajectory_ptr_[i].y;
      }

      GADGET_DEBUG2("Trajectory calculated: maxk.x = %f, maxk.y = %f, mink.x = %f, mink.y = %f\n",  maxk.x, maxk.y, mink.x, mink.y);

      delete [] x_trajectory;
      delete [] y_trajectory;
      delete [] weighting;
      delete [] xgrad;
      delete [] ygrad;
    }
    
    profiles_per_frame_ = Interleaves_;
    shared_profiles_ = 0;
  }

  return GPUCGGadget::process_config(mb);
}

int  GPUCGGoldenSpiralGadget::copy_samples_for_profile(float* host_base_ptr,
						       std::complex<float>* data_base_ptr,
						       int profile_no,
						       int channel_no)
{
  
  memcpy(host_base_ptr + 
	 (channel_no*allocated_samples_ + profile_no*samples_per_profile_) * 2,
	 data_base_ptr + channel_no*SamplesPerADC_, 
	 sizeof(float)*samples_per_profile_*2);
  
  return GADGET_OK;
}

int GPUCGGoldenSpiralGadget::calculate_trajectory()
{
   cudaError_t err;
  if (allocated_dev_traj_samples_ != Interleaves_*samples_per_profile_) {
     if (trajectory_dev_ptr_) {
       cudaFree(trajectory_dev_ptr_);
       trajectory_dev_ptr_ = 0;
     }

     if (dcw_dev_ptr_) {
       cudaFree(dcw_dev_ptr_);
       dcw_dev_ptr_ = 0;
     }

     cudaMalloc( (void**) &trajectory_dev_ptr_, 
		 samples_per_profile_*Interleaves_*sizeof(float2) );

     err = cudaGetLastError();
     if( err != cudaSuccess ){
       GADGET_DEBUG2("Failed to allocate memory for trajectory: %s\n",
		     cudaGetErrorString(err));
      return GADGET_FAIL;
    }


     cudaMalloc( (void**) &dcw_dev_ptr_, 
		 samples_per_profile_*Interleaves_*sizeof(float) );
     
     err = cudaGetLastError();
     if( err != cudaSuccess ){
       GADGET_DEBUG2("Failed to allocate memory for density weights: %s\n",
		     cudaGetErrorString(err));
      return GADGET_FAIL;
    }

     allocated_dev_traj_samples_ = samples_per_profile_*Interleaves_;

  }

  if (Reordering_ != 16) { //If this is not golden angle    
    //This one is easy, we just copy the stuff
    cudaMemcpy(trajectory_dev_ptr_,
	       host_trajectory_ptr_,
	       host_allocated_traj_samples_*sizeof(float2),
	       cudaMemcpyHostToDevice);

  } else {
    //This must be a golden angle....we need to deal with that

    float2* tmp_traj = new float2[samples_per_profile_*Interleaves_];
    if (!tmp_traj) {
      GADGET_DEBUG1("Failed to allocate temporary host memory for trajectory\n");
      return GADGET_FAIL;
    }


    for (int i = 0; i < Interleaves_; i++) {
      float rotation = (current_profile_offset_ + i) * 2* (PI / ((sqrt(5.0) + 1)/2));

      float cos_rot = cos(rotation);
      float sin_rot = sin(rotation);
      for (int s = 0; s < samples_per_profile_; s++) {
	tmp_traj[i*samples_per_profile_ + s].x = 
	  host_trajectory_ptr_[s].x * cos_rot +
	  host_trajectory_ptr_[s].y * sin_rot;

	tmp_traj[i*samples_per_profile_ + s].y = 
	  -host_trajectory_ptr_[s].x * sin_rot +
	  host_trajectory_ptr_[s].y * cos_rot;
      }
    }

    cudaMemcpy(trajectory_dev_ptr_,
	       tmp_traj,
	       host_allocated_traj_samples_*sizeof(float2),
	       cudaMemcpyHostToDevice);


    delete [] tmp_traj;
  }

  err = cudaGetLastError();
  if( err != cudaSuccess ){
    GADGET_DEBUG2("Failed to upload trajectory: %s\n",
		  cudaGetErrorString(err));
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

int GPUCGGoldenSpiralGadget::calculate_density_compensation()
{
  if (!dcw_dev_ptr_) {
    GADGET_DEBUG1("Memory has not been allocated for density compensation\n");
    return GADGET_FAIL;
  }

  cudaMemcpy(dcw_dev_ptr_,
	     host_density_weight_ptr_,
	     host_allocated_traj_samples_*sizeof(float),
	     cudaMemcpyHostToDevice);
  
  cudaError_t err = cudaGetLastError();
  if( err != cudaSuccess ){
    GADGET_DEBUG2("Failed to upload density compensation weights: %s\n",
		  cudaGetErrorString(err));
    return GADGET_FAIL;
  }


  return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GPUCGGoldenSpiralGadget)
