#include "GPUCGGoldenRadial.h"


int GPUCGGoldenRadialGadget::calculate_trajectory()
{
  if (trajectory_dev_ptr_) {
    cudaFree(trajectory_dev_ptr_);
    trajectory_dev_ptr_ = 0;
  }
  
  if (!trajectory_dev_ptr_) {
    trajectory_dev_ptr_ = compute_trajectory_radial_2d<1>(matrix_size_.x, 
							  matrix_size_os_.x, 
							  samples_per_profile_, 
							  profiles_per_frame_, 
							  current_profile_offset_, 
							  1,
							  gc_factor_);
  }

  if (!trajectory_dev_ptr_) {
    GADGET_DEBUG1("Failed to allocate trajectory");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

int GPUCGGoldenRadialGadget::calculate_density_compensation()
{

  //TODO: add check to see if we really need to recalculate this.
  if (dcw_dev_ptr_) {
    cudaFree(dcw_dev_ptr_);
    dcw_dev_ptr_ = 0;
  }

  if (!dcw_dev_ptr_) {
    dcw_dev_ptr_ = compute_dcw_radial_2d<1>( matrix_size_.x, 
					     matrix_size_os_.x, 
					     samples_per_profile_, 
					     profiles_per_frame_, 
					     current_profile_offset_, 
					     1,
					     gc_factor_);
  }

  if (!dcw_dev_ptr_) {
    GADGET_DEBUG1("Failed to calculate density compensation weights\n");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GPUCGGoldenRadialGadget)
