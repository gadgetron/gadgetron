#include "GPUCGGoldenRadial.h"
#include "radial_utilities.h"

int GPUCGGoldenRadialGadget::calculate_trajectory()
{
	// Define trajectories
	traj_ = compute_radial_trajectory_golden_ratio_2d<float>
		( samples_per_profile, profiles_per_frame, frames_per_reconstruction, iteration*profiles_per_reconstruction );

  if (!traj_.get()) {
    GADGET_DEBUG1("Failed to allocate trajectory");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

int GPUCGGoldenRadialGadget::calculate_density_compensation()
{

	dcw  = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile, profiles_per_frame, (float)matrix_size_os_.vec[0]/(float)matrix_size_.vec[0], 
      1.0f/((float)samples_per_profile/(float)max(matrix_size_.vec[0],matrix_size_.vec[1])) );
  
  if (!dcw_.get()) {
    GADGET_DEBUG1("Failed to calculate density compensation weights\n");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GPUCGGoldenRadialGadget)
