#include "GPUCGGoldenRadial.h"
#include "radial_utilities.h"
#include "Gadgetron.h"

int GPUCGGoldenRadialGadget::calculate_trajectory()
{
	// Define trajectories
	traj_ = compute_radial_trajectory_golden_ratio_2d<float>
		( samples_per_profile_, profiles_per_frame_, 1, current_profile_offset_ );

  if (!traj_.get()) {
    GADGET_DEBUG1("Failed to allocate trajectory");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

int GPUCGGoldenRadialGadget::calculate_density_compensation()
{

	dcw_ = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile_, profiles_per_frame_, (float)matrix_size_os_.vec[0]/(float)matrix_size_.vec[0], 
      1.0f/((float)samples_per_profile_/(float)std::max(matrix_size_.vec[0],matrix_size_.vec[1])) );
  
  if (!dcw_.get()) {
    GADGET_DEBUG1("Failed to calculate density compensation weights\n");
    return GADGET_FAIL;
  }

  return GADGET_OK;
}

GADGET_FACTORY_DECLARE(GPUCGGoldenRadialGadget)
