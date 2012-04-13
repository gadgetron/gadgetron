#include "GPUCGGoldenRadial.h"
#include "radial_utilities.h"
#include "Gadgetron.h"

boost::shared_ptr< cuNDArray<floatd2> >
GPUCGGoldenRadialGadget::calculate_trajectory()
{
	// Define trajectories
	boost::shared_ptr< cuNDArray<floatd2> > traj = compute_radial_trajectory_golden_ratio_2d<float>
		( samples_per_profile_, profiles_per_frame_, 1, current_profile_offset_ );

  if (!traj.get()) {
    GADGET_DEBUG1("Failed to allocate trajectory");
    return boost::shared_ptr< cuNDArray<floatd2> >();
  }

  return traj;
}

boost::shared_ptr< cuNDArray<float> >
GPUCGGoldenRadialGadget::calculate_density_compensation()
{

	boost::shared_ptr< cuNDArray<float> > dcw = compute_radial_dcw_golden_ratio_2d
    ( samples_per_profile_, profiles_per_frame_, (float)matrix_size_os_.vec[0]/(float)matrix_size_.vec[0], 
      1.0f/((float)samples_per_profile_/(float)std::max(matrix_size_.vec[0],matrix_size_.vec[1])) );
  
  if (!dcw.get()) {
    GADGET_DEBUG1("Failed to calculate density compensation weights\n");
    return boost::shared_ptr< cuNDArray<float> >();
  }

  return dcw;
}

GADGET_FACTORY_DECLARE(GPUCGGoldenRadialGadget)
