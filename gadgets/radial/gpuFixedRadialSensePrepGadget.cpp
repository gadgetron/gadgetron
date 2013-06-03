#include "gpuFixedRadialSensePrepGadget.h"
#include "radial_utilities.h"

#include <math.h>
#define _USE_MATH_DEFINES

namespace Gadgetron{

  boost::shared_ptr< cuNDArray<floatd2> > 
  gpuFixedRadialSensePrepGadget::calculate_trajectory_for_buffer(long profile_offset)
  {    
    return compute_radial_trajectory_fixed_angle_2d<float>
      ( this->samples_per_profile_, this->max_profiles_in_buffer_, 1, profile_offset );
  }

  boost::shared_ptr< hoNDArray<floatd2> > 
  gpuFixedRadialSensePrepGadget::calculate_trajectory_for_reconstruction(long profile_offset)
  {
    float angular_offset;

    if(this->use_multiframe_grouping_) 
      angular_offset = 0.0f;
    else{
      long local_frame = (profile_offset/this->profiles_per_frame_)%this->aceleration_factor_;
      float angular_offset = M_PI/float(this->profiles_per_frame_)*float(local_frame)/float(this->acceleration_factor_);
    }

    return compute_radial_trajectory_fixed_angle_2d<float>
      ( this->samples_per_profile_, this->profiles_per_frame_, 
	this->frames_per_reconstruction_, angular_offset )->to_host();
  }

  boost::shared_ptr< cuNDArray<float> >   
  gpuFixedRadialSensePrepGadget::calculate_density_compensation_for_buffer()
  {
    // The weights do not change even if the trajectory rotates
    return compute_radial_dcw_fixed_angle_2d<float>
      ( this->samples_per_profile_, this->max_profiles_in_buffer_, this->oversampling_factor_, 
	1.0f/(float(this->samples_per_profile_)/float(this->image_dimensions_[0])) );
  }

  boost::shared_ptr< hoNDArray<float> >   
  gpuFixedRadialSensePrepGadget::calculate_density_compensation_for_reconstruction()
  {
    // The weights do not change even if the trajectory rotates
    return compute_radial_dcw_fixed_angle_2d<float>
      ( this->samples_per_profile_, this->profiles_per_frame_, this->oversampling_factor_, 
	1.0f/(float(this->samples_per_profile_)/float(this->image_dimensions_[0])) )->to_host();
  }
  
  GADGET_FACTORY_DECLARE(gpuFixedRadialSensePrepGadget)
}
