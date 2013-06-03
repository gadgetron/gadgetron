#include "gpuGoldenRadialSensePrepGadget.h"
#include "radial_utilities.h"

namespace Gadgetron{

  gpuGoldenRadialSensePrepGadget::gpuGoldenRadialSensePrepGadget() 
    : gpuRadialSensePrepBase() 
  {
    set_parameter(std::string("profiles_per_frame").c_str(), "32");
    this->use_golden_ratio_ = true;
  };

  int gpuGoldenRadialSensePrepGadget::process_config(ACE_Message_Block* mb)
  {
    profiles_per_frame_ = get_int_value(std::string("profiles_per_frame").c_str());
    return gpuRadialSensePrepBase::process_config(mb);
  }

  boost::shared_ptr< cuNDArray<floatd2> > 
  gpuGoldenRadialSensePrepGadget::calculate_trajectory_for_buffer(long profile_offset)
  {
    boost::shared_ptr< cuNDArray<floatd2> > tmp =
      compute_radial_trajectory_golden_ratio_2d<float>( this->samples_per_profile_, this->max_profiles_in_buffer_, 1, profile_offset );

    // The temporary array should be "permuted" to fit the order of the profiles in the buffer:
    // - the first sample in 'tmp' corresponds to the profile at position 'profile_offset' in the buffer 
    // - and 'tmp' then "wraps around" the buffer.

    boost::shared_ptr< cuNDArray<floatd2> > result(new cuNDArray<floatd2>(tmp.get_dimensions()));    
    
    long profile_idx_in_buffer = profile_offset%this->max_profiles_in_buffer_;

    std::vector<unsigned int> head_out_dimensions, tail_out_dimensions;
    head_out_dimensions.push_back(profile_idx_in_buffer*this->samples_per_profile_);
    tail_out_dimensions.push_back((this->max_profiles_in_buffer_-profile_idx_in_buffer)*this->samples_per_profile_);

    cuNDArray<floatd2> head_out( &head_out_dimensions, result.get_data_ptr() );
    cuNDArray<floatd2> tail_out( &tail_out_dimensions, result.get_data_ptr()+profile_idx_in_buffer*this->samples_per_profile_ );

    cuNDArray<floatd2> head_in( &tail_out_dimensions, tmp.get_data_ptr() );
    cuNDArray<floatd2> tail_in( &head_out_dimensions, tmp.get_data_ptr()+(this->max_profiles_in_buffer_-profile_idx_in_buffer)*this->samples_per_profile_ );

    head_out = tail_in;
    tail_out = head_in;

    return result;
  }

  boost::shared_ptr< hoNDArray<floatd2> > 
  gpuGoldenRadialSensePrepGadget::calculate_trajectory_for_reconstruction(long profile_offset)
  {
    return compute_radial_trajectory_fixed_angle_2d<float>
      ( this->samples_per_profile_, this->profiles_per_frame_, 
	(this->use_multiframe_grouping_) ? this->acceleration_factor_ : 1, profile_offset )->to_host();
  }

  boost::shared_ptr< cuNDArray<float> >   
  gpuGoldenRadialSensePrepGadget::calculate_density_compensation_for_buffer()
  {
    // The weights do not change even if the trajectory rotates
    return compute_radial_dcw_fixed_angle_2d<float>
      ( this->samples_per_profile_, this->max_profiles_in_buffer_, this->oversampling_factor_, 
	1.0f/(float(this->samples_per_profile_)/float(this->image_dimensions_[0])) );
  }

  boost::shared_ptr< hoNDArray<float> >   
  gpuGoldenRadialSensePrepGadget::calculate_density_compensation_for_reconstruction()
  {
    // The weights do not change even if the trajectory rotates
    return compute_radial_dcw_fixed_angle_2d<float>
      ( this->samples_per_profile_, this->profiles_per_frame_, this->oversampling_factor_, 
	1.0f/(float(this->samples_per_profile_)/float(this->image_dimensions_[0])) )->to_host();
  }
  
  GADGET_FACTORY_DECLARE(gpuGoldenRadialSensePrepGadget)
}
