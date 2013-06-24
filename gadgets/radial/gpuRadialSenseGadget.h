#pragma once

#include "gadgetron_radial_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "vector_td.h"
#include "cuNFFT.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgSolver.h"
#include "cuCgPreconditioner.h"
#include "cuSenseBufferCg.h"
#include "ismrmrd.h"

#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

/*
  ------------------------------------------
  Trajectory modes for radial reconstruction
  ------------------------------------------
  
  Mode 0 and Mode 1 are variants of 'fixed' radial trajectories with interframe rotation.
  Mode 2 denotes a radial trajectory with an angular profile spacing based on the golden ratio (~111,25 degrees).
  
  Let 
  'i' denote the number of profiles per (undersampled) frame
  'j' denote the number of frames per trajectory rotation (to obtain a fully sampled acquisition)
  'h' denote a variable of type ISMRMRD::AcquisitionHeader

  It is possible to explicitly set 'i' and 'j' in the Gadgetron configuration file.
  For some modes this is (partly) required, 
  for others they will be automatically determined from the incoming profile headers.
  
  Mode 0:
  -------
  For each rotation cycle profiles are numbered using the scheme

    0+0*j,0+1*j,0+2*j,...,0+(i-1)*j, (1st frame)
    1+0*j,1+1*j,1+2*j,...,1+(i-1)*j, (2nd frame)
    2+0*j,2+1*j,2+2*j,...,2+(i-1)*j, (3rd frame)
    ...,
    (j-1)+0*j,(j-1)+1*j,(j-1)+2*j,...,(j-1)+(i-1)*j

  as given in h.idx.kspace_encode_step_1.
  Both 'i' and 'j' are automatically derived and thus need not be explicitly specified in a configuration file.
  For mode 0 both 'i' and 'j' can be changed dynamically as desired e.g. for real-time imaging.

  Mode 1:
  -------
  Profiles are numbered 0,1,2,...,i-1, 0,1,2,...,i-1, ... as given in h.idx.kspace_encode_step_1.
  'j' is estimated as 'matrix_size'/'i' and should be explicitly set in the configuration file if this is not the case, e.g.:
  <property><name>frames_per_rotation</name><value>8</value></property>
      

  Mode 2:
  -------
  Profiles are numbered 
  0,1,2,...,i-1, 0,1,2,...,i-1, 0,1,2,...,i-1, ...
  or
  0,1,2,...,i-1, i,i+1,i+2,...,2*i-1, 2*i,2*i+1,2*i+2,3*i-1, ...
  as given in h.idx.kspace_encode_step_1.
  'i' should be explicitly specified in the Gadgetron configuration file, e.g.:
  <property><name>profiles_per_frame</name><value>32</value></property>
  If not it defaults to i=32.
  'j' is explicitly set to '1' even if specified in the configuration file.
*/

namespace Gadgetron{

  class EXPORTGADGETS_RADIAL gpuRadialSenseGadget :
    public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
  {

  public:
    GADGET_DECLARE(gpuRadialSenseGadget);

    gpuRadialSenseGadget();
    virtual ~gpuRadialSenseGadget();

  protected:
    
    virtual int process_config(ACE_Message_Block* mb);

    virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2);

  private:

    inline bool position_equal(float* position) {
      for (unsigned int i = 0; i < 3; i++) {
	if (position_[i] != position[i]) return false;
      }
      return true;
    }
    
    inline bool read_dir_equal(float* cosines) {
      for (unsigned int i = 0; i < 3; i++) {
	if (read_dir_[i] != cosines[i]) return false;
      }
      return true;
    }
    
    inline bool phase_dir_equal(float* cosines) {
      for (unsigned int i = 0; i < 3; i++) {
	if (phase_dir_[i] != cosines[i]) return false;
      }
      return true;
    }

    inline bool slice_dir_equal(float* cosines) {
      for (unsigned int i = 0; i < 3; i++) {
	if (slice_dir_[i] != cosines[i]) return false;
      }
      return true;
    }

    bool reconfigure_;
    virtual void reconfigure();

    GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* duplicate_profile
      ( GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *profile );

    boost::shared_ptr< hoNDArray<float_complext> > extract_samples_from_queue
      ( ACE_Message_Queue<ACE_MT_SYNCH> *queue, bool acknowledge_sliding_window );

    // Compute trajectory/dcw for a reconstruction (to store internally)
    int calculate_trajectory_for_reconstruction(long profile_offset);
    int calculate_density_compensation_for_reconstruction();

    // Compute trajectory/dcw for adding frames to the accumulation buffer
    boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory_for_frame(long profile_offset);
    boost::shared_ptr< cuNDArray<float> > calculate_density_compensation_for_frame();

    // Compute trajectory/dcw for the fully sampled accumulation buffer (iterative buffer mode only)
    boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory_for_rhs(long profile_offset);
    boost::shared_ptr< cuNDArray<float> > calculate_density_compensation_for_rhs();

    int slices_;
    int sets_;
    int image_counter_;
    int image_series_;
    int device_number_;
    int mode_; // See note above

    long samples_per_profile_;
    long profiles_per_frame_;           // for an undersampled frame
    long frames_per_rotation_;          // representing a fully sampled frame
    long rotations_per_reconstruction_; // the number of rotations to batch per reconstruction. Set to '0' to reconstruct frames individually.
    long buffer_frames_per_rotation_;   // the number of buffer subcycles
    long buffer_length_in_rotations_;   // the number of buffer cycles

    long *previous_profile_;
    long *profiles_counter_frame_;
    long *profiles_counter_global_;

    long sliding_window_profiles_;
    long sliding_window_rotations_;

    float kernel_width_;
    float oversampling_factor_;

    unsigned int num_coils_;

    float position_[3];
    float read_dir_[3];
    float phase_dir_[3];
    float slice_dir_[3];

    bool output_timing_;
    bool buffer_using_solver_;
    bool buffer_update_needed_;

    boost::shared_ptr< hoNDArray<floatd2> > host_traj_recon_;
    boost::shared_ptr< hoNDArray<float> > host_weights_recon_;
    
    boost::shared_ptr< hoNDArray<float_complext> > csm_host_;
    boost::shared_ptr< hoNDArray<float_complext> > reg_host_;
    
    boost::shared_ptr< cuSenseBuffer<float,2> > acc_buffer_;

    std::vector<unsigned int> fov_;
    std::vector<unsigned int> image_dimensions_;
    std::vector<unsigned int> image_dimensions_recon_;
    uintd2 image_dimensions_recon_os_;

    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > frame_queue_;
    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > recon_queue_;
  };
}
