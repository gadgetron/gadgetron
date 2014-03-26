#pragma once

#include "gadgetron_radial_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "vector_td.h"
#include "cuNFFT.h"
#include "cuCgPreconditioner.h"
#include "cuSenseBufferCg.h"

#include <ismrmrd.h>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

/*
  Prep gadget for retrospectively gated Sense based on golden ratio sampling.
  Thus only radial modes 2-3 are supported.  
*/

namespace Gadgetron{

  class EXPORTGADGETS_RADIAL gpuRetroGatedSensePrepGadget :
    public Gadget2< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> > >
  {

  public:

    gpuRetroGatedSensePrepGadget();
    virtual ~gpuRetroGatedSensePrepGadget();

  protected:
    
    virtual int process_config(ACE_Message_Block *mb);

    virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader > *m1,
			GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2);

  private:

    inline bool vec_equal(float *in1, float *in2) {
      for (unsigned int i = 0; i < 3; i++) {
        if (in1[i] != in2[i]) return false;
      }
      return true;
    }
    
    boost::shared_array<bool> reconfigure_;
    virtual void reconfigure(unsigned int set, unsigned int slice);

    GadgetContainerMessage< hoNDArray< std::complex<float> > >*
      duplicate_profile( GadgetContainerMessage< hoNDArray< std::complex<float> > > *profile );

    boost::shared_ptr< hoNDArray<float_complext> > extract_samples_from_buffer_queue( unsigned int set, unsigned int slice );

    int extract_samples_and_trajectory_from_recon_queue
      ( unsigned int set, unsigned int slice, boost::shared_ptr< hoNDArray<float_complext> > samples, boost::shared_ptr< hoNDArray<floatd2> > trajectory );

    int calculate_density_compensation_for_reconstruction(unsigned int set, unsigned int slice);

    boost::shared_ptr< cuNDArray<floatd2> > 
      calculate_trajectory_for_buffer(long profile_offset, unsigned int set, unsigned int slice);

    boost::shared_ptr< cuNDArray<float> >
      calculate_density_compensation_for_buffer(unsigned int set, unsigned int slice);

    boost::shared_ptr< cuNDArray<floatd2> > 
      calculate_trajectory_for_rhs(long profile_offset, unsigned int set, unsigned int slice);

    boost::shared_ptr< cuNDArray<float> > 
      calculate_density_compensation_for_rhs(unsigned int set, unsigned int slice);

    int slices_;
    int sets_;
    int device_number_;
    int mode_;

    unsigned short phys_time_index_;

    long samples_per_profile_;
    long profiles_per_frame_;
    long frames_per_cardiac_cycle_;

    // The number of buffer cycles
    long profiles_per_buffer_frame_;
    long num_buffer_frames_inner_; 
    long num_buffer_frames_outer_;

    // Internal book-keeping
    boost::shared_array<unsigned int> first_profile_acq_time_;
    boost::shared_array<unsigned int> first_profile_phys_time_;
    boost::shared_array<unsigned int> previous_timestamp_;
    boost::shared_array<long> profiles_counter_global_;

    // We will discard profiles until the first R-wave is encountered
    boost::shared_array<bool> Rw_reached_;
    boost::shared_array<unsigned int> Rw_offset_;

    // For the buffer
    float kernel_width_;
    float oversampling_factor_;

    boost::shared_array<long> image_counter_;
    boost::shared_array<unsigned int> num_coils_;

    boost::shared_array<float[3]> position_;
    boost::shared_array<float[3]> read_dir_;
    boost::shared_array<float[3]> phase_dir_;
    boost::shared_array<float[3]> slice_dir_;

    bool output_timing_;
    bool buffer_using_solver_;

    boost::shared_array<bool> buffer_update_needed_;

    boost::shared_array< hoNDArray<float> > host_weights_recon_;
    
    boost::shared_array< hoNDArray<float_complext> > csm_host_;
    boost::shared_array< hoNDArray<float_complext> > reg_host_;
    
    boost::shared_array< cuSenseBuffer<float,2> > acc_buffer_;
    boost::shared_array< cuSenseBufferCg<float,2> > acc_buffer_cg_;

    std::vector<size_t> fov_;
    std::vector<size_t> image_dimensions_;
    std::vector<size_t> image_dimensions_recon_;
    uint64d2 image_dimensions_recon_os_;

    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > buffer_profiles_queue_;
    boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> > recon_profiles_queue_;
  };
}
