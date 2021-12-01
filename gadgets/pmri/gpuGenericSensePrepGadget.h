#pragma once

#include "gadgetron_gpupmri_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"
#include "vector_td.h"
#include "cuNFFT.h"
#include "cuCgPreconditioner.h"
#include "cuSenseBufferCg.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <queue>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace Gadgetron{

  class EXPORTGADGETS_GPUPMRI gpuGenericSensePrepGadget :
    public Gadget3< ISMRMRD::AcquisitionHeader, hoNDArray< std::complex<float> >, hoNDArray<float> >
  {
    
  public:
    using ReadoutMessagePtr = GadgetContainerMessage<hoNDArray<std::complex<float>>> *;
    using TrajectoryMessagePtr = GadgetContainerMessage<hoNDArray<float>> *;
    using ImageHeaderMessagePtr = GadgetContainerMessage<ISMRMRD::ImageHeader> *;

    GADGET_DECLARE(gpuGenericSensePrepGadget);

    gpuGenericSensePrepGadget();
    virtual ~gpuGenericSensePrepGadget();

  protected:
    GADGET_PROPERTY(deviceno, int, "GPU device number", 0);
    GADGET_PROPERTY(buffer_length_in_rotations, int, "Number of rotations in a buffer", 1);
    GADGET_PROPERTY(buffer_using_solver, bool, "Use solver for buffer", false);
    GADGET_PROPERTY(buffer_convolution_kernel_width, float, "Convolution kernel width for buffer", 5.5);
    GADGET_PROPERTY(buffer_convolution_oversampling_factor, float, "Oversampling used in buffer convolution", 1.25);
    GADGET_PROPERTY(reconstruction_os_factor_x, float, "Oversampling for reconstruction in x-direction", 1.0);
    GADGET_PROPERTY(reconstruction_os_factor_y, float, "Oversampling for reconstruction in y-direction", 1.0);
    GADGET_PROPERTY(rotations_per_reconstruction, int, "Number of rotations per reconstruction", 0);
    GADGET_PROPERTY(propagate_csm_from_set, int, "Which set to use for CSM", -1);
    GADGET_PROPERTY(output_timing, bool, "Output timing information", false);
    GADGET_PROPERTY(sliding_window_readouts, int, "Number of readouts in sliding window", 0);
    GADGET_PROPERTY(sliding_window_rotations, int, "Number of rotations in sliding window", 0);
    GADGET_PROPERTY(readouts_per_frame, int, "Readouts per frame", 0);
    GADGET_PROPERTY(frames_per_rotation, int, "Frames per rotation", 0);
    GADGET_PROPERTY(buffer_frames_per_rotation, int, "Frames per rotation in buffer", 1);


    virtual int process_config(ACE_Message_Block *mb);

    virtual int process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader > *m1,        // header
			GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2,  // data
			GadgetContainerMessage< hoNDArray<float> > *m3 );                // traj/dcw

  private:

    inline bool vec_equal(float *in1, float *in2) {
      for (unsigned int i = 0; i < 3; i++) {
	if (in1[i] != in2[i]) return false;
      }
      return true;
    }
    
    boost::shared_array<bool> reconfigure_;
    virtual void reconfigure(unsigned int set, unsigned int slice);

    template<class T> GadgetContainerMessage< hoNDArray<T> >* 
      duplicate_array( GadgetContainerMessage< hoNDArray<T> > *array );
    
    boost::shared_ptr< hoNDArray<float_complext> > 
      extract_samples_from_queue ( std::queue<ReadoutMessagePtr> &queue,
				   bool sliding_window, unsigned int set, unsigned int slice );
    
    boost::shared_ptr< hoNDArray<float> > 
      extract_trajectory_from_queue ( std::queue<TrajectoryMessagePtr> &queue,
				      bool sliding_window, unsigned int set, unsigned int slice );
      
    void extract_trajectory_and_dcw_from_queue
      ( std::queue<TrajectoryMessagePtr> &queue, bool sliding_window, unsigned int set, unsigned int slice,
	unsigned int samples_per_frame, unsigned int num_frames,
	cuNDArray<floatd2> *traj, cuNDArray<float> *dcw );
    
    int slices_;
    int sets_;
    int device_number_;
    long samples_per_readout_;

    boost::shared_array<long> image_counter_;
    boost::shared_array<long> readouts_per_frame_;  // for an undersampled frame
    boost::shared_array<long> frames_per_rotation_; // representing a fully sampled frame

    // The number of rotations to batch per reconstruction. 
    // Set to '0' to reconstruct frames individually.
    long rotations_per_reconstruction_; 

    // The number of buffer cycles
    long buffer_length_in_rotations_; 

    boost::shared_array<long> buffer_frames_per_rotation_; // the number of buffer subcycles

    // Internal book-keping
    boost::shared_array<long> previous_readout_no_;
    boost::shared_array<long> acceleration_factor_;
    boost::shared_array<long> readout_counter_frame_;
    boost::shared_array<long> readout_counter_global_;

    long sliding_window_readouts_;
    long sliding_window_rotations_;

    float kernel_width_;
    float oversampling_factor_;

    boost::shared_array<unsigned int> num_coils_;

    boost::shared_array<float[3]> position_;
    boost::shared_array<float[3]> read_dir_;
    boost::shared_array<float[3]> phase_dir_;
    boost::shared_array<float[3]> slice_dir_;

    bool output_timing_;
    bool buffer_using_solver_;

    int propagate_csm_from_set_;
    boost::shared_ptr< cuNDArray<float_complext> > csm_;

    boost::shared_array<bool> buffer_update_needed_;

    boost::shared_array< hoNDArray<float_complext> > csm_host_;
    boost::shared_array< hoNDArray<float_complext> > reg_host_;
    
    boost::shared_array< cuSenseBuffer<float,2> > acc_buffer_;
    boost::shared_array< cuSenseBufferCg<float,2> > acc_buffer_cg_;

    std::vector<size_t> fov_;
    std::vector<size_t> image_dimensions_;
    std::vector<size_t> image_dimensions_recon_;
    uint64d2 image_dimensions_recon_os_;

    std::map<unsigned int, std::queue<ReadoutMessagePtr>> frame_readout_queue_;
    std::map<unsigned int, std::queue<ReadoutMessagePtr>> recon_readout_queue_;
    std::map<unsigned int, std::queue<TrajectoryMessagePtr>> frame_traj_queue_;
    std::map<unsigned int, std::queue<TrajectoryMessagePtr>> recon_traj_queue_;
    std::map<unsigned int, std::queue<ImageHeaderMessagePtr>> image_headers_queue_;
  };
}
