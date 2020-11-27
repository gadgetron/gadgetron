#include "CartesianToGenericGadget.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

  CartesianToGenericGadget::CartesianToGenericGadget() 
  {
  }

  CartesianToGenericGadget::~CartesianToGenericGadget() {}
  
  int CartesianToGenericGadget::process_config(ACE_Message_Block* mb)
  {
    // Get the Ismrmrd header
    //
    ISMRMRD::IsmrmrdHeader h;
    ISMRMRD::deserialize(mb->rd_ptr(),h);
    
    
    if (h.encoding.size() != 1) {
      GDEBUG("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }
    
    // Get the encoding space and trajectory description
    ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
    ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
    ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

    // Enforcement of the matrix size being a multiple of the "warp size"
    warp_size_ = matrix_size_as_multiple_of.value();

    matrix_size_.push_back( (e_space.matrixSize.x+warp_size_-1)/warp_size_*warp_size_);
    matrix_size_.push_back( (e_space.matrixSize.y+warp_size_-1)/warp_size_*warp_size_);

    center_phase_ = e_limits.kspace_encoding_step_1 ? e_limits.kspace_encoding_step_1->center : 0;

    return GADGET_OK;
  }

  int CartesianToGenericGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2)
  {
    // Noise should have been consumed by the noise adjust, but just in case...
    //

    bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
    if (is_noise) {
      m1->release();
      return GADGET_OK;
    }

    // Make a new array as continuation of m1, and pass along
    //

    size_t samples_per_readout = m1->getObjectPtr()->number_of_samples;
    size_t center_sample = m1->getObjectPtr()->center_sample;
    size_t offset_readout = (matrix_size_[0]>>1)-center_sample; // In case of partial Fourier
    size_t offset_phase = (matrix_size_[1]>>1)-center_phase_; // In case of partial Fourier
    size_t phase_encode_step = m1->getObjectPtr()->idx.kspace_encode_step_1;

    std::vector<size_t> trajectory_dimensions;
    trajectory_dimensions.push_back(3);
    trajectory_dimensions.push_back(samples_per_readout);
    
    GadgetContainerMessage< hoNDArray<float> > *cont = new GadgetContainerMessage< hoNDArray<float> >();
    cont->getObjectPtr()->create(trajectory_dimensions);
    m2->cont(cont);

    float *traj_ptr = cont->getObjectPtr()->get_data_ptr();

    for( size_t sample=0; sample<samples_per_readout; sample++ ){

      // trajectory x (normalized to [-0.5;0.5])
      traj_ptr[sample*3+0] = float(sample+offset_readout)/float(matrix_size_[0])-0.5f;

      // trajectory y (normalized to [-0.5;0.5])
      traj_ptr[sample*3+1] = float(phase_encode_step+offset_phase)/float(matrix_size_[1])-0.5f;

      // dcw
      traj_ptr[sample*3+2] = 1.0f;
    }
        
    if (this->next()->putq(m1) < 0) {
      GDEBUG("Failed to put job on queue.\n");
      return GADGET_FAIL;
    }
    
    return GADGET_OK;
  }
  
  GADGET_FACTORY_DECLARE(CartesianToGenericGadget)
}
