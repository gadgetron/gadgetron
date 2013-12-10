#include "CartesianToGenericGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

namespace Gadgetron{

  CartesianToGenericGadget::CartesianToGenericGadget() 
  {
    set_parameter(std::string("matrix_size_as_a_multipluple_of").c_str(), "1");
  }

  CartesianToGenericGadget::~CartesianToGenericGadget() {}
  
  int CartesianToGenericGadget::process_config(ACE_Message_Block* mb)
  {
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    if( cfg.get() == 0x0 ){
      GADGET_DEBUG1("Unable to parse Ismrmrd header\n");
      return GADGET_FAIL;
    }

    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();

    if (e_seq.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    // Enforcement of the matrix size being a multiple of the "warp size"
    warp_size_ = get_int_value(std::string("matrix_size_as_a_multipluple_of").c_str());

    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    matrix_size_.push_back( (e_space.matrixSize().x()+warp_size_-1)/warp_size_*warp_size_);
    matrix_size_.push_back( (e_space.matrixSize().y()+warp_size_-1)/warp_size_*warp_size_);

    center_phase_ = e_limits.kspace_encoding_step_1().get().center();

    return GADGET_OK;
  }

  int CartesianToGenericGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2)
  {
    // Noise should have been consumed by the noise adjust, but just in case...
    //

    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
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
    cont->getObjectPtr()->create(&trajectory_dimensions);
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
      GADGET_DEBUG1("Failed to put job on queue.\n");
      return GADGET_FAIL;
    }
    
    return GADGET_OK;
  }
  
  GADGET_FACTORY_DECLARE(CartesianToGenericGadget)
}
