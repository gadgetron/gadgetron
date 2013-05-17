#include "cuSbFixedRadial.h"
#include "radial_utilities.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "Gadgetron.h"

namespace Gadgetron{

  cuSbFixedRadialGadget::cuSbFixedRadialGadget()
    : cuSbGadget()
    , previous_projection_(-1)
    , dynamic_acceleration_factor_(1)
  {
  }

  boost::shared_ptr< cuNDArray<floatd2> >
  cuSbFixedRadialGadget::calculate_trajectory()
  {
    // Define trajectories
    boost::shared_ptr< cuNDArray<floatd2> > traj = compute_radial_trajectory_fixed_angle_2d<float>
      ( samples_per_profile_, profiles_per_frame_, 1 );

    if (!traj.get()) {
      GADGET_DEBUG1("Failed to allocate trajectory");
      return boost::shared_ptr< cuNDArray<floatd2> >();
    }

    return traj;
  }

  int cuSbFixedRadialGadget::process(
				      GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
				      GadgetContainerMessage<hoNDArray<std::complex<float> > >* m2)
  {
    if (previous_projection_ < 0) {
      previous_projection_ = m1->getObjectPtr()->idx.kspace_encode_step_1;
    } else {
      if (m1->getObjectPtr()->idx.kspace_encode_step_1 > previous_projection_) {
	dynamic_acceleration_factor_ = m1->getObjectPtr()->idx.kspace_encode_step_1 - previous_projection_;
      }
      previous_projection_ = m1->getObjectPtr()->idx.kspace_encode_step_1;
    }

    profiles_per_frame_ = total_projections_/dynamic_acceleration_factor_;

    return cuSbGadget::process(m1,m2);
  }

  int cuSbFixedRadialGadget::process_config(ACE_Message_Block* mb)
  {
    int ret = cuSbGadget::process_config(mb);

    if (ret != GADGET_OK) {
      GADGET_DEBUG1("Base class process_config failed\n");
      return GADGET_FAIL;
    }

    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    std::vector<long> dims;
    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    if (e_seq.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    total_projections_ = e_limits.kspace_encoding_step_1().get().maximum()+1;

    return GADGET_OK;
  }


  boost::shared_ptr< cuNDArray<float> >
  cuSbFixedRadialGadget::calculate_density_compensation()
  {
    boost::shared_ptr< cuNDArray<float> > dcw = compute_radial_dcw_fixed_angle_2d
      ( samples_per_profile_, profiles_per_frame_, (float)matrix_size_os_.vec[0]/(float)matrix_size_.vec[0], 
	1.0f/((float)samples_per_profile_/(float)std::max(matrix_size_.vec[0],matrix_size_.vec[1])) );
  
    if (!dcw.get()) {
      GADGET_DEBUG1("Failed to calculate density compensation weights\n");
      return boost::shared_ptr< cuNDArray<float> >();
    }

    return dcw;
  }

  GADGET_FACTORY_DECLARE(cuSbFixedRadialGadget)
}
