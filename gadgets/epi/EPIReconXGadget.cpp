#include "EPIReconXGadget.h"

#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

namespace Gadgetron{

  EPIReconXGadget::EPIReconXGadget() {}
  EPIReconXGadget::~EPIReconXGadget() {}

int EPIReconXGadget::process_config(const mrd::Header& header)
{
  auto&h = header;
  
  verboseMode_ = verboseMode.value();

  if (h.encoding.size() == 0) {
    GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
    GDEBUG("This Gadget needs an encoding description\n");
    return GADGET_FAIL;
  }

  GDEBUG("Number of encoding spaces = %d\n", h.encoding.size());

  // Get the encoding space and trajectory description
  mrd::EncodingSpaceType e_space = h.encoding[0].encoded_space;
  mrd::EncodingSpaceType r_space = h.encoding[0].recon_space;
  mrd::EncodingLimitsType e_limits = h.encoding[0].encoding_limits;
  mrd::TrajectoryDescriptionType traj_desc;

  if (h.encoding[0].trajectory_description) {
    traj_desc = *h.encoding[0].trajectory_description;
  } else {
    GDEBUG("Trajectory description missing");
    return GADGET_FAIL;
  }

  if (traj_desc.identifier != "ConventionalEPI") {
    GDEBUG("Expected trajectory description identifier 'ConventionalEPI', not found.");
    return GADGET_FAIL;
  }

  // Primary encoding space is for EPI
  reconx.encodeNx_  = e_space.matrix_size.x;
  reconx.encodeFOV_ = e_space.field_of_view_mm.x;
  reconx.reconNx_   = r_space.matrix_size.x;
  reconx.reconFOV_  = r_space.field_of_view_mm.x;
  
  // TODO: we need a flag that says it's a balanced readout.
  for (std::vector<mrd::UserParameterLongType>::iterator i (traj_desc.user_parameter_long.begin()); i != traj_desc.user_parameter_long.end(); ++i) {
    if (i->name == "rampUpTime") {
      reconx.rampUpTime_ = i->value;
    } else if (i->name == "rampDownTime") {
      reconx.rampDownTime_ = i->value;
    } else if (i->name == "flatTopTime") {
      reconx.flatTopTime_ = i->value;
    } else if (i->name == "acqDelayTime") {
      reconx.acqDelayTime_ = i->value;
    } else if (i->name == "numSamples") {
      reconx.numSamples_ = i->value;
    }
  }

  for (std::vector<mrd::UserParameterDoubleType>::iterator i (traj_desc.user_parameter_double.begin()); i != traj_desc.user_parameter_double.end(); ++i) {
    if (i->name == "dwellTime") {
      reconx.dwellTime_ = i->value;
    }
  }

  // If the flat top time is not set in the header, then we assume that rampSampling is off
  // and we set the flat top time from the number of samples and the dwell time.
  if (reconx.flatTopTime_ == 0) {
      reconx.flatTopTime_ = reconx.dwellTime_ * reconx.numSamples_;
  }

  // Compute the trajectory
  reconx.computeTrajectory();

  // Second encoding space is an even readout for PAT REF e.g. FLASH
  if ( h.encoding.size() > 1 ) {
    mrd::EncodingSpaceType e_space2 = h.encoding[1].encoded_space;
    mrd::EncodingSpaceType r_space2 = h.encoding[1].recon_space;
    reconx_other.encodeNx_  = r_space2.matrix_size.x;
    reconx_other.encodeFOV_ = r_space2.field_of_view_mm.x;
    reconx_other.reconNx_   = r_space2.matrix_size.x;
    reconx_other.reconFOV_  = r_space2.field_of_view_mm.x;
    reconx_other.numSamples_ = e_space2.matrix_size.x;
    oversamplng_ratio2_ = (float)e_space2.matrix_size.x / r_space2.matrix_size.x;
    reconx_other.dwellTime_ = 1.0;
    reconx_other.computeTrajectory();
  }

#ifdef USE_OMP
  omp_set_num_threads(1);
#endif // USE_OMP

  return 0;
}

int EPIReconXGadget::process(GadgetContainerMessage<mrd::Acquisition>* m1)
{
  auto& hdr_in = m1->getObjectPtr()->head;
  auto& data_in = m1->getObjectPtr()->data;

  mrd::AcquisitionHeader hdr_out;
  hoNDArray<std::complex<float> > data_out;

  data_out.create(reconx.reconNx_, data_in.get_size(1));

  // Switch the reconstruction based on the encoding space (e.g. for FLASH Calibration)
  if (hdr_in.encoding_space_ref == 0) {
    reconx.apply(hdr_in, data_in, hdr_out, data_out);
  }
  else
  {
    if(reconx_other.encodeNx_>data_in.get_size(0)/ oversamplng_ratio2_)
    {
        reconx_other.encodeNx_ = (int)(data_in.get_size(0) / oversamplng_ratio2_);
        reconx_other.computeTrajectory();
    }

    if (reconx_other.reconNx_>data_in.get_size(0) / oversamplng_ratio2_)
    {
        reconx_other.reconNx_ = (int)(data_in.get_size(0) / oversamplng_ratio2_);
    }

    if(reconx_other.numSamples_>data_in.get_size(0))
    {
        reconx_other.numSamples_ = data_in.get_size(0);
        reconx_other.computeTrajectory();
    }

    reconx_other.apply(hdr_in, data_in, hdr_out, data_out);
  }

  // Replace the contents of m1 with the new header and the contentes of m2 with the new data
  m1->getObjectPtr()->head = hdr_out;
  m1->getObjectPtr()->data = data_out;

  // It is enough to put the first one, since they are linked
  if (this->next()->putq(m1) == -1) {
    m1->release();
    GERROR("EPIReconXGadget::process, passing data on to next gadget");
    return -1;
  }

  return 0;
}

GADGET_FACTORY_DECLARE(EPIReconXGadget)
}


