#include "EPIReconXGadget.h"
#include "ismrmrd/xml.h"

#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

namespace Gadgetron{

  EPIReconXGadget::EPIReconXGadget() {}
  EPIReconXGadget::~EPIReconXGadget() {}

int EPIReconXGadget::process_config(ACE_Message_Block* mb)
{
  ISMRMRD::IsmrmrdHeader h;
  ISMRMRD::deserialize(mb->rd_ptr(),h);
  
  
  verboseMode_ = verboseMode.value();

  if (h.encoding.size() == 0) {
    GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
    GDEBUG("This Gadget needs an encoding description\n");
    return GADGET_FAIL;
  }

  GDEBUG("Number of encoding spaces = %d\n", h.encoding.size());

  // Get the encoding space and trajectory description
  ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
  ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
  ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;
  ISMRMRD::TrajectoryDescription traj_desc;

  if (h.encoding[0].trajectoryDescription) {
    traj_desc = *h.encoding[0].trajectoryDescription;
  } else {
    GDEBUG("Trajectory description missing");
    return GADGET_FAIL;
  }

  if (traj_desc.identifier != "ConventionalEPI") {
    GDEBUG("Expected trajectory description identifier 'ConventionalEPI', not found.");
    return GADGET_FAIL;
  }

  // Primary encoding space is for EPI
  reconx.encodeNx_  = e_space.matrixSize.x;
  reconx.encodeFOV_ = e_space.fieldOfView_mm.x;
  reconx.reconNx_   = r_space.matrixSize.x;
  reconx.reconFOV_  = r_space.fieldOfView_mm.x;
  
  // TODO: we need a flag that says it's a balanced readout.
  for (std::vector<ISMRMRD::UserParameterLong>::iterator i (traj_desc.userParameterLong.begin()); i != traj_desc.userParameterLong.end(); ++i) {
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

  for (std::vector<ISMRMRD::UserParameterDouble>::iterator i (traj_desc.userParameterDouble.begin()); i != traj_desc.userParameterDouble.end(); ++i) {
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
    ISMRMRD::EncodingSpace e_space2 = h.encoding[1].encodedSpace;
    ISMRMRD::EncodingSpace r_space2 = h.encoding[1].reconSpace;
    reconx_other.encodeNx_  = r_space2.matrixSize.x;
    reconx_other.encodeFOV_ = r_space2.fieldOfView_mm.x;
    reconx_other.reconNx_   = r_space2.matrixSize.x;
    reconx_other.reconFOV_  = r_space2.fieldOfView_mm.x;
    reconx_other.numSamples_ = e_space2.matrixSize.x;
    oversamplng_ratio2_ = (float)e_space2.matrixSize.x / r_space2.matrixSize.x;
    reconx_other.dwellTime_ = 1.0;
    reconx_other.computeTrajectory();
  }

#ifdef USE_OMP
  omp_set_num_threads(1);
#endif // USE_OMP

  return 0;
}

int EPIReconXGadget::process(
          GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  ISMRMRD::AcquisitionHeader hdr_in = *(m1->getObjectPtr());
  ISMRMRD::AcquisitionHeader hdr_out;
  hoNDArray<std::complex<float> > data_out;

  data_out.create(reconx.reconNx_, m2->getObjectPtr()->get_size(1));

  // Switch the reconstruction based on the encoding space (e.g. for FLASH Calibration)
  if (hdr_in.encoding_space_ref == 0) {
    reconx.apply(*m1->getObjectPtr(), *m2->getObjectPtr(), hdr_out, data_out);
  }
  else
  {
    if(reconx_other.encodeNx_>m2->getObjectPtr()->get_size(0)/ oversamplng_ratio2_)
    {
        reconx_other.encodeNx_ = (int)(m2->getObjectPtr()->get_size(0) / oversamplng_ratio2_);
        reconx_other.computeTrajectory();
    }

    if (reconx_other.reconNx_>m2->getObjectPtr()->get_size(0) / oversamplng_ratio2_)
    {
        reconx_other.reconNx_ = (int)(m2->getObjectPtr()->get_size(0) / oversamplng_ratio2_);
    }

    if(reconx_other.numSamples_>m2->getObjectPtr()->get_size(0))
    {
        reconx_other.numSamples_ = m2->getObjectPtr()->get_size(0);
        reconx_other.computeTrajectory();
    }

    reconx_other.apply(*m1->getObjectPtr(), *m2->getObjectPtr(), hdr_out, data_out);
  }

  // Replace the contents of m1 with the new header and the contentes of m2 with the new data
  *m1->getObjectPtr() = hdr_out;
  *m2->getObjectPtr() = data_out;

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


