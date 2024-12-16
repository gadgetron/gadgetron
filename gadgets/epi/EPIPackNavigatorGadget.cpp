#include "EPIPackNavigatorGadget.h"

namespace Gadgetron{

  EPIPackNavigatorGadget::EPIPackNavigatorGadget() {}
  EPIPackNavigatorGadget::~EPIPackNavigatorGadget() {}

int EPIPackNavigatorGadget::process_config(const mrd::Header& header)
{
  auto& h = header;
  if (h.encoding.size() == 0) {
    GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
    GDEBUG("This Gadget needs an encoding description\n");
    return GADGET_FAIL;
  }

  // Get the encoding space and trajectory description
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


  for (std::vector<mrd::UserParameterLongType>::iterator i (traj_desc.user_parameter_long.begin()); i != traj_desc.user_parameter_long.end(); ++i) {
    if (i->name == "numberOfNavigators") {
      numNavigators_ = i->value;
    }
  }

  verboseMode_ = verboseMode.value();

  navNumber_ = -1;
  epiEchoNumber_ = -1;

  return 0;
}

int EPIPackNavigatorGadget::process(GadgetContainerMessage<mrd::Acquisition>* m1)
{

  //GDEBUG_STREAM("Nav: " << navNumber_ << "    " << "Echo: " << epiEchoNumber_ << std::endl);

  // Get a reference to the acquisition header
  auto& hdr = m1->getObjectPtr()->head;

  // Pass on the non-EPI data (e.g. FLASH Calibration)
  if (hdr.encoding_space_ref > 0) {
    // It is enough to put the first one, since they are linked
    if (this->next()->putq(m1) == -1) {
      m1->release();
      GERROR("EPIPackNavigatorGadget::process, passing data on to next gadget");
      return -1;
    }
    return 0;
  }

  // We have data from encoding space 0.

  // Make an armadillo matrix of the data
  arma::cx_fmat adata = as_arma_matrix(m1->getObjectPtr()->data);

  // Check to see if the data is a navigator line or an imaging line
  if (hdr.flags.HasFlags(mrd::AcquisitionFlags::kIsPhasecorrData))
  {

    // Increment the navigator counter
    navNumber_ += 1;

    // If the number of navigators per shot is exceeded, then
    // we are at the beginning of the next shot
    if (navNumber_ == numNavigators_) {
      navNumber_ = 0;
      epiEchoNumber_ = -1;
    }

    // If we are at the beginning of a shot, then initialize
    if (navNumber_==0) {
      // Set the size of the storage array
      navdata_.set_size( adata.n_rows, m1->getObjectPtr()->Coils(), numNavigators_);
    }

    // Store the navigator data
    navdata_.slice(navNumber_) = adata;
  }
  else // imaging line
  {
      // Increment the echo number
      epiEchoNumber_ += 1;

      // Pack navigators into the integer component (c.f. msb) of early epi echo float data
      // I am assuming for early epi echoes that max(abs(adata)) < 0.5e-2
      // Matlab eventually unpacks with the following code
      //            navs=round(dat(:,:,1:3,:,:)); %% EPI data component < 0.5
      //            dat(:,:,1:3,:,:)=(dat(:,:,1:3,:,:)-navs)*1e-2;
      //            navs=navs*1e-6;

      if(epiEchoNumber_ < numNavigators_)
      {
	adata=1.0e2*adata+arma::round(1.0e6*navdata_.slice(epiEchoNumber_));
      }

  }

  // Pass on the imaging data
  // TODO: this should be controlled by a flag
  if (hdr.flags.HasFlags(mrd::AcquisitionFlags::kIsPhasecorrData)) {
    m1->release();
  }
  else {
    // It is enough to put the first one, since they are linked
    if (this->next()->putq(m1) == -1) {
      m1->release();
      GERROR("EPIPackNavigatorGadget::process, passing data on to next gadget");
      return -1;
    }
  }

  return 0;
}

GADGET_FACTORY_DECLARE(EPIPackNavigatorGadget)
}


