#include "GadgetIsmrmrdReadWrite.h"
#include "EPIReconXGadget.h"
#include "Gadgetron.h"

namespace Gadgetron{

int EPIReconXGadget::process_config(ACE_Message_Block* mb)
{
  boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));
  ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
  ISMRMRD::ismrmrdHeader::acquisitionSystemInformation_optional a_seq = cfg->acquisitionSystemInformation();

    verboseMode_ = this->get_bool_value("verboseMode");

    // Get the encoding space and trajectory description
    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();
    ISMRMRD::trajectoryDescriptionType traj_desc = (*e_seq.begin()).trajectoryDescription().get();

    if (std::strcmp(traj_desc.identifier().c_str(), "ConventionalEPI")) {
      GADGET_DEBUG1("Expected trajectory description identifier 'ConventionalEPI', not found.");
      return GADGET_FAIL;
    }
    
    reconx.encodeNx_  = e_space.matrixSize().x();
    reconx.encodeFOV_ = e_space.fieldOfView_mm().x();
    reconx.reconNx_   = r_space.matrixSize().x();
    reconx.reconFOV_  = r_space.fieldOfView_mm().x();

    // TODO: we need a flag that says it's a balanced readout.

    for (ISMRMRD::trajectoryDescriptionType::userParameterLong_sequence::iterator i (traj_desc.userParameterLong().begin ()); i != traj_desc.userParameterLong().end(); ++i) {
      if (std::strcmp(i->name().c_str(),"rampUpTime") == 0) {
	reconx.rampUpTime_ = i->value();
      } else if (std::strcmp(i->name().c_str(),"rampDownTime") == 0) {
	reconx.rampDownTime_ = i->value();
      } else if (std::strcmp(i->name().c_str(),"flatTopTime") == 0) {
	reconx.flatTopTime_ = i->value();
      } else if (std::strcmp(i->name().c_str(),"acqDelayTime") == 0) {
	reconx.acqDelayTime_ = i->value();
      } else if (std::strcmp(i->name().c_str(),"numSamples") == 0) {
	reconx.numSamples_ = i->value();
      } else {
	GADGET_DEBUG2("WARNING: unused trajectory parameter %s found\n", i->name().c_str());
      }
    }

    for (ISMRMRD::trajectoryDescriptionType::userParameterDouble_sequence::iterator i (traj_desc.userParameterDouble().begin ()); i != traj_desc.userParameterDouble().end(); ++i) {
      if (std::strcmp(i->name().c_str(),"dwellTime") == 0) {
	reconx.dwellTime_ = i->value();
      } else {
	GADGET_DEBUG2("WARNING: unused trajectory parameter %s found\n", i->name().c_str());
      }
    }

    // Compute the trajectory
    reconx.computeTrajectory();

  return 0;
}

int EPIReconXGadget::process(
          GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  ISMRMRD::AcquisitionHeader hdr_out;
  hoNDArray<std::complex<float> > data_out;
  data_out.create(reconx.reconNx_,m2->getObjectPtr()->get_size(1));

  reconx.apply(*m1->getObjectPtr(), *m2->getObjectPtr(), hdr_out, data_out);

  // Replace the contents of m1 with the new header and the contentes of m2 with the new data
  *m1->getObjectPtr() = hdr_out;
  *m2->getObjectPtr() = data_out;

  // It is enough to put the first one, since they are linked
  if (this->next()->putq(m1) == -1) {
    m1->release();
    ACE_ERROR_RETURN( (LM_ERROR,
		       ACE_TEXT("%p\n"),
		       ACE_TEXT("EPIReconXGadget::process, passing data on to next gadget")),
		      -1);
  }

  return 0;
}

GADGET_FACTORY_DECLARE(EPIReconXGadget)
}


