#include "GadgetIsmrmrdReadWrite.h"
#include "EPICorrGadget.h"
#include "Gadgetron.h"

namespace Gadgetron{

  EPICorrGadget::EPICorrGadget() {}
  EPICorrGadget::~EPICorrGadget() {}

int EPICorrGadget::process_config(ACE_Message_Block* mb)
{
  boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));
  ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
  ISMRMRD::trajectoryDescriptionType traj_desc = (*e_seq.begin()).trajectoryDescription().get();

  if (std::strcmp(traj_desc.identifier().c_str(), "ConventionalEPI")) {
    GADGET_DEBUG1("Expected trajectory description identifier 'ConventionalEPI', not found.");
    return GADGET_FAIL;
  }

  for (ISMRMRD::trajectoryDescriptionType::userParameterLong_sequence::iterator i (traj_desc.userParameterLong().begin ()); i != traj_desc.userParameterLong().end(); ++i) {
    if (std::strcmp(i->name().c_str(),"numberOfNavigators") == 0) {
      numNavigators_ = i->value();
    } else if (std::strcmp(i->name().c_str(),"etl") == 0) {
      etl_ = i->value();
    }
  }

  verboseMode_ = this->get_bool_value("verboseMode");

  corrComputed_ = false;
  navNumber_ = -1;
  epiEchoNumber_ = -1;

  return 0;
}

int EPICorrGadget::process(
          GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  //std::cout << "Nav: " << navNumber_ << "    " << "Echo: " << epiEchoNumber_ << std::endl;

  // Get a reference to the acquisition header
  ISMRMRD::AcquisitionHeader &hdr = *m1->getObjectPtr();

  // Pass on the non-EPI data (e.g. FLASH Calibration)
  if (hdr.encoding_space_ref > 0) {
    // It is enough to put the first one, since they are linked
    if (this->next()->putq(m1) == -1) {
      m1->release();
      ACE_ERROR_RETURN( (LM_ERROR,
			 ACE_TEXT("%p\n"),
			 ACE_TEXT("EPICorrGadget::process, passing data on to next gadget")),
			-1);
    }
    return 0;
  }

  // We have data from encoding space 0.

  // Make an armadillo matrix of the data
  arma::cx_fmat adata = as_arma_matrix(m2->getObjectPtr());

  // Check to see if the data is a navigator line or an imaging line
  if (ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PHASECORR_DATA).isSet(hdr.flags)) {

    // Increment the navigator counter
    navNumber_ += 1;

    // If the number of navigators per shot is exceeded, then
    // we are at the beginning of the next shot
    if (navNumber_ == numNavigators_) {
      corrComputed_ = false;
      navNumber_ = 0;
      epiEchoNumber_ = -1;
    }
    
    // If we are at the beginning of a shot, then initialize
    if (navNumber_==0) {
      // Set the size of the corrections and storage arrays
      corrpos_.set_size( adata.n_rows);
      corrneg_.set_size( adata.n_rows );
      navdata_.set_size( adata.n_rows, hdr.active_channels, numNavigators_);
      // Store the first navigator's polarity
      startNegative_ = ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_REVERSE).isSet(hdr.flags);
    }

    // Store the navigator data
    navdata_.slice(navNumber_) = adata;

    // If this is the last of the navigators for this shot, then
    // compute the correction operator
    if (navNumber_ == (numNavigators_-1)) {
      arma::cx_fvec ctemp =  arma::zeros<arma::cx_fvec>(adata.n_rows);    // temp column complex
      arma::fvec tvec = arma::zeros<arma::fvec>(adata.n_rows);            // temp column real
      arma::fvec x = arma::linspace<arma::fvec>(-0.5, 0.5, adata.n_rows); // Evenly spaced x-space locations
      int p; // counter
      
      // Accumulate over navigator triplets and sum over coils
      // this is the average phase difference between odd and even navigators
      for (p=0; p<numNavigators_-2; p=p+2) {
	ctemp += arma::sum(arma::conj(navdata_.slice(p)+navdata_.slice(p+2)) % navdata_.slice(p+1),1);
      }
      
      // TODO: Add a configuration toggle to switch between correction types

      // Point-wise phase estimate
      //for (p=0; p<adata.n_rows; p++) {
      //  tvec[p] = std::arg(ctemp[p]);
      //}

      // Robust fit to a straight line
      float slope = ctemp.n_rows * std::arg(arma::cdot(ctemp.rows(0,ctemp.n_rows-2), ctemp.rows(1,ctemp.n_rows-1)));
      ctemp = ctemp % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(x.n_rows), -slope*x));
      float intercept = std::arg(arma::sum(ctemp));
      //std::cout << "Slope = " << slope << std::endl;
      //std::cout << "Intercept = " << intercept << std::endl;
      tvec = slope*x + intercept;
      
      // Odd and even phase corrections
      if (!startNegative_) {
	// if the first navigator is a positive readout, we need to flip the sign of our correction
	tvec = -1.0*tvec;
      }
      corrpos_ = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(x.n_rows), -0.5*tvec));
      corrneg_ = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(x.n_rows), +0.5*tvec));
      corrComputed_ = true;
    }

  }
  else {
    // Increment the echo number
    epiEchoNumber_ += 1;
    // TODO: use this to apply the B0 correction

    // Apply the correction
    // We use the armadillo notation that loops over all the columns
    if (ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_REVERSE).isSet(hdr.flags)) {
      // Negative readout
      for (int p=0; p<adata.n_cols; p++) {
	adata.col(p) %= corrneg_;
      }
      // Now that we have corrected we set the readout direction to positive
      hdr.flags &= !(ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_REVERSE).bitmask_);
    } 
    else {
      // Positive readout
      for (int p=0; p<adata.n_cols; p++) {
	adata.col(p) %= corrpos_;
      }
    }
  }

  // Pass on the imaging data
  // TODO: this should be controlled by a flag
  if (ISMRMRD::FlagBit(ISMRMRD::ACQ_IS_PHASECORR_DATA).isSet(hdr.flags)) {
    m1->release();
  } 
  else {
    // It is enough to put the first one, since they are linked
    if (this->next()->putq(m1) == -1) {
      m1->release();
      ACE_ERROR_RETURN( (LM_ERROR,
			 ACE_TEXT("%p\n"),
			 ACE_TEXT("EPICorrGadget::process, passing data on to next gadget")),
			-1);
    }
  }

  return 0;
}

GADGET_FACTORY_DECLARE(EPICorrGadget)
}


