#include "EPICorrGadget.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

  EPICorrGadget::EPICorrGadget() {}
  EPICorrGadget::~EPICorrGadget() {}

int EPICorrGadget::process_config(ACE_Message_Block* mb)
{
  ISMRMRD::IsmrmrdHeader h;
  ISMRMRD::deserialize(mb->rd_ptr(),h);

  if (h.encoding.size() == 0) {
    GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
    GDEBUG("This Gadget needs an encoding description\n");
    return GADGET_FAIL;
  }

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


  for (std::vector<ISMRMRD::UserParameterLong>::iterator i (traj_desc.userParameterLong.begin()); i != traj_desc.userParameterLong.end(); ++i) {
    if (i->name == "numberOfNavigators") {
      numNavigators_ = i->value;
    } else if (i->name == "etl") {
      etl_ = i->value;
    }
  }

  verboseMode_ = verboseMode.value();

  corrComputed_ = false;
  navNumber_ = -1;
  epiEchoNumber_ = -1;

  return 0;
}

int EPICorrGadget::process(
          GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  //GDEBUG_STREAM("Nav: " << navNumber_ << "    " << "Echo: " << epiEchoNumber_ << std::endl);

  // Get a reference to the acquisition header
  ISMRMRD::AcquisitionHeader &hdr = *m1->getObjectPtr();

  // Pass on the non-EPI data (e.g. FLASH Calibration)
  if (hdr.encoding_space_ref > 0) {
    // It is enough to put the first one, since they are linked
    if (this->next()->putq(m1) == -1) {
      m1->release();
      GERROR("EPICorrGadget::process, passing data on to next gadget");
      return -1;
    }
    return 0;
  }

  // We have data from encoding space 0.

  // Make an armadillo matrix of the data
  arma::cx_fmat adata = as_arma_matrix(m2->getObjectPtr());

  // Check to see if the data is a navigator line or an imaging line
  if (hdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PHASECORR_DATA)) {

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
      corrB0_.set_size( adata.n_rows );
      corrpos_.set_size( adata.n_rows);
      corrneg_.set_size( adata.n_rows );
      navdata_.set_size( adata.n_rows, hdr.active_channels, numNavigators_);
      // Store the first navigator's polarity
      startNegative_ = hdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_REVERSE);
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
      
      
      /////////////////////////////////////                                                                             
      //////      B0 correction      //////                                                                             
      /////////////////////////////////////                                                                             

      if (B0CorrectionMode.value() > 0)
      {
          // Accumulate over navigator pairs and sum over coils                                                         
          // this is the average phase difference between consecutive odd or even navigators                            
          for (p=0; p<numNavigators_-2; p=p+1)
          {
              ctemp += arma::sum(arma::conj(navdata_.slice(p)) % navdata_.slice(p+2),1);
          }

          // Perform the fit:                                                                                           
          float slope = 0.;
          float intercept = 0.;
          if (B0CorrectionMode.value() >= 1)
          {
              // If a linear term is requested, compute it first (in the complex domain):                               
              if (B0CorrectionMode.value() == 2)
              {          // Robust fit to a straight line:                                                              
                  slope = (ctemp.n_rows-1) * std::arg(arma::cdot(ctemp.rows(0,ctemp.n_rows-2), ctemp.rows(1,ctemp.n_rows-1)));
                  //GDEBUG_STREAM("Slope = " << slope << std::endl);                                                    

                  // Correct for the slope, to be able to compute the average phase:                                    
                  ctemp = ctemp % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>( ctemp.n_rows ), -slope*x));
              }   // end of the B0CorrectionMode==2                                                                     

              // Now, compute the mean phase:                                                                           
              intercept = std::arg(arma::sum(ctemp));
              //GDEBUG_STREAM("Intercept = " << intercept << std::endl);                                                

              // Then, our estimate of the phase:                                                                       
              tvec = slope*x + intercept;

          }       // end of B0CorrectionMode >= 1                                                                       

          // The B0 Correction:                                                                                         
          // 0.5* because what we have calculated was the phase difference between every other navigator                
          corrB0_ = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(ctemp.n_rows), -0.5*tvec));

      }        // end of B0CorrectionMode > 0                                                                           
      else
      {      // No B0 correction:                                                                                       
          corrB0_.ones();
      }


      ////////////////////////////////////////////////////                                                              
      //////      Odd-Even correction -- Phase      //////                                                              
      ////////////////////////////////////////////////////                                                              

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
      float slope = (ctemp.n_rows-1) * std::arg(arma::cdot(ctemp.rows(0,ctemp.n_rows-2), ctemp.rows(1,ctemp.n_rows-1)));
      ctemp = ctemp % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(x.n_rows), -slope*x));
      float intercept = std::arg(arma::sum(ctemp));
      //GDEBUG_STREAM("Slope = " << slope << std::endl);
      //GDEBUG_STREAM("Intercept = " << intercept << std::endl);
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

    if (epiEchoNumber_ == 0)
    {
        // For now, we will correct the phase evolution of each EPI line, with respect
        //   to the first line in the EPI readout train (echo 0), due to B0 inhomogeneities.
        //   That is, the reconstructed images will have the phase that the object had at
        //   the beginning of the EPI readout train (excluding the phase due to encoding),
        //   multiplied by the coil phase.
        // Later, we could add the time between the excitation and echo 0, or between one
        //   of the navigators and echo 0, to correct for phase differences from shot to shot.
        //   This will be important for multi-shot EPI acquisitions.
        RefNav_to_Echo0_time_ES_ = 0;
    }

     // Apply the correction
    // We use the armadillo notation that loops over all the columns
    if (hdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_REVERSE)) {
      // Negative readout
      for (int p=0; p<adata.n_cols; p++) {
        adata.col(p) %= (arma::pow(corrB0_,epiEchoNumber_+RefNav_to_Echo0_time_ES_) % corrneg_);
      }
      // Now that we have corrected we set the readout direction to positive
      hdr.clearFlag(ISMRMRD::ISMRMRD_ACQ_IS_REVERSE);
    } 
    else {
      // Positive readout
      for (int p=0; p<adata.n_cols; p++) {
        adata.col(p) %= (arma::pow(corrB0_,epiEchoNumber_+RefNav_to_Echo0_time_ES_) % corrpos_);
      }
    }
  }

  // Pass on the imaging data
  // TODO: this should be controlled by a flag
  if (hdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PHASECORR_DATA)) {
    m1->release();
  } 
  else {
    // It is enough to put the first one, since they are linked
    if (this->next()->putq(m1) == -1) {
      m1->release();
      GERROR("EPICorrGadget::process, passing data on to next gadget");
      return -1;
    }
  }

  return 0;
}

GADGET_FACTORY_DECLARE(EPICorrGadget)
}


