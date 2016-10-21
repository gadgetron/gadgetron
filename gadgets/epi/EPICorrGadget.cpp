#include "EPICorrGadget.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

  #define OE_PHASE_CORR_POLY_ORDER 4
  
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
    
    int Nx_ = adata.n_rows;

    // If we are at the beginning of a shot, then initialize
    if (navNumber_==0) {
      // Set the size of the corrections and storage arrays
      corrB0_.set_size(  Nx_ );
      corrpos_.set_size( Nx_ );
      corrneg_.set_size( Nx_ );
      navdata_.set_size( Nx_, hdr.active_channels, numNavigators_);
      navdata_.zeros();
      // Store the first navigator's polarity
      startNegative_ = hdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_REVERSE);
    }

    // Store the navigator data
    navdata_.slice(navNumber_) = adata;

    // If this is the last of the navigators for this shot, then
    // compute the correction operator
    if (navNumber_ == (numNavigators_-1)) {
      arma::cx_fvec ctemp =  arma::zeros<arma::cx_fvec>(Nx_);    // temp column complex
      arma::fvec tvec = arma::zeros<arma::fvec>(Nx_);            // temp column real
      arma::fvec x = arma::linspace<arma::fvec>(-0.5, 0.5, Nx_); // Evenly spaced x-space locations
      arma::fmat X;
      if ( OEPhaseCorrectionMode.value().compare("polynomial")==0 )
      {
          X  = arma::zeros<arma::fmat>( Nx_ ,OE_PHASE_CORR_POLY_ORDER+1);
          X.col(0) = arma::ones<arma::fvec>( Nx_ );
          X.col(1) = x;                       // x
          X.col(2) = arma::square(x);         // x^2
          X.col(3) = x % X.col(2);            // x^3
          X.col(4) = arma::square(X.col(2));  // x^4
      }
      int p; // counter
      
      
      /////////////////////////////////////
      //////      B0 correction      //////
      /////////////////////////////////////

      if ( B0CorrectionMode.value().compare("none")!=0 )    // If B0 correction is requested
      {
          // Accumulate over navigator pairs and sum over coils
          // this is the average phase difference between consecutive odd or even navigators
          for (p=0; p<numNavigators_-2; p++)
          {
              ctemp += arma::sum(arma::conj(navdata_.slice(p)) % navdata_.slice(p+2),1);
          }

          // Perform the fit:
          float slope = 0.;
          float intercept = 0.;
          if ( (B0CorrectionMode.value().compare("mean")==0)  ||
               (B0CorrectionMode.value().compare("linear")==0) )
          {
              // If a linear term is requested, compute it first (in the complex domain):
              if (B0CorrectionMode.value().compare("linear")==0)
              {          // Robust fit to a straight line:
                  slope = (Nx_-1) * std::arg(arma::cdot(ctemp.rows(0,Nx_-2), ctemp.rows(1,Nx_-1)));
                  //GDEBUG_STREAM("Slope = " << slope << std::endl);

                  // Correct for the slope, to be able to compute the average phase:
                  ctemp = ctemp % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>( Nx_ ), -slope*x));
              }   // end of the B0CorrectionMode == "linear"

              // Now, compute the mean phase:
              intercept = std::arg(arma::sum(ctemp));
              //GDEBUG_STREAM("Intercept = " << intercept << std::endl);

              // Then, our estimate of the phase:
              tvec = slope*x + intercept;

          }       // end of B0CorrectionMode == "mean" or "linear"

          // The B0 Correction:
          // 0.5* because what we have calculated was the phase difference between every other navigator
          corrB0_ = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(ctemp.n_rows), -0.5*tvec));

      }        // end of B0CorrectionMode != "none"
      else
      {      // No B0 correction:
          corrB0_.ones();
      }


      ////////////////////////////////////////////////////
      //////      Odd-Even correction -- Phase      //////
      ////////////////////////////////////////////////////

      if (OEPhaseCorrectionMode.value().compare("none")!=0)    // If Odd-Even phase correction is requested
      {
          // Accumulate over navigator triplets and sum over coils
          // this is the average phase difference between odd and even navigators
          // Note: we have to correct for the B0 evolution between navigators before
          ctemp.zeros();      // set all elements to zero
          for (p=0; p<numNavigators_-2; p=p+2)
          {
              ctemp += arma::sum( arma::conj( navdata_.slice(p)/repmat(corrB0_,1,navdata_.n_cols) + navdata_.slice(p+2)%repmat(corrB0_,1,navdata_.n_cols) ) % navdata_.slice(p+1),1);
          }

          float slope = 0.;
          float intercept = 0.;
          if ( (OEPhaseCorrectionMode.value().compare("mean")==0      ) ||
               (OEPhaseCorrectionMode.value().compare("linear")==0    ) ||
               (OEPhaseCorrectionMode.value().compare("polynomial")==0) )
          {
              // If a linear term is requested, compute it first (in the complex domain):                                                                 
              if ( (OEPhaseCorrectionMode.value().compare("linear")==0    ) ||
                   (OEPhaseCorrectionMode.value().compare("polynomial")==0) )
              {          // Robust fit to a straight line:                                                                                                
                  slope = (Nx_-1) * std::arg(arma::cdot(ctemp.rows(0,Nx_-2), ctemp.rows(1,Nx_-1)));

                  // Now correct for the slope, to be able to compute the average phase:
                  ctemp = ctemp % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>( Nx_ ), -slope*x));
              }   // end of the OEPhaseCorrectionMode == "linear" or "polynomial"

              // Now, compute the mean phase:
              intercept = std::arg(arma::sum(ctemp));

              // Then, our estimate of the phase:
              tvec = slope*x + intercept;

              // If a polynomial fit is requested:
              if (OEPhaseCorrectionMode.value().compare("polynomial")==0)
              {
                  // Fit the residuals (i.e., after removing the linear trend) to a polynomial.
                  // You cannot fit the phase directly to the polynomial because it doesn't work
                  //   in cases that the phase wraps across the image.
                  // Since we have already removed the slope (in the if OEPhaseCorrectionMode
                  //   == "linear" or "polynomial" step), just remove the constant phase:
                  ctemp = ctemp % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>( Nx_ ), -intercept*arma::ones<arma::fvec>( Nx_ )));

                  // Use the magnitude of the average odd navigator as weights:
                  arma::fvec ctemp_odd  = arma::zeros<arma::fvec>(Nx_);    // temp column complex for odd  magnitudes
                  for (int p=0; p<numNavigators_-2; p=p+2)
                  {
                      ctemp_odd  += ( arma::sqrt(arma::sum(arma::square(arma::abs(navdata_.slice(p))),1)) + arma::sqrt(arma::sum(arma::square(arma::abs(navdata_.slice(p+2))),1)) )/2;
                  }

                  arma::fmat WX     = arma::diagmat(ctemp_odd) * X;   // Weighted polynomial matrix
                  arma::fvec Wctemp( Nx_ );                           // Weighted phase residual
                  for (int p=0; p<Nx_; p++)
                  {
                      Wctemp(p) = ctemp_odd(p) * std::arg(ctemp(p));
                  }

                  // Solve for the polynomial coefficients:
                  arma::fvec phase_poly_coef = arma::solve( WX , Wctemp );

                  // Then, update our estimate of the phase correction:
                  tvec += X * phase_poly_coef;     // ( Note the "+=" )

              }   // end of OEPhaseCorrectionMode == "polynomial"

          }       // end of OEPhaseCorrectionMode == "mean", "linear" or "polynomial"

          if (!startNegative_) {
            // if the first navigator is a positive readout, we need to flip the sign of our correction
            tvec = -1.0*tvec;
          }
      }    // end of OEPhaseCorrectionMode != "none"
      else
      {      // No OEPhase correction:
          tvec.zeros();
      }

      // Odd and even phase corrections
      corrpos_ = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(Nx_), -0.5*tvec));
      corrneg_ = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(Nx_), +0.5*tvec));
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
