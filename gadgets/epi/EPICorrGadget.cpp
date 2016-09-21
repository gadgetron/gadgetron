/** \file       EPICorrGadget.cpp
    \brief      NYU CBI's version of the EPI navigator correction
    \authors    Pablo Velasco and Souheil Inati
*/                                                                                                                                   

#include "EPICorrGadget.h"
#include "ismrmrd/xml.h"

namespace Gadgetron{

  #define OE_PHASE_CORR_POLY_ORDER 4

  EPICorrGadget::EPICorrGadget()
  {
    Nx_ = 0;
    echoSpacing_ms_ = 1.;   // some dummy value > 0
  }

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

  for (std::vector<ISMRMRD::UserParameterDouble>::iterator i (traj_desc.userParameterDouble.begin()); i != traj_desc.userParameterDouble.end(); ++i) {
    if (i->name == "echoSpacing") {
      echoSpacing_ms_ = i->value;
    }
  }
  GDEBUG_STREAM("ES_ms = " << echoSpacing_ms_);

  // Initialize arrays needed for temporal filtering, if filtering is requested:
  if (navigatorParameterFilterLength.value() > 1)
  {
    // TO DO: Take into account the acceleration along E2:
    E2_  = e_limits.kspace_encoding_step_2 ? e_limits.kspace_encoding_step_2->maximum - e_limits.kspace_encoding_step_2->minimum + 1 : 1;
    size_t REP = e_limits.repetition ? e_limits.repetition->maximum - e_limits.repetition->minimum + 1 : 1;
    size_t SET = e_limits.set        ? e_limits.set->maximum        - e_limits.set->minimum        + 1 : 1;
    size_t SLC = e_limits.slice      ? e_limits.slice->maximum      - e_limits.slice->minimum      + 1 : 1;
    // NOTE: For EPI sequences, "segment" indicates odd/even readout, so we don't need a separate dimension for it.
    //GDEBUG_STREAM("E2: " << E2_ << "; SLC: " << SLC << "; REP: " << REP << "; SET: " << SET);

    // For 3D sequences, the e2 index in the navigator is always 0 (there is no phase encoding in
    //   the navigator), so we keep track of the excitation number for each slice and set) to do
    //   the filtering:
    excitNo_.resize(SLC);
    for (size_t i = 0; i < SLC; ++i)
    {
      excitNo_[i].resize( SET, size_t(0) );
    }

    // For 3D sequences, all e2 phase encoding steps excite the whole volume, so the
    //   navigators should be the same.  So when we filter across repetitions, we have
    //   to do it also through e2.  Bottom line: e2 and repetition are equivalent.
    Nav_mag_.create(         E2_*REP, SET, SLC);
    B0_intercept_.create(    E2_*REP, SET, SLC);
    if (B0CorrectionMode.value() == 2)
    {
      B0_slope_.create(        E2_*REP, SET, SLC);
    }
    OE_mag_intercept_.create(E2_*REP, SET, SLC);
    if (OEMagintudeCorrectionMode.value() >= 2)
    {
      OE_mag_slope_.create(    E2_*REP, SET, SLC);
      if (OEMagintudeCorrectionMode.value() == 3)
      {
        OE_mag_curvature_.create(    E2_*REP, SET, SLC);
      }
    }
    OE_phi_intercept_.create(E2_*REP, SET, SLC);
    if (OEPhaseCorrectionMode.value() >= 2)
    {
      OE_phi_slope_.create(    E2_*REP, SET, SLC);
      if (OEPhaseCorrectionMode.value() == 3)
      {
        OE_phi_poly_coef_.resize( OE_PHASE_CORR_POLY_ORDER+1 );
        for (size_t i = 0; i < OE_phi_poly_coef_.size(); ++i)
        {
          OE_phi_poly_coef_[i].create( E2_*REP, SET, SLC);
        }
      }
    }

    // Armadillo vector of evenly-spaced timepoints to filter navigator parameters:
    t_ = arma::linspace<arma::fvec>( 0, navigatorParameterFilterLength.value()-1, navigatorParameterFilterLength.value() );
  }

  verboseMode_ = verboseMode.value();

  corrComputed_ = false;
  navNumber_ = -1;
  epiEchoNumber_ = -1;

  //GDEBUG_STREAM("EPICorrGadget configured");                                                                                         
  return 0;
}

int EPICorrGadget::process(
          GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  //GDEBUG_STREAM("Nav: " << navNumber_ << "    " << "Echo: " << epiEchoNumber_);

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
      // ... but only if the size along the RO direction has changed:
      if ( adata.n_rows != Nx_ ){
        // Set the size of the corrections and storage arrays
        Nx_ = adata.n_rows;
        corrpos_.set_size( Nx_ );
        corrneg_.set_size( Nx_ );
        corrB0_.set_size(  Nx_ );
        navdata_.set_size( Nx_ , hdr.active_channels, numNavigators_);
        navdata_.fill(0.);
      }
      // Store the first navigator's polarity
      startNegative_ = hdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_REVERSE);
    }

    if (navNumber_==referenceNavigatorNumber.value())
    {
      // get the time stamp for this first navigator:
      RefNav_time_ = hdr.acquisition_time_stamp;
      //GDEBUG_STREAM("RefNav acq time: " << RefNav_time_ );
    }

    // Store the navigator data
    navdata_.slice(navNumber_) = adata;

    // If this is the last of the navigators for this shot, then
    // compute the correction operator
    if (navNumber_ == (numNavigators_-1)) {
      arma::cx_fvec ctemp =  arma::zeros<arma::cx_fvec>( Nx_ );    // temp column complex
      arma::fvec tvec = arma::zeros<arma::fvec>( Nx_ );            // temp column real
      // if the size along the RO direction has changed, reset:
      if ( x_.n_rows != Nx_ ){
        x_ = arma::linspace<arma::fvec>(-0.5, 0.5, Nx_); // Evenly spaced x-space locations
        if ( OEMagintudeCorrectionMode.value() == 3 ||  OEPhaseCorrectionMode.value() == 3 )
        {
          X_  = arma::zeros<arma::fmat>( Nx_ ,OE_PHASE_CORR_POLY_ORDER+1);
          X_.col(0) = arma::ones<arma::fvec>( Nx_ );
          X_.col(1) = x_;                       // x
          X_.col(2) = arma::square(x_);         // x^2
          X_.col(3) = x_ % X_.col(2);           // x^3
          X_.col(4) = arma::square(X_.col(2));  // x^4
        }
      }

      // mean of the reference navigator (across RO and channels):
      std::complex<float> navMean = arma::mean( arma::vectorise( navdata_.slice(referenceNavigatorNumber.value()) ) );
      //GDEBUG_STREAM("navMean = " << navMean);

      // for clarity, we'll use the following when filtering navigator parameters:
      size_t set, slc, exc;
      if (navigatorParameterFilterLength.value() > 1)
      {
        set = hdr.idx.set;
        slc = hdr.idx.slice;
        // Careful: kspace_encode_step_2 for a navigator is always 0, and at this point we
        //          don't have access to the kspace_encode_step_2 for the next line.  Instead,
        //          keep track of the excitation number for this set and slice:
        exc = excitNo_[slc][set];   // excitation number with this same specific set and slc
        //GDEBUG_STREAM("Excitation number:" << exc << "; slice: " << slc);

        // If, for whatever reason, we are getting more repetitions than the
        //   header specified, increase the size of the array to accomodate:
        if ( exc >= (Nav_mag_.get_size(0)/E2_) )
        {
          this->increase_no_repetitions( 100 );     // add 100 volumes more, to be safe
        }
        Nav_mag_(exc, set, slc) = std::abs(navMean);
      }

      /////////////////////////////////////
      //////      B0 correction      //////
      /////////////////////////////////////

      if (B0CorrectionMode.value() > 0)
      {
        // Accumulate over navigator pairs and sum over coils
        // this is the average phase difference between consecutive odd or even navigators
        for (int p=0; p<numNavigators_-2; p=p+1)
        {
          ctemp += arma::sum(arma::conj(navdata_.slice(p)) % navdata_.slice(p+2),1);
        }
        // scale it, so that abs(ctemp) is close to the navigator signal (abs(navdata_.slice(p)))
        ctemp = ctemp / arma::sqrt(arma::abs(ctemp));

        float slope = 0.;
        float intercept = 0.;
        if (B0CorrectionMode.value() >= 1)
        {
          // If a linear term is requested, compute it first (in the complex domain):
          if (B0CorrectionMode.value() == 2)
          {            // Robust fit to a straight line:
            slope = (Nx_-1) * std::arg(arma::cdot(ctemp.rows(0,Nx_-2), ctemp.rows(1,Nx_-1)));
            //GDEBUG_STREAM("Slope = " << slope);
            // If we need to filter the estimate:
            if (navigatorParameterFilterLength.value() > 1)
            {
              // (Because to estimate the intercept (constant term) we need to use the slope estimate,
              //   we want to filter it first):
              //   - Store the value in the corresponding array (we want to store it before filtering)
              B0_slope_(exc, set, slc) = slope;
              //   - Filter parameter:
              slope = filter_nav_correction_parameter( B0_slope_, Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value() );
            }

            // Now correct for the slope, to be able to compute the average phase:
            ctemp = ctemp % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>( Nx_ ), -slope*x_));

          }   // end of the B0CorrectionMode.value()==2


          // Now, compute the mean phase:
          intercept = std::arg(arma::sum(ctemp));
          //GDEBUG_STREAM("Intercept = " << intercept);
          if (navigatorParameterFilterLength.value() > 1)
          {
            //   - Store the value found in the corresponding array:
            B0_intercept_(exc, set, slc) = intercept;
            //   - Filter parameters:
            // Filter in the complex domain (last arg:"true"), to avoid smoothing across phase wraps:
            intercept = filter_nav_correction_parameter( B0_intercept_, Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value(), true );
          }

          // Then, our estimate of the phase:
          tvec = slope*x_ + intercept;

        }       // end of B0CorrectionMode.value() >= 1

        // 0.5* because what we have calculated was the phase difference between every other navigator
        corrB0_ = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(Nx_), -0.5*tvec));

      }        // end of B0CorrectionMode.value() > 0
      else
      {      // No B0 correction:
        corrB0_.ones();
      }


      ////////////////////////////////////////////////////////
      //////      Odd-Even correction -- Magnitude      //////
      ////////////////////////////////////////////////////////

      arma::fvec tvec_mag( Nx_ );
      arma::fvec ctemp_odd( Nx_ );     // temp column complex for odd  magnitudes
      arma::fvec ctemp_even( Nx_ );    // temp column complex for even magnitudes
      if (OEMagintudeCorrectionMode.value() > 0)
      {                                                                                                                                
        // Accumulate over navigator triplets and sum over coils
        // this is the average amplitudes of odd and even navigators
        ctemp_odd.zeros();
        ctemp_even.zeros();
        for (int p=0; p<numNavigators_-2; p=p+2)
        {
          ctemp_odd  += ( arma::sqrt(arma::sum(arma::square(arma::abs(navdata_.slice(p))),1)) + arma::sqrt(arma::sum(arma::square(arma::abs(navdata_.slice(p+2))),1)) )/2;
          ctemp_even += arma::sqrt(arma::sum(arma::square(arma::abs(navdata_.slice(p+1))),1));
        }

        float intercept = 1.;
        float slope     = 0.;
        if (OEMagintudeCorrectionMode.value() == 1)
        {
          // Just mean term:
          float intercept = arma::as_scalar( arma::solve( ctemp_odd, ctemp_even ) );
          if (navigatorParameterFilterLength.value() > 1)
          {
            //   - Store the values found in the corresponding arrays:
            OE_mag_intercept_(exc, set, slc) = intercept;
            //   - Filter parameters:
            intercept = filter_nav_correction_parameter( OE_mag_intercept_, Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value() );
          }
          // Then, our estimate of the correction
          tvec_mag.fill(intercept);
        }
        else if (OEMagintudeCorrectionMode.value() == 2)
        {
          // Fit to a straight line + offset:  even/odd = ax + c = A*b:
          arma::fvec mag_coef = arma::solve( arma::join_horiz( ctemp_odd, ctemp_odd % x_ ),  ctemp_even );
          if (navigatorParameterFilterLength.value() > 1)
          {
            //   - Store the values found in the corresponding arrays:
            OE_mag_intercept_(exc, set, slc) = mag_coef(0);
            OE_mag_slope_(    exc, set, slc) = mag_coef(1);
            //   - Filter parameters:
            intercept = filter_nav_correction_parameter( OE_mag_intercept_, Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value() );
            slope = filter_nav_correction_parameter( OE_mag_slope_    , Nav_mag_,exc, set, slc, navigatorParameterFilterLength.value() );
          }
          else
          {
            intercept = mag_coef(0);
            slope     = mag_coef(1);
          }
          // Then, our estimate of the correction
          tvec_mag = intercept + slope*x_;
        }      // end of if OEMagintudeCorrectionMode.value() == 2
        else if (OEMagintudeCorrectionMode.value() == 3)
        {
          // Fit to a polynomial of degree 2:  even/odd = a0 + a1x + a2x^2 = X*A:
          arma::fvec mag_coef = arma::solve( arma::diagmat( ctemp_odd ) * X_.head_cols(3) ,  ctemp_even );
          if (navigatorParameterFilterLength.value() > 1)
          {
            //   - Store the values found in the corresponding arrays:
            OE_mag_intercept_(exc, set, slc) = mag_coef(0);
            OE_mag_slope_(    exc, set, slc) = mag_coef(1);
            OE_mag_curvature_(exc, set, slc) = mag_coef(2);
            //   - Filter parameters:
            mag_coef(0) = filter_nav_correction_parameter( OE_mag_intercept_, Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value() );
            mag_coef(1) = filter_nav_correction_parameter( OE_mag_slope_    , Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value() );
            mag_coef(2) = filter_nav_correction_parameter( OE_mag_curvature_, Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value() );
          }
          // Then, our estimate of the correction
          tvec_mag = X_.head_cols(3) * mag_coef;
        }      // end of if OEMagintudeCorrectionMode.value() == 3

      }        // end of OEMagintudeCorrectionMode > 0
      else
      {      // No OE magnitude correction:
        tvec_mag.ones();     // to turn it off, set to 1.
      }



      ////////////////////////////////////////////////////
      //////      Odd-Even correction -- Phase      //////
      ////////////////////////////////////////////////////

      if (OEPhaseCorrectionMode.value() > 0)
      {
        // Accumulate over navigator triplets and sum over coils
        // this is the average phase difference between odd and even navigators
        // Note: we have to correct for the B0 evolution between navigators before
        ctemp.zeros();      // set all elements to zero
        for (int p=0; p<numNavigators_-2; p=p+2)
        {
          ctemp += arma::sum( arma::conj( navdata_.slice(p)/repmat(corrB0_,1,navdata_.n_cols) + navdata_.slice(p+2)%repmat(corrB0_,1,navdata_.n_cols) ) % navdata_.slice(p+1),1);
        }
        // scale, so that abs(ctemp) is close to the navigator signal (abs(navdata_.slice(p)))
        ctemp = ctemp / arma::sqrt(arma::abs(ctemp));

        float slope = 0.;
        float intercept = 0.;
        if (OEPhaseCorrectionMode.value() >= 1)
        {
          // If a linear term is requested, compute it first (in the complex domain):
          if (OEPhaseCorrectionMode.value() >= 2)
          {            // Robust fit to a straight line:
            slope = (Nx_-1) * std::arg(arma::cdot(ctemp.rows(0,Nx_-2), ctemp.rows(1,Nx_-1)));
            // If we need to filter the estimate:
            if (navigatorParameterFilterLength.value() > 1)
            {
              // (Because to estimate the intercept (constant term) we need to use the slope estimate,
              //   we want to filter it first):
              //   - Store the value in the corresponding array (we want to store it before filtering)
              OE_phi_slope_(exc, set, slc) = slope;
              //   - Filter parameter:                                                                                                 
              slope = filter_nav_correction_parameter( OE_phi_slope_, Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value() );
            }

            // Now correct for the slope, to be able to compute the average phase:
            ctemp = ctemp % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>( Nx_ ), -slope*x_));

          }   // end of the OEPhaseCorrectionMode==2


          // Now, compute the mean phase:
          intercept = std::arg(arma::sum(ctemp));
          //GDEBUG_STREAM("Intercept = " << intercept);
          if (navigatorParameterFilterLength.value() > 1)
          {
            //   - Store the value found in the corresponding array:
            OE_phi_intercept_(exc, set, slc) = intercept;
            //   - Filter parameters:
            // Filter in the complex domain ("true"), to avoid smoothing across phase wraps:
            intercept = filter_nav_correction_parameter( OE_phi_intercept_, Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value(), true );
          }

          // Then, our estimate of the phase:
          tvec = slope*x_ + intercept;

          // If a polynomial fit is requested:
          if (OEPhaseCorrectionMode.value() == 3)
          {
            // Fit the residuals (i.e., after removing the linear trend) to a polynomial.
            //   Use the magnitude of the average odd navigator as weights.
            // You cannot fit the phase directly to the polynomial because it doesn't work
            //   in cases that the phase wraps across the image.
            // Since we have already removed the slope (in the if OEPhaseCorrectionMode >= 2 step),
            //   just remove the constant phase:
            ctemp = ctemp % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>( Nx_ ), -intercept*arma::ones<arma::fvec>( Nx_ )));

            arma::fmat WX     = arma::diagmat(ctemp_odd) * X_;
            arma::fvec Wctemp( Nx_ );
            for (int p=0; p<Nx_; p++)
            {
              Wctemp(p) = ctemp_odd(p) * std::arg(ctemp(p));
            }
            arma::fvec phase_poly_coef = arma::solve( WX , Wctemp );
            if (navigatorParameterFilterLength.value() > 1)
            {
              for (size_t i = 0; i < OE_phi_poly_coef_.size(); ++i)
              {
                //   - Store the value found in the corresponding array:
                OE_phi_poly_coef_[i](exc, set, slc) = phase_poly_coef(i);

                //   - Filter parameters:
                phase_poly_coef(i) = filter_nav_correction_parameter( OE_phi_poly_coef_[i], Nav_mag_, exc, set, slc, navigatorParameterFilterLength.value() );
              }
              //GDEBUG_STREAM("OE_phi_poly_coef size: " << OE_phi_poly_coef_.size());
            }

            // Then, update our estimate of the phase correction:
            tvec += X_ * phase_poly_coef;     // ( Note the "+=" )

          }   // end of OEPhaseCorrectionMode.value() == 3

        }       // end of OEPhaseCorrectionMode.value() >= 1


        // Point-wise phase estimate
        //for (int p=0; p<adata.n_rows; p++) {
        //  tvec[p] = std::arg(ctemp[p]);
        //}

      }    // end of OEPhaseCorrectionMode >=0
      else
      {      // No OEPhase correction:
        tvec.ones();
      }

      // Odd and even phase corrections
      if (!startNegative_)
      {
        // if the first navigator is a positive readout, we need to flip the sign of our phase correction
        tvec = -1.0*tvec;
        // not sure if this is correct!!  maybe we need to do:
        //        mag_coef = arma::solve( arma::join_horiz( ctemp_even, ctemp_even % x ),  ctemp_odd );
        //        tvec_mag = mag_coef(0) + mag_coef(1)*x;
        // and invert (element-wise) the magnitude correction
        tvec_mag = 1/tvec_mag;
      }

      // if (navigatorParameterFilterLength.value() > 1 &&
      //     slc == (Nav_mag_.get_size(2)-1) &&
      //     set == (Nav_mag_.get_size(1)-1) &&
      //     exc == (Nav_mag_.get_size(0)/E2_ - 1) )
      // {
      //   write_nd_array< float >( &Nav_mag_, "/tmp/navs/Nav_mag.real" );
      //   write_nd_array< float >( &B0_slope_, "/tmp/navs/B0_slope.real" );
      //   write_nd_array< float >( &B0_intercept_, "/tmp/navs/B0_intercept.real" );
      //   write_nd_array< float >( &OE_mag_slope_, "/tmp/navs/OE_mag_slope.real" );
      //   write_nd_array< float >( &OE_mag_intercept_, "/tmp/navs/OE_mag_intercept.real" );
      //   write_nd_array< float >( &OE_mag_curvature_, "/tmp/navs/OE_mag_curvature.real" );
      //   write_nd_array< float >( &OE_phi_slope_, "/tmp/navs/OE_phi_slope.real" );
      //   write_nd_array< float >( &OE_phi_intercept_, "/tmp/navs/OE_phi_intercept.real" );
      // }
  


      ////////////////////////////////////////
      //////      Total correction      //////
      ////////////////////////////////////////
  
      // Combine all corrections:
      if (correctNavPhase.value())
      {
        // include the mean shot-to-shot phase variation (important for 3D and segmented sequences):
        corrpos_ = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(Nx_), -0.5*tvec - std::arg(navMean) ));
        corrneg_ = tvec_mag % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(Nx_), +0.5*tvec - std::arg(navMean)));
      }
      else
      {
        corrpos_ = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(Nx_), -0.5*tvec ));
        corrneg_ = tvec_mag % arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>(Nx_), +0.5*tvec));
      }
  
      corrComputed_ = true;
  
      if (navigatorParameterFilterLength.value() > 1) {
        excitNo_[slc][set]++;
      }
  
    }    // end of  "if (navNumber_ == (numNavigators_-1))"                                                                          
  }      // end of  "if (hdr.isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PHASECORR_DATA))"                                                     
  else {
    // Not navigator data, but image data

    // Increment the echo number
    epiEchoNumber_ += 1;

    if (epiEchoNumber_ == 0)
    {
      // get the time-stamp for this first echo, compare it to the reference navigator, and
      //   save it in units of echo spacing, to take it into the account in the B0 correction:
      uint32_t acqTimeEcho0 = hdr.acquisition_time_stamp;
      // Note: the acquisition time stamps are in 2.5 ms units:
      //RefNav_to_Echo0_time_ES_ = 2.5 * float( acqTimeEcho0 - RefNav_time_ ) / echoSpacing_ms_ ;
      //GDEBUG_STREAM("RefNav_to_Echo0_time_ES_: " << RefNav_to_Echo0_time_ES_ << "; approx: " << numNavigators_-referenceNavigatorNumber.value() );
      // because the raster time for the acquisition time stamp is 2.5 ms, you can't accurately
      //    measure the time between the reference navigator and the first readout.  We'll use the
      //    number of navigators acquired after the reference navigator as an approximation:
      RefNav_to_Echo0_time_ES_ = numNavigators_-referenceNavigatorNumber.value();
      // Note: we could add "a little bit more", something like 0.1 (in units of ES), but it is
      //   going to depend on the Ky_max and the maximum gradient amplitude and slew rate availabe.
    }

    // Apply the correction
    // We are correcting for the mean phase of the second navigator (it is included in corrneg_
    //   and corrpos_), so we have to apply the B0 correction not from the first RO echo, but
    //   from the second navigator.  We can approximate the time between the second navigator and
    //   the first redout by the number of navigators that are collected after the second one:
    //   numNavigators_ - 2
    // This approximation doesn't take into account the EPI readout pre-winders, but it is probably
    //   not that long compared to the echo spacing.
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


////////////////////////////////////////////////////
//
//  filter_nav_correction_parameter
//
//    funtion to filter (over e2/repetition number) a navigator parameter.
//    - nav_corr_param_array: array of navigator parameters
//    - weights_array       : array with weights for the filtering
//    - exc                 : current excitation number (for this set and slice)
//    - set                 : set of the array to filter (current one)
//    - slc                 : slice of the array to filter (current one)
//    - Nt                  : number of e2/timepoints/repetitions to filter
//    - filter_in_complex_domain : whether to filter in the complex domain, to avoid +/- pi wraps (default: false)

  float EPICorrGadget::filter_nav_correction_parameter( hoNDArray<float>& nav_corr_param_array,
                                hoNDArray<float>& weights_array,
                                size_t exc,
                                size_t set,
                                size_t slc,
                                size_t Nt,
                                bool   filter_in_complex_domain )
  {
    // If the array to be filtered doesn't have 3 dimensions, we are in big trouble:
    if ( nav_corr_param_array.get_number_of_dimensions() != 3 )
    {
      GERROR("EPICorrGadget::filter_nav_correction_parameter, incorrect number of dimensions of the array.\n");
      return -1;
    }

    // The dimensions of the weights array should be the same as the parameter array:
    if ( !nav_corr_param_array.dimensions_equal( &weights_array ) )
    {
      GERROR("EPICorrGadget::filter_nav_correction_parameter, dimensions of the parameter and weights arrays don't match.\n");
      return -1;
    }

    // If this repetition number is less than then number of repetitions to exclude...
    if ( exc < navigatorParameterFilterExcludeVols.value()*E2_ )
    {
      //   no filtering is needed, just return the corresponding value:
      return nav_corr_param_array(exc, set, slc );
    }

    // for now, just to a simple (robust) linear fit to the previous Nt timepoints:
    // TO DO: do we want to do something fancier?

    //
    // extract the timeseries (e2 phase encoding steps and repetitions)
    // of parameters and weights corresponding to the requested indices:

    // make sure we don't use more timepoints (e2 phase encoding steps and repetitions)
    //    that the currently acquired (minus the ones we have been asked to exclude
    //    from the beginning of the run):
    Nt = std::min( Nt, exc - (navigatorParameterFilterExcludeVols.value()*E2_) + 1 );

    // create armadillo vectors, and stuff them in reverse order (from the
    // current timepoint, looking backwards). This way, the filtered value
    // we want would be simply the intercept):
    arma::fvec weights =  arma::zeros<arma::fvec>( Nt );
    arma::fvec params  =  arma::zeros<arma::fvec>( Nt );
    for (size_t t = 0; t < Nt; ++t)
    {
        weights(t) = weights_array(        exc-t, set, slc );
        params( t) = nav_corr_param_array( exc-t, set, slc );
    }

    /////     weighted fit:          b = (W*[1 t_])\(W*params);    /////

    float filtered_param;

    // if we need to filter in the complex domain:
    if (filter_in_complex_domain)
    {
        arma::cx_fvec zparams = arma::exp(arma::cx_fvec(arma::zeros<arma::fvec>( Nt ), params));            // zparams = exp( i*params );
        arma::cx_fvec B = arma::solve( arma::cx_fmat( arma::join_horiz( weights, weights % t_.head(Nt) ), arma::zeros<arma::fmat>( Nt,2) ),   weights % zparams );
        filtered_param = std::arg( arma::as_scalar(B(0)) );
    }
    else
    {
        arma::fvec B = arma::solve( arma::join_horiz( weights, weights % t_.head(Nt) ),  weights % params );
        filtered_param = arma::as_scalar(B(0));
    }

    //if ( current_timepoint==(weights_array.get_size(0)-1) && slc==(weights_array.get_size(2)-1) ){
    //    write_nd_array< float >( &weights_array, "/tmp/nav_weights.real" );
    //    write_nd_array< float >( &nav_corr_param_array, "/tmp/nav_param_array.real" );
    //    GDEBUG_STREAM("filtered parameter: " << filtered_param );
    //}

    return filtered_param;
  }


////////////////////////////////////////////////////
//
//  increase_no_repetitions
//
//    funtion to increase the size of the navigator parameter arrays used for filtering
//    - delta_rep: how many more repetitions to add

  void EPICorrGadget::increase_no_repetitions( size_t delta_rep )
  {
      size_t REP     = Nav_mag_.get_size(0)/E2_;   // current maximum number of repetitions
      size_t new_REP = REP + delta_rep;
      size_t SET     = Nav_mag_.get_size(1);
      size_t SLC     = Nav_mag_.get_size(2);

      // create a new temporary array:
      hoNDArray<float> tmpArray( E2_*new_REP, SET, SLC);
      tmpArray.fill(float(0.));

      // For each navigator parameter array, copy what we have so far to the temporary array, and then copy back:

      // Nav_mag_ :
      for (size_t slc = 0; slc < SLC; ++slc)
      {
          for (size_t set = 0; set < SET; ++set)
          {
              memcpy( &tmpArray(0,set,slc), &Nav_mag_(0,set,slc), Nav_mag_.get_number_of_bytes()/SET/SLC );
          }
      }
      Nav_mag_ = tmpArray;

      // B0_intercept_ :
      for (size_t slc = 0; slc < SLC; ++slc)
      {
          for (size_t set = 0; set < SET; ++set)
          {
              memcpy( &tmpArray(0,set,slc), &B0_intercept_(0,set,slc), B0_intercept_.get_number_of_bytes()/SET/SLC );
          }
      }
      B0_intercept_ = tmpArray;

      // B0_slope_ :
      if (B0CorrectionMode.value() == 2)
      {
          for (size_t slc = 0; slc < SLC; ++slc)
          {
              for (size_t set = 0; set < SET; ++set)
              {
                  memcpy( &tmpArray(0,set,slc), &B0_slope_(0,set,slc), B0_slope_.get_number_of_bytes()/SET/SLC );
              }
          }
          B0_slope_ = tmpArray;
      }

      // OE_mag_intercept_ :
      for (size_t slc = 0; slc < SLC; ++slc)
      {
          for (size_t set = 0; set < SET; ++set)
          {
              memcpy( &tmpArray(0,set,slc), &OE_mag_intercept_(0,set,slc), OE_mag_intercept_.get_number_of_bytes()/SET/SLC );
          }
      }
      OE_mag_intercept_ = tmpArray;

      // OE_mag_slope_ and OE_mag_curvature_ :
      if (OEMagintudeCorrectionMode.value() >= 2)
      {
          for (size_t slc = 0; slc < SLC; ++slc)
          {
              for (size_t set = 0; set < SET; ++set)
              {
                  memcpy( &tmpArray(0,set,slc), &OE_mag_slope_(0,set,slc), OE_mag_slope_.get_number_of_bytes()/SET/SLC );
              }
          }
          OE_mag_slope_ = tmpArray;

          if (OEMagintudeCorrectionMode.value() == 3)
          {
              for (size_t slc = 0; slc < SLC; ++slc)
              {
                  for (size_t set = 0; set < SET; ++set)
                  {
                      memcpy( &tmpArray(0,set,slc), &OE_mag_curvature_(0,set,slc), OE_mag_curvature_.get_number_of_bytes()/SET/SLC );
                  }
              }
              OE_mag_curvature_ = tmpArray;
          }
      }

      // OE_phi_intercept_ :
      for (size_t slc = 0; slc < SLC; ++slc)
      {
          for (size_t set = 0; set < SET; ++set)
          {
              memcpy( &tmpArray(0,set,slc), &OE_phi_intercept_(0,set,slc), OE_phi_intercept_.get_number_of_bytes()/SET/SLC );
          }
      }
      OE_phi_intercept_ = tmpArray;

      // OE_phi_slope_ :
      if (OEPhaseCorrectionMode.value() >= 2)
      {
          for (size_t slc = 0; slc < SLC; ++slc)
          {
              for (size_t set = 0; set < SET; ++set)
              {
                  memcpy( &tmpArray(0,set,slc), &OE_phi_slope_(0,set,slc), OE_phi_slope_.get_number_of_bytes()/SET/SLC );
              }
          }
          OE_phi_slope_ = tmpArray;

          // OE_phi_poly_coef_ :
          if (OEPhaseCorrectionMode.value() == 3)
          {
              for (size_t i = 0; i < OE_phi_poly_coef_.size(); ++i)
              {
                  for (size_t slc = 0; slc < SLC; ++slc)
                  {
                      for (size_t set = 0; set < SET; ++set)
                      {
                          memcpy( &tmpArray(0,set,slc), &OE_phi_poly_coef_[i](0,set,slc), OE_phi_poly_coef_[i].get_number_of_bytes()/SET/SLC );
                      }
                  }
                  OE_phi_poly_coef_[i] = tmpArray;
              }
          }
      }

      GDEBUG_STREAM("EPICorrGadget WARNING: repetition number larger than what specified in header");


  }



  GADGET_FACTORY_DECLARE(EPICorrGadget)
}
