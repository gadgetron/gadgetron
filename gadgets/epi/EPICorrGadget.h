/** \file       EPICorrGadget.h                                                                                                      
    \brief      NYU-CBI's version of the EPI navigator correction                                                                    
    \author     Pablo Velasco, Souheil Inati
*/                                                                                                                                   

#ifndef EPICORRGADGET_H
#define EPICORRGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "hoArmadillo.h"
#include "gadgetron_epi_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Gadgetron{

  class  EXPORTGADGETS_EPI EPICorrGadget :
  public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      EPICorrGadget();
      virtual ~EPICorrGadget();

    protected:
      GADGET_PROPERTY(verboseMode, bool, "Verbose output", false);
      GADGET_PROPERTY(referenceNavigatorNumber, int, "Navigator number to be used as reference, both for phase correction and weig\
hts for filtering (default=1 -- second navigator)", 1);
      GADGET_PROPERTY(correctNavPhase, bool, "Correct the data for the phase of the reference navigator (default=false)", false);
      GADGET_PROPERTY(B0CorrectionMode, size_t, "B0 correction mode: 0=none, 1=mean b0 (default), 2=mean+linear term", 1);
      GADGET_PROPERTY(OEMagintudeCorrectionMode, size_t, "Odd-Even magnitude-correction mode: 0=none, 1=mean difference (default),\
 2=mean+linear term", 1);
      GADGET_PROPERTY(OEPhaseCorrectionMode, size_t, "Odd-Even phase-correction mode: 0=none, 1=mean phase difference, 2=mean+line\
ar term (default), 3=polynomial", 2);
      GADGET_PROPERTY(navigatorParameterFilterLength, size_t, "Number of repetitions to use to filter the navigator parameters (se\
t to 0 or negative for no filtering)", 5);
      GADGET_PROPERTY(navigatorParameterFilterExcludeVols, size_t, "Number of volumes/repetitions to exclude from the beginning of\
 the run when filtering the navigator parameters (default: 1)", 2);

      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
              GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      float filter_nav_correction_parameter( hoNDArray<float>& nav_corr_param_array,
                         hoNDArray<float>& weights_array,
                         size_t exc,  // current excitation number (for this set and slice)
                         size_t set,  // set of the array to filter (current one)
                         size_t slc,  // slice of the array to filter (current one)
                         size_t Nt,   // number of e2/timepoints/repetitions to filter
                         bool   filter_in_complex_domain = false);
      void increase_no_repetitions( size_t delta_rep );

      // in verbose mode, more info is printed out
      bool verboseMode_;


      // --------------------------------------------------
      // variables for navigator parameter computation
      // --------------------------------------------------

      int           Nx_;       // Length of the acquisition data along the readout direction
      uint32_t      RefNav_time_;             // Acquisition_time_stamp for the reference navigator (used for B0 correction)
      float         echoSpacing_ms_;          // Echo spacing in ms
      float         RefNav_to_Echo0_time_ES_; // Time (in echo-spacing uints) between the reference navigator and the first RO echo \
(used for B0 correction)
      arma::fvec    x_;        // Evenly spaced x-space locations
      arma::fmat    X_;        // Polynomial matrix: (1 x x^2 x^3...)

      arma::cx_fvec corrB0_;   // B0 correction
      arma::cx_fvec corrpos_;  // total correction for positive readouts
      arma::cx_fvec corrneg_;  //                  for negative readouts
      arma::cx_fcube navdata_; // 3D arma cube to store the navigator data

      // epi parameters
      int numNavigators_;
      int etl_;

      // for a given shot
      bool corrComputed_;
      int navNumber_;
      int epiEchoNumber_;
      bool startNegative_;


      // --------------------------------------------------
      // variables for navigator parameter filtering
      // --------------------------------------------------

      arma::fvec    t_;        // vector with repetition numbers, for navigator filtering
      size_t        E2_;       // number of kspace_encoding_step_2
      std::vector< std::vector<size_t> >  excitNo_;  // Excitation number (for each set and slice)

      // arrays for navigator parameter filtering:

      hoNDArray<float>    Nav_mag_;      // array to store the average navigator magnitude
      hoNDArray<float>    B0_slope_;     // array to store the B0-correction linear   term (for filtering)
      hoNDArray<float>    B0_intercept_; // array to store the B0-correction constant term (for filtering)
      hoNDArray<float>    OE_mag_slope_;     // array to store the Odd-Even magnitude-correction linear   term (for filtering)
      hoNDArray<float>    OE_mag_intercept_; // array to store the Odd-Even magnitude-correction constant term (for filtering)
      hoNDArray<float>    OE_mag_curvature_; // array to store the Odd-Even magnitude-correction quadratic term (for filtering)
      hoNDArray<float>    OE_phi_slope_;     // array to store the Odd-Even phase-correction linear   term (for filtering)
      hoNDArray<float>    OE_phi_intercept_; // array to store the Odd-Even phase-correction constant term (for filtering)
      std::vector< hoNDArray<float> > OE_phi_poly_coef_;   // vector of arrays to store the polynomial coefficients for Odd-Even pha\
se correction

    };
}
#endif //EPICORRGADGET_H

