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
      GADGET_PROPERTY(B0CorrectionMode, size_t, "B0 correction mode: 0=none, 1=mean b0 (default), 2=mean+linear term", 1);

      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
              GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      // in verbose mode, more info is printed out
      bool verboseMode_;

      // --------------------------------------------------                                                             
      // variables for navigator parameter computation                                                                  
      // --------------------------------------------------                                                             

      float         RefNav_to_Echo0_time_ES_; // Time (in echo-spacing uints) between the reference navigator and the f\
irst RO echo (used for B0 correction)                                                                                   
      arma::cx_fvec corrB0_;      // B0 correction                                                                       
      arma::cx_fvec corrpos_;     // Odd-Even correction -- positive readouts
      arma::cx_fvec corrneg_;     // Odd-Even correction -- negative readouts
      arma::cx_fcube navdata_;

      // epi parameters
      int numNavigators_;
      int etl_;

      // for a given shot
      bool corrComputed_;
      int navNumber_;
      int epiEchoNumber_;
      bool startNegative_;

    };
}
#endif //EPICORRGADGET_H
