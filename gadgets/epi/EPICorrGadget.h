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
      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
              GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      // in verbose mode, more info is printed out
      bool verboseMode_;

      arma::cx_fvec corrpos_;
      arma::cx_fvec corrneg_;  
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
