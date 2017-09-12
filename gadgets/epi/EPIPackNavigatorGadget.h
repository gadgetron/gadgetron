#ifndef EPIPACKNAVIGATORGADGET_H
#define EPIPACKNAVIGATORGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "hoArmadillo.h"
#include "gadgetron_epi_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Gadgetron{

  class  EXPORTGADGETS_EPI EPIPackNavigatorGadget :
  public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      EPIPackNavigatorGadget();
      virtual ~EPIPackNavigatorGadget();

    protected:
      GADGET_PROPERTY(verboseMode, bool, "Verbose output", false);

      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
              GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      // in verbose mode, more info is printed out
      bool verboseMode_;

      arma::cx_fcube navdata_;

      // epi parameters
      int numNavigators_;

      // for a given shot
      int navNumber_;
      int epiEchoNumber_;

    };
}
#endif //EPIPACKNAVIGATORGADGET_H
