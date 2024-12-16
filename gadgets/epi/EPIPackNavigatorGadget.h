#ifndef EPIPACKNAVIGATORGADGET_H
#define EPIPACKNAVIGATORGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "hoArmadillo.h"

#include <complex>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Gadgetron{

  class  EPIPackNavigatorGadget :
  public Gadget1<mrd::Acquisition>
    {
    public:
      EPIPackNavigatorGadget();
      virtual ~EPIPackNavigatorGadget();

    protected:
      GADGET_PROPERTY(verboseMode, bool, "Verbose output", false);

      virtual int process_config(const mrd::Header& header);
      virtual int process(GadgetContainerMessage<mrd::Acquisition>* m1);

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
