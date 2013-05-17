#ifndef CUSBGOLDENRADIALGADGET_H
#define CUSBGOLDENRADIALGADGET_H
#pragma once

#include "cuSbGadget.h"

namespace Gadgetron{

  class EXPORTGADGETSSBSENSE cuSbGoldenRadialGadget : public cuSbGadget
  {

  public:
    GADGET_DECLARE(cuSbGoldenRadialGadget);

  protected:
    virtual boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory();
    virtual boost::shared_ptr< cuNDArray<float> > calculate_density_compensation();
  };
}

#endif //CUSBGOLDENRADIALGADGET_H
