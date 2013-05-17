#ifndef CUSBFIXEDRADIALGADGET_H
#define CUSBFIXEDRADIALGADGET_H
#pragma once

#include "cuSbGadget.h"

namespace Gadgetron{
  class EXPORTGADGETSSBSENSE cuSbFixedRadialGadget : public cuSbGadget
  {
  public:
    GADGET_DECLARE(cuSbFixedRadialGadget);

    cuSbFixedRadialGadget();

    virtual int process( GadgetContainerMessage< ISMRMRD::AcquisitionHeader >* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > > * m2 );
    virtual int process_config( ACE_Message_Block* mb );

  protected:
    virtual boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory();
    virtual boost::shared_ptr< cuNDArray<float> > calculate_density_compensation();

    unsigned int total_projections_;
    unsigned int dynamic_acceleration_factor_;
    int previous_projection_;
  };
}

#endif //CUSBGOLDENRADIALGADGET_H
