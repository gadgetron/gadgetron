#pragma once

#include "gadgetron_hyper_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "cuNDArray.h"
#include "complext.h"
#include "cuNFFTOperator.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>

namespace Gadgetron {

  class EXPORTGADGETSHYPER CMRT3DGadget : 
    public Gadget2< ISMRMRD::ImageHeader, hoNDArray< std::complex<float> > >
  {  
  public:
    CMRT3DGadget() : images_received_(0), images_used_(0) {};
    ~CMRT3DGadget() {};

  protected:

    virtual int process_config(ACE_Message_Block* mb);
    virtual int process(GadgetContainerMessage< ISMRMRD::ImageHeader > *m1,
                        GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2);

    virtual boost::shared_ptr< cuNDArray<floatd2> > calculate_trajectory();
    virtual boost::shared_ptr< cuNDArray<float> > calculate_density_compensation();

    boost::shared_ptr< cuNDArray< complext<float> > > buffer_;
    boost::shared_ptr< cuNFFTOperator<float,2> > E_;
    std::vector<size_t> image_space_dimensions_3D_;
    unsigned int num_projections_expected_;
    unsigned int num_projections_to_use_;
    unsigned int projections_percentage_;
    unsigned int images_received_;
    unsigned int images_used_;
  };
}
