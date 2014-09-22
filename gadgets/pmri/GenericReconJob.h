#pragma once

#include "hoNDArray.h"
#include "vector_td.h"

#include <ismrmrd/ismrmrd.h>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

namespace Gadgetron{
  
  class GenericReconJob
  {
  public:
    
    GenericReconJob() {}
    ~GenericReconJob() {}

    boost::shared_array<ISMRMRD::ImageHeader> image_headers_;

    boost::shared_ptr< hoNDArray<float_complext> >  dat_host_;
    boost::shared_ptr< hoNDArray<floatd2>        >  tra_host_;
    boost::shared_ptr< hoNDArray<float>          >  dcw_host_;
    boost::shared_ptr< hoNDArray<float_complext> >  csm_host_;
    boost::shared_ptr< hoNDArray<float_complext> >  reg_host_;
  };
}
