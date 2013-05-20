#pragma once

#include "hoNDArray.h"
#include "vector_td.h"

namespace Gadgetron{
  
  class SenseJob
  {
  public:
    
    SenseJob() {}
    ~SenseJob() {}
    
    boost::shared_ptr< hoNDArray<float_complext> > csm_host_;
    boost::shared_ptr< hoNDArray<float_complext> > reg_host_;
    boost::shared_ptr< hoNDArray<floatd2>        > tra_host_;
    boost::shared_ptr< hoNDArray<float>          > dcw_host_;
    boost::shared_ptr< hoNDArray<float_complext> > dat_host_;
  };
}
