#pragma once

#include "hoNDArray.h"
#include "cpureg_export.h"

namespace Gadgetron{
  
  // Downsample array to half size by averaging
  template<class REAL, unsigned int D> EXPORTCPUREG boost::shared_ptr< hoNDArray<REAL> > downsample( hoNDArray<REAL> *data );
  
  // Linear interpolation upsampling to array of doubled dimensions
  template<class REAL, unsigned int D> EXPORTCPUREG boost::shared_ptr< hoNDArray<REAL> > upsample( hoNDArray<REAL> *data );
}
