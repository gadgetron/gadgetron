#pragma once

#include "cuNDArray.h"
#include "gpureg_export.h"

namespace Gadgetron{
  
  // Downsample array to half size by averaging
  template<class REAL, unsigned int D> EXPORTGPUREG boost::shared_ptr< cuNDArray<REAL> > downsample( cuNDArray<REAL> *data );
  
  // Linear interpolation upsampling to array of doubled dimensions
  template<class REAL, unsigned int D> EXPORTGPUREG boost::shared_ptr< cuNDArray<REAL> > upsample( cuNDArray<REAL> *data );
}
