/** \file b1_map.h
    \brief Utility to estimate b1 maps (MRI coil sensitivities), GPU based. 
*/

#pragma once

#include "gpupmri_export.h"
#include "cuNDArray.h"
#include "vector_td.h"
#include "complext.h"

#include <boost/smart_ptr.hpp>

namespace Gadgetron{

  /** 
   * \brief Estimate b1 map (coil sensitivities) of single or double precision according to REAL and of dimensionality D.
   * \param data Reconstructed reference images from the individual coils. Dimensionality is D+1 where the latter dimensions denotes the coil images.
   * \param taget_coils Denotes the number of target coils. Cannot exceed the size of dimension D of the data. A negative value indicates that sensitivity maps are computed for the full coil image dimension.
   */
  template<class REAL, unsigned long long D> EXPORTGPUPMRI boost::shared_ptr< cuNDArray<complext<REAL> > >
  estimate_b1_map( cuNDArray<complext<REAL> > *data, int target_coils = -1 );
}
