/** \file b1_map.h
    \brief Utility to estimate b1 maps (MRI coil sensitivities), GPU based. 
*/

#pragma once

#include "gpupmri_export.h"
#include "cuNDArray.h"
#include "vector_td.h"
#include "complext.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace Gadgetron{

  /** 
   * \brief Estimate b1 map (coil sensitivities) of single or double precision according to REAL and of dimensionality D.
   * \param data Reconstructed reference images from the individual coils. Dimensionality is D+1 where the latter dimensions denotes the coil images.
   * \param taget_coils Denotes the number of target coils. Cannot exceed the size of dimension D of the data. A negative value indicates that sensitivity maps are computed for the full coil image dimension.
   */
  template<class REAL, unsigned int D> EXPORTGPUPMRI  cuNDArray<complext<REAL> >
  estimate_b1_map(const cuNDArray<complext<REAL>>& data, int target_coils = -1 );

    /** 
   * \brief Estimate b1 map (coil sensitivities) of single or double precision using the NIH Souheil method
   * \param data [RO E1 CHA] for single 2D or [RO E1 N CHA] for multiple 2D reconstructed reference images from the individual coils. 
   */
  template<class REAL> EXPORTGPUPMRI bool
  estimate_b1_map_2D_NIH_Souheil( cuNDArray<complext<REAL> >* data, cuNDArray<complext<REAL> >* csm, size_t ks, size_t power,
                                  cuNDArray<complext<REAL> >& D, cuNDArray<complext<REAL> >& DH_D, 
                                  cuNDArray<complext<REAL> >& V1, cuNDArray<complext<REAL> >& U1 );
}
