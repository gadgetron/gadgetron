#pragma once

#include "linearOperator.h"
#include "hoCuNDArray.h"
#include <vector>
#include "GPUTimer.h"

#include "PS_Geometry.h"
#include "PS_BinningData.h"

namespace Gadgetron{
template<class REAL> class hoCudaConebeamProjectionOperator 
  : public linearOperator<hoCuNDArray<REAL> >
{
public:
  hoCudaConebeamProjectionOperator() : linearOperator<hoCuNDArray<REAL> >() {
    is_weights = 0x0;
    ps_weights = 0x0;
    geometry = 0x0;
    binning = 0x0;
    numProjections = 0x0;
    is_scale = 1.0;
    ps_scale = 1.0;
    preprocessed = false;
  }

  virtual ~hoCudaConebeamProjectionOperator() {}

  virtual void mult_M( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate = false );
  virtual void mult_MH( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate = false );
  virtual void mult_MH_M( hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, bool accumulate = false );

  virtual void setup(PS_Geometry* geometry,
                    PS_BinningData* binning,
                    std::vector<float>& angles, std::vector<float>& offsetx, std::vector<float>& offsety, unsigned int ppb,
                    floatd3 is_spacing_in_mm, uintd2 ps_dims_in_pixels,
                    unsigned int numSamplesPerRay,
                    bool use_fbp) { 

    this->geometry = geometry;
    this->binning = binning;
      
    this->angles = angles;
    this->offsetx = offsetx;
    this->offsety = offsety;
    this->ppb = ppb; // projections per batch
    this->is_spacing_in_mm = is_spacing_in_mm;
    this->ps_dims_in_pixels = ps_dims_in_pixels;

    this->use_fbp = use_fbp;
    this->numSamplesPerRay = numSamplesPerRay;

    this->preprocessed = true;
  }

  void set_is_scale(REAL scale) {
    is_scale = scale;
  }

  void set_ps_scale(REAL scale) {
    ps_scale = scale;
  }

  void set_is_weights(hoCuNDArray<REAL>* isw) {
    is_weights = isw;
  }

  void set_ps_weights(hoCuNDArray<REAL>* psw) {
    ps_weights = psw;
  }

  virtual boost::shared_ptr< linearOperator< hoCuNDArray<REAL> > > clone() {
    return linearOperator< hoCuNDArray<REAL> >::clone(this);
  }

protected:
  //std::vector<unsigned int> projections_per_bin;
  std::vector<float> angles;
  std::vector<float> offsetx;
  std::vector<float> offsety;
  bool use_fbp;
  unsigned int ppb;
  unsigned int numSamplesPerRay;
  floatd3 is_spacing_in_mm;
  uintd2 ps_dims_in_pixels;
  unsigned int numProjections;
  REAL is_scale, ps_scale;
  hoCuNDArray<REAL> *ps_weights, *is_weights;
  PS_Geometry* geometry;
  PS_BinningData* binning;
  bool preprocessed;
};
}
