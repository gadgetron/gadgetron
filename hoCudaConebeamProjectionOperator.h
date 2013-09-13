#pragma once
#include "hoCuNDArray_math.h"
#include "CBCT_acquisition.h"
#include "CBCT_binning.h"
#include "hoCuNDArray_operators.h"

#include "linearOperator.h"

#include <math_constants.h>
#include <vector>

namespace Gadgetron{
  
  class hoCudaConebeamProjectionOperator : public linearOperator< hoCuNDArray<float> >
  {
  public:
    hoCudaConebeamProjectionOperator() : linearOperator< hoCuNDArray<float> >()
    {
      samples_per_pixel_ = 1.5;      
      max_angle_ = 360.0f;
      use_fbp_ = false;
      use_oversampling_in_fbp_ = false;
      projections_per_batch_ = 20;
      preprocessed_ = false;
    }

    virtual ~hoCudaConebeamProjectionOperator() {}

    virtual void mult_M( hoCuNDArray<float> *in, hoCuNDArray<float> *out, bool accumulate = false );
    virtual void mult_MH( hoCuNDArray<float> *in, hoCuNDArray<float> *out, bool accumulate = false );

    virtual void setup( boost::shared_ptr<CBCT_acquisition> acquisition,
			boost::shared_ptr<CBCT_binning> binning,
			floatd3 is_dims_in_mm )
    {      
      acquisition_ = acquisition;
      binning_ = binning;
      is_dims_in_mm_ = is_dims_in_mm;
      preprocessed_ = true;
    }

    inline void use_filtered_backprojections( bool use_fbp ){
      use_fbp_ = use_fbp;      
    }

    inline void use_oversampling_in_filtered_backprojection( bool use_os ){
      use_oversampling_in_fbp_ = use_os;
    }

    inline void set_num_projections_per_batch( unsigned int projections_per_batch ){
      projections_per_batch_ = projections_per_batch;
    }

    inline void set_num_samples_per_pixel( float samples_per_pixel ){
      samples_per_pixel_ = samples_per_pixel;
    }

    inline void set_short_scan_maximum_angle( float angle ){
      max_angle_ = angle;
    }

    virtual boost::shared_ptr< linearOperator< hoCuNDArray<float> > > clone() {
      return linearOperator< hoCuNDArray<float> >::clone(this);
    }
    
  protected:
    boost::shared_ptr<CBCT_acquisition> acquisition_;
    boost::shared_ptr<CBCT_binning> binning_;
    floatd3 is_dims_in_mm_;
    float samples_per_pixel_;
    float max_angle_;
    bool use_fbp_;
    bool use_oversampling_in_fbp_;
    unsigned int projections_per_batch_;
    bool preprocessed_;
    boost::shared_ptr<hoCuNDArray<float> > variance_;
  };
}
