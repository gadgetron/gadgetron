#pragma once

#include "CBCT_acquisition.h"
#include "CBCT_binning.h"
#include "hoNDArray_operators.h"
#include "linearOperator.h"

#include <math_constants.h>
#include <vector>

namespace Gadgetron{
  
  class hoCudaConebeamProjectionOperator : public linearOperator< hoNDArray<float> >
  {
  public:
    hoCudaConebeamProjectionOperator() : linearOperator< hoNDArray<float> >() 
    {
      preprocessed = false;
    }

    virtual ~hoCudaConebeamProjectionOperator() {}

    virtual void mult_M( hoNDArray<float> *in, hoNDArray<float> *out, bool accumulate = false );
    virtual void mult_MH( hoNDArray<float> *in, hoNDArray<float> *out, bool accumulate = false );
    //virtual void mult_MH_M( hoNDArray<float> *in, hoNDArray<float> *out, bool accumulate = false );

    virtual void setup( boost::shared_ptr<CBCT_acquisition> acquisition,
			boost::shared_ptr<CBCT_binning> binning,
			unsigned int projections_per_batch,
			unsigned int num_samples_per_ray,
			floatd3 is_spacing_in_mm,
			bool use_fbp = false,
			bool use_oversampling_in_fbp = false, 
			float maximum_angle = 2.0f*CUDART_PI_F ) 
    {      
      this->acquisition = acquisition;
      this->binning = binning;      
      this->projections_per_batch = projections_per_batch;
      this->num_samples_per_ray = num_samples_per_ray;
      this->is_spacing_in_mm = is_spacing_in_mm;
      this->use_fbp = use_fbp;
      this->use_oversampling_in_fbp = use_oversampling_in_fbp;
      this->maximum_angle = maximum_angle;
      this->preprocessed = true;
    }

    virtual boost::shared_ptr< linearOperator< hoNDArray<float> > > clone() {
      return linearOperator< hoNDArray<float> >::clone(this);
    }

  protected:
    boost::shared_ptr<CBCT_acquisition> acquisition;
    boost::shared_ptr<CBCT_binning> binning;
    unsigned int projections_per_batch;
    unsigned int num_samples_per_ray;
    floatd3 is_spacing_in_mm;
    bool use_fbp;
    bool use_oversampling_in_fbp;
    float maximum_angle;
    bool preprocessed;
  };
}
