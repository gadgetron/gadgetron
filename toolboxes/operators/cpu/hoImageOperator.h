/** \file hoImageOperator.h
    \brief Image regularization operator, CPU based.
*/

#pragma once

#include "hoNDArray_math.h"
#include "complext.h"
#include "imageOperator.h"

#include <cmath>
#include <algorithm>

namespace Gadgetron{

  template <class T> class hoImageOperator : public imageOperator< hoNDArray<typename realType<T>::Type >, hoNDArray<T> >
  {
  public:

    hoImageOperator() : imageOperator< hoNDArray<typename realType<T>::Type >, hoNDArray<T> >() {}
    virtual ~hoImageOperator() {}    

    typedef typename imageOperator< hoNDArray<typename realType<T>::Type>, hoNDArray<T> >::REAL REAL;

    virtual boost::shared_ptr< linearOperator< hoNDArray<T> > > clone() {
      return linearOperator< hoNDArray<T> >::clone(this);
    }

  protected:

    // Estimate offset to the regularization image
    virtual REAL estimate_offset()
    {
      // Estimation based on simple histogram analysis:
      // Returns an estimation of the "average" intensity of the 'sigma' proportion of the image with the smallest intensities.
      //
      
      const unsigned int granularity = 50000; 
      std::vector<unsigned int> histogram(granularity,0);
      REAL max_value = this->image_->at(amax(this->image_.get()));
      REAL *d = this->image_->get_data_ptr();

      for( unsigned int i=0; i<this->image_->get_number_of_elements(); i++) {
	unsigned int bin = std::min(static_cast<unsigned int>(std::floor((d[i]/max_value)*granularity)), granularity-1);
	histogram[bin]++;
      }
      
      //Find 1th percentile
      //

      unsigned int cumsum = 0, counter = 0;
      while (cumsum < (unsigned int)(REAL(0.01)*this->image_->get_number_of_elements())) {
	cumsum += histogram[counter++];
      }      
      return REAL(counter+1)*max_value/granularity;
    }
  };
}
