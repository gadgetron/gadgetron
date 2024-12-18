#pragma once

#include "hoNDArray_math.h"
#include "resampleOperator.h"
#include "complext.h"
#include "hoArmadillo.h"

namespace Gadgetron{

  template <class T, unsigned int D>
  class hoLinearResampleOperator : public resampleOperator<hoNDArray<typename realType<T>::Type>, hoNDArray<T> >
  {
  public:

    hoLinearResampleOperator() : resampleOperator<hoNDArray<typename realType<T>::Type>, hoNDArray<T> >() {}
    virtual ~hoLinearResampleOperator() {}

    virtual void mult_M( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate = false);
    virtual void mult_MH( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate = false);
    virtual void set_displacement_field( boost::shared_ptr< hoNDArray<typename realType<T>::Type> > offsets );
    virtual void reset();

  private:
    inline bool is_border_pixel( typename reald<typename realType<T>::Type,D>::Type co, typename uint64d<D>::Type dims );
    inline unsigned int get_num_neighbors();

  protected:
    arma::SpMat<typename realType<T>::Type> R_T_; //Contains the TRANSPOSED resampling matrix.
  };
}
