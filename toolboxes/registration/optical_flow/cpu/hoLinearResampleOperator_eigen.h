#pragma once

#include "hoNDArray_math.h"

#include "resampleOperator.h"
#include "complext.h"
#include "cpureg_export.h"

#include <armadillo>
#include <Eigen/Sparse>

namespace Gadgetron{

  template <class T, unsigned int D>
  class EXPORTCPUREG hoLinearResampleOperator_eigen : public resampleOperator<hoNDArray<typename realType<T>::Type>, hoNDArray<T> >
  {  
  public:
  
    hoLinearResampleOperator_eigen() : resampleOperator<hoNDArray<typename realType<T>::Type>, hoNDArray<T> >() {}
    virtual ~hoLinearResampleOperator_eigen() {}
  
    virtual void mult_M( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate = false);
    virtual void mult_MH( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate = false);
    virtual void set_displacement_field( boost::shared_ptr< hoNDArray<typename realType<T>::Type> > offsets );
  
    virtual unsigned int get_temporal_dimension_size() { return temporal_dim_size_; }
  
    virtual boost::shared_ptr< linearOperator< hoNDArray<T> > > clone() {
      return linearOperator< hoNDArray<T> >::clone(this);
    }
  
  private:
    inline bool is_border_pixel( typename reald<typename realType<T>::Type,D>::Type co, typename uintd<D>::Type dims );
    inline unsigned int get_num_neighbors();
  
  protected:
    boost::shared_ptr< Eigen::SparseMatrix<typename realType<T>::Type> > R_;
    unsigned int temporal_dim_size_;
  };
}
