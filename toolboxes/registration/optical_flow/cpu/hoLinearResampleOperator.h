#pragma once

#include "hoResampleOperator.h"
#include <Eigen/Sparse>

template <class REAL, class T, unsigned int D>
class EXPORTGPUREG hoLinearResampleOperator : public hoResampleOperator<REAL,T,D>
{
  
 public:
  
  hoLinearResampleOperator() : hoResampleOperator<REAL,T,D>() {}
  virtual ~hoLinearResampleOperator() {}
  
  virtual int mult_M( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate = false);
  virtual int mult_MH( hoNDArray<T> *in, hoNDArray<T> *out, bool accumulate = false);

  //
  // If 'identity_dimension' is non-negative it indicates that a frame of "zero displacements" 
  // should be inserted at this dimension. A negative value (deafault) is ignored.
  //

  virtual bool set_displacement_field( boost::shared_ptr< hoNDArray<REAL> > offsets );

  virtual unsigned int get_temporal_dimension_size() { return temporal_dim_size_; }

  virtual boost::shared_ptr< linearOperator< REAL, hoNDArray<T> > > clone() {
    return linearOperator< REAL, hoNDArray<T> >::clone(this);
  }

 private:
  inline bool is_border_pixel( typename reald<REAL,D>::Type co, typename uintd<D>::Type dims );
  inline unsigned int get_num_neighbors();
  
 protected:
  boost::shared_ptr< Eigen::SparseMatrix<REAL> > R_;
  unsigned int temporal_dim_size_;
};
