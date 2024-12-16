#pragma once

#include "cuSenseBuffer.h"
#include "cuCgSolver.h"
#include "cuCgPreconditioner.h"

namespace Gadgetron{

  template<class REAL, unsigned int D>
  class cuSenseBufferCg : public cuSenseBuffer<REAL,D>
  {
  public:

    typedef complext<REAL> _complext;
    typedef typename uint64d<D>::Type _uint64d;
    typedef typename reald<REAL,D>::Type _reald;

    cuSenseBufferCg() : cuSenseBuffer<REAL,D>() {}
    virtual ~cuSenseBufferCg() {}

    inline void set_dcw_for_rhs( boost::shared_ptr< cuNDArray<REAL> > dcw ){
      this->E_->set_dcw(sqrt(dcw.get()));
    }

    virtual void preprocess( cuNDArray<_reald> *traj );

    virtual void setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W,
			unsigned int num_coils, unsigned int num_cycles, unsigned int num_sub_cycles );

    virtual boost::shared_ptr< cuNDArray<_complext> > get_combined_coil_image();

  protected:
    cuCgSolver<_complext> cg_;
    boost::shared_ptr< cuCgPreconditioner<_complext> > D_;
  };

}
