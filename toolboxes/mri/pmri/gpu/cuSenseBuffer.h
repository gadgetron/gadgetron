#pragma once

#include "cuBuffer.h"
#include "cuNonCartesianSenseOperator.h"

#include <stdio.h>

namespace Gadgetron{

  template<class REAL, unsigned int D>
  class EXPORTGPUPMRI cuSenseBuffer : public cuBuffer<REAL,D>
  {
  public:
    
    typedef typename cuBuffer<REAL,D>::_complext _complext;
    typedef typename cuBuffer<REAL,D>::_uint64d  _uint64d;
    typedef typename cuBuffer<REAL,D>::_reald    _reald;

    cuSenseBuffer() : cuBuffer<REAL,D>() {}
    virtual ~cuSenseBuffer() {}

    virtual void setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W, 
                        unsigned int num_coils, unsigned int num_cycles, unsigned int num_sub_cycles );

    virtual void set_csm( boost::shared_ptr< cuNDArray<_complext> > csm ){
      csm_ = csm;
    }
    
    virtual boost::shared_ptr< cuNDArray< complext<REAL> > > get_combined_coil_image();

  protected:
    boost::shared_ptr< cuNDArray<_complext> > csm_;
    boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D> > E_;
  };
  
}
