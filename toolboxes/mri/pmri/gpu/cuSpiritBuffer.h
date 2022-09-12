#pragma once

#include "cuBuffer.h"
#include "cuCgSolver.h"
#include "../../../nfft/NFFTOperator.h"

namespace Gadgetron{

      template<class REAL, unsigned int D = false>
  class EXPORTGPUPMRI cuSpiritBuffer : public cuBuffer<REAL,D>
  {
  public:
    
    typedef typename cuBuffer<REAL,D>::_complext _complext;
    typedef typename cuBuffer<REAL,D>::_uint64d  _uint64d;
    typedef typename cuBuffer<REAL,D>::_reald    _reald;

    cuSpiritBuffer() : cuBuffer<REAL,D>() {
      E_ = boost::make_shared< NFFTOperator<cuNDArray,REAL,D> >();
    }
    
    virtual ~cuSpiritBuffer() {}
    
    inline void set_dcw_for_rhs( boost::shared_ptr< cuNDArray<REAL> > dcw ){
      this->E_->set_dcw(dcw);
    }

    virtual void setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W, 
                        unsigned int num_coils, unsigned int num_cycles, unsigned int num_sub_cycles );
    
    virtual void preprocess( cuNDArray<_reald> *traj );

    virtual boost::shared_ptr< cuNDArray< complext<REAL> > > get_accumulated_coil_images();
    virtual boost::shared_ptr< cuNDArray< complext<REAL> > > get_combined_coil_image();
    
  protected:
    cuCgSolver<_complext> cg_;
    boost::shared_ptr< NFFTOperator<cuNDArray,REAL,D> > E_;
  };
  
}
