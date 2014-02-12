#pragma once

#include "cuBuffer.h"
#include "cuCgSolver.h"
#include "cuNFFTOperator.h"

namespace Gadgetron{

  template<class REAL, unsigned int D, bool ATOMICS = false> 
  class EXPORTGPUPMRI cuSpiritBuffer : public cuBuffer<REAL,D,ATOMICS>
  {
  public:
    
    typedef typename cuBuffer<REAL,D,ATOMICS>::_complext _complext;
    typedef typename cuBuffer<REAL,D,ATOMICS>::_uint64d  _uint64d;
    typedef typename cuBuffer<REAL,D,ATOMICS>::_reald    _reald;

    cuSpiritBuffer() : cuBuffer<REAL,D,ATOMICS>() {
      E_ = boost::shared_ptr< cuNFFTOperator<REAL,D> >(new cuNFFTOperator<REAL,D>() );
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
    boost::shared_ptr< cuNFFTOperator<REAL,D> > E_;
  };
  
  // To prevent the use of atomics with doubles.
  template<unsigned int D> class EXPORTGPUPMRI cuSpiritBuffer<double,D,true>{};  
}
