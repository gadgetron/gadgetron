#pragma once

#include "cuNonCartesianSenseOperator.h"
#include "vector_td_utilities.h"
#include "complext.h"
#include "gpupmri_export.h"

namespace Gadgetron{
  
  template<class REAL, unsigned int D, bool ATOMICS = false> class EXPORTGPUPMRI cuSenseBuffer
  {
  public:
    
    typedef complext<REAL> _complext;
    typedef typename uintd<D>::Type _uintd;
    typedef typename reald<REAL,D>::Type _reald;

    cuSenseBuffer();
    virtual ~cuSenseBuffer() {}
    
    virtual void set_csm( boost::shared_ptr< cuNDArray<_complext> > csm ){
      csm_ = csm;
    }

    virtual void set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw ){
      dcw_ = dcw;
    }
    
    virtual void clear();
    inline void setup( _uintd matrix_size, _uintd matrix_size_os, unsigned int num_coils, unsigned int num_cycles, unsigned int num_sub_cycles );
    virtual void add_frame_data( cuNDArray<_complext> *samples, cuNDArray<_reald> *trajectory );
    virtual boost::shared_ptr< cuNDArray<_complext> > get_acc_coil_images( bool normalize = false );
    
  private:
    _uintd matrix_size_, matrix_size_os_;
    unsigned int num_coils_;
    unsigned int cycle_length_, sub_cycle_length_;
    unsigned int cur_idx_, cur_sub_idx_;
    bool acc_buffer_empty_;
    cuNDArray<_complext> acc_buffer_, cyc_buffer_;
    boost::shared_ptr< cuNDArray<_complext> > csm_;
    boost::shared_ptr< cuNDArray<REAL> dcw_;
    boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D,ATOMICS> > E_;
  };
  
  // To prevent the use of atomics with doubles.
  template<unsigned int D> class EXPORTGPUPMRI cuSenseBuffer<double,D,true>{};
}
