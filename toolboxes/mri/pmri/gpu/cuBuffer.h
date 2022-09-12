#pragma once

#include "vector_td_utilities.h"
#include "complext.h"
#include "cuNDArray.h"
#include "cuNFFT.h"
#include "gpupmri_export.h"

#include <boost/shared_ptr.hpp>

namespace Gadgetron{
  
  template<class REAL, unsigned int D> class EXPORTGPUPMRI cuBuffer
  {
  public:
    
    typedef complext<REAL> _complext;
    typedef typename uint64d<D>::Type _uint64d;
    typedef typename reald<REAL,D>::Type _reald;

    cuBuffer();
    virtual ~cuBuffer() {}
    
    virtual void set_dcw( boost::shared_ptr< cuNDArray<REAL> > dcw ){
      dcw_ = dcw;
    }
    
    inline REAL get_normalization_factor(){
      return REAL(1)/(((REAL)cycle_length_-REAL(1))*(REAL)sub_cycle_length_);
    }
    
    virtual void clear();

    virtual void setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W, 
                        unsigned int num_coils, unsigned int num_cycles, unsigned int num_sub_cycles );

    // Boolean return value indicates whether the accumulation buffer has changed (i.e. a cycle has been completed)
    virtual bool add_frame_data( cuNDArray<_complext> *samples, cuNDArray<_reald> *trajectory ); 

    virtual boost::shared_ptr< cuNDArray< complext<REAL> > > get_accumulated_coil_images();

    // Workaround for weird boost/g++ error
    virtual boost::shared_ptr< cuNDArray< complext<REAL> > > get_combined_coil_image() = 0;
    
  protected:
    _uint64d matrix_size_, matrix_size_os_;
    REAL W_;
    unsigned int num_coils_;
    unsigned int cycle_length_, sub_cycle_length_;
    unsigned int cur_idx_, cur_sub_idx_;
    bool acc_buffer_empty_;
    boost::shared_ptr< cuNDArray<_complext> > acc_buffer_;
    boost::shared_ptr< cuNDArray<_complext> > cyc_buffer_;
    boost::shared_ptr< cuNDArray<_complext> > acc_image_;
    boost::shared_ptr< cuNDArray<REAL> > dcw_;
    boost::shared_ptr< cuNFFT_plan<REAL,D> > nfft_plan_;
  };

  // To prevent the use of atomics with doubles.
}
