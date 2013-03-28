#pragma once

#include "cuNonCartesianSenseOperator.h"
#include "vector_td.h"
#include "complext.h"
#include "gpupmri_export.h"

namespace Gadgetron{

  template<class REAL, unsigned int D, bool ATOMICS = false> class EXPORTGPUPMRI cuSenseRHSBuffer
  {
  public:

    typedef complext<REAL> _complext;
    typedef typename uintd<D>::Type _uintd;
    typedef typename reald<REAL,D>::Type _reald;

    cuSenseRHSBuffer() 
      {
	num_coils_ = 0;
	cur_idx_ = cur_sub_idx_ = 0;
    
	cycle_length_ = 5; sub_cycle_length_ = 8; // This gives a data buffer of (cycle_length_-1)*sub_cycle_length_ frames
	acc_buffer_empty_ = true;
      }
  
    virtual ~cuSenseRHSBuffer() {}
  
    inline void set_num_coils( unsigned int num_coils ){ num_coils_ = num_coils; }
    inline void set_cycle_lengths( unsigned int num_cycles, unsigned int num_sub_cycles ){
      cycle_length_ = num_cycles; sub_cycle_length_ = num_sub_cycles;
    }
  
    virtual void clear();
    virtual void set_sense_operator( boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D,ATOMICS> > op );
    virtual void add_frame_data( cuNDArray<_complext> *samples, cuNDArray<_reald> *trajectory );
    virtual boost::shared_ptr< cuNDArray<_complext> > get_acc_coil_images( bool normalize = false );

  private:
    unsigned int num_coils_;
    unsigned int cycle_length_, sub_cycle_length_;
    unsigned int cur_idx_, cur_sub_idx_;
    bool acc_buffer_empty_;
    cuNDArray<_complext> acc_buffer_, cyc_buffer_;
    boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D,ATOMICS> > sense_op_;
  };

  //To prevent the use of Atomics with doubles.
  template<unsigned int D> class EXPORTGPUPMRI cuSenseRHSBuffer<double,D,true>{};
}
