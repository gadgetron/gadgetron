#pragma once

#include "gadgetron_export.h"
#include "cuNonCartesianSenseOperator.h"
#include "vector_td.h"

template<class REAL, unsigned int D>
class EXPORTGPUPMRI cuSenseRHSBuffer
{
 public:
  
  typedef typename complext<REAL>::Type _complext;
  typedef typename uintd<D>::Type _uintd;
  typedef typename reald<REAL,D>::Type _reald;

  cuSenseRHSBuffer( int device = -1 )
  {
    if( device<0 ){
      if( cudaGetDevice( &device_ ) != cudaSuccess ){
	    std::cerr << "cuSenseRHSBuffer:: unable to get current device." << std::endl;
	    device_ = 0;
      }
	}
    else
      device_ = device;
 
	cur_idx_ = cur_sub_idx_ = 0;
	cycle_length_ = 4; sub_cycle_length_ = 8;
	acc_buffer_empty_ = true;
  }

  virtual ~cuSenseRHSBuffer() {}

  virtual int set_sense_operator( boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D> > op );
  virtual int add_frame_data( cuNDArray<_complext> *samples, cuNDArray<_reald> *trajectory );
  virtual boost::shared_ptr< cuNDArray<_complext> > get_acc_coil_images();
  virtual boost::shared_ptr< cuNDArray<_complext> > get_cur_rhs();

private:
	unsigned int cycle_length_, sub_cycle_length_;
	unsigned int cur_idx_, cur_sub_idx_;
	bool acc_buffer_empty_;
	int device_;
	cuNDArray<_complext> acc_buffer_, cyc_buffer_;
	boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D> > sense_op_;
};
