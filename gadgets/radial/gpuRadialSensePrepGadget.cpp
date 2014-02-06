#include "gpuRadialSensePrepGadget.h"
#include "b1_map.h"
#include "cuSenseBufferCg.h"

namespace Gadgetron{

  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuRadialSensePrepGadget::compute_csm( unsigned int idx )
  {    
    // Estimate and update csm related data structures
    //
  
    cuBuffer<float,2> *acc_buffer = &this->acc_buffer_[idx];
  
    boost::shared_ptr< cuNDArray<float_complext> > csm_data = 
      acc_buffer->get_accumulated_coil_images();
    
    if( !csm_data.get() ){
      GADGET_DEBUG1("Error during accumulation buffer computation\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }
    
    boost::shared_ptr< cuNDArray<float_complext> > csm = 
      estimate_b1_map<float,2>( csm_data.get() );
  
    if( !csm.get() ){
      GADGET_DEBUG1("Error during coil estimation\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }            
    
    static_cast<cuSenseBuffer<float,2>*>( acc_buffer )->set_csm(csm);
    return csm->to_host(); 
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuRadialSensePrepGadget::compute_reg( unsigned int set, unsigned int slice, bool new_frame )
  {    
    // Estimate and update regularization image related data structures
    //
    
    cuBuffer<float,2> *acc_buffer = &this->acc_buffer_[set*this->slices_+slice];

    if( buffer_using_solver_ && ( mode_ == 2 || mode_ == 3 ) ){
      static_cast<cuSenseBufferCg<float,2>*>( acc_buffer )->preprocess
        ( calculate_trajectory_for_rhs( this->profiles_counter_global_[set*this->slices_+slice] - ((new_frame) ? 1 : 0), set, slice).get());
    }
    
    boost::shared_ptr< cuNDArray<float_complext> > reg_image = 
      acc_buffer->get_combined_coil_image();
    
    if( !reg_image.get() ){
      GADGET_DEBUG1("Error computing regularization image\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }            
    
    return reg_image->to_host();
  }

  void 
  gpuRadialSensePrepGadget::allocate_accumulation_buffer( unsigned int size )
  {    
    // Allocate accumulation buffer
    //
  
    if( this->buffer_using_solver_ )
      this->acc_buffer_ = boost::shared_array< cuBuffer<float,2> >(new cuSenseBufferCg<float,2>[size]);
    else
      this->acc_buffer_ = boost::shared_array< cuBuffer<float,2> >(new cuSenseBuffer<float,2>[size]);
  }

  void gpuRadialSensePrepGadget::reconfigure(unsigned int set, unsigned int slice)
  {    
    gpuRadialPrepGadget::reconfigure(set, slice);
    
    cuBuffer<float,2> *acc_buffer = &this->acc_buffer_[set*this->slices_+slice];
    
    if( buffer_using_solver_ ){
      ((cuSenseBufferCg<float,2>*) acc_buffer)->set_dcw_for_rhs(calculate_density_compensation_for_rhs(set, slice));
      ((cuSenseBufferCg<float,2>*) acc_buffer)->preprocess(calculate_trajectory_for_rhs(0, set, slice).get());
    }    
  }
}
