#include "gpuRadialSensePrepGadget.h"
#include "b1_map.h"
#include "cuSenseBufferCg.h"

namespace Gadgetron{

  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuRadialSensePrepGadget::compute_csm( unsigned int idx )
  {    
    // Estimate and update csm related data structures
    //
  
    cuSenseBuffer<float,2> *acc_buffer = 
      (this->buffer_using_solver_) ? &this->acc_buffer_sense_cg_[idx] : &this->acc_buffer_sense_[idx];
  
    boost::shared_ptr< cuNDArray<float_complext> > csm_data = 
      acc_buffer->get_accumulated_coil_images();
    
    if( !csm_data.get() ){
      GDEBUG("Error during accumulation buffer computation\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }
    
    auto csm = boost::make_shared<cuNDArray<float_complext>>(estimate_b1_map<float,2>( *csm_data ));
  
    if( !csm.get() ){
      GDEBUG("Error during coil estimation\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }            
    
    acc_buffer->set_csm(csm);
    return csm->to_host(); 
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuRadialSensePrepGadget::compute_reg( unsigned int set, unsigned int slice, bool new_frame )
  {    
    // Estimate and update regularization image related data structures
    //
    
    cuSenseBuffer<float,2> *acc_buffer = (this->buffer_using_solver_) ? 
      &this->acc_buffer_sense_cg_[set*this->slices_+slice] : &this->acc_buffer_sense_[set*this->slices_+slice];

    if( buffer_using_solver_ && ( mode_ == 2 || mode_ == 3 ) ){
      static_cast<cuSenseBufferCg<float,2>*>( acc_buffer )->preprocess
        ( calculate_trajectory_for_rhs( this->profiles_counter_global_[set*this->slices_+slice] - ((new_frame) ? 1 : 0), set, slice).get());
    }
    
    boost::shared_ptr< cuNDArray<float_complext> > reg_image = 
      acc_buffer->get_combined_coil_image();
    
    if( !reg_image.get() ){
      GDEBUG("Error computing regularization image\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }            
    
    return reg_image->to_host();
  }

  void 
  gpuRadialSensePrepGadget::allocate_accumulation_buffer( unsigned int size )
  {    
    // Allocate accumulation buffer
    //
  
    if( this->buffer_using_solver_ ){
      this->acc_buffer_sense_cg_ = boost::shared_array< cuSenseBufferCg<float,2> >(new cuSenseBufferCg<float,2>[size]);
    }
    else{
      this->acc_buffer_sense_ = boost::shared_array< cuSenseBuffer<float,2> >(new cuSenseBuffer<float,2>[size]);
    }
  }

  void gpuRadialSensePrepGadget::reconfigure(unsigned int set, unsigned int slice, bool use_dcw)
  {    
    gpuRadialPrepGadget::reconfigure(set, slice, use_dcw);
    
    if( buffer_using_solver_ ){

      if(use_dcw) 
        this->acc_buffer_sense_cg_[set*this->slices_+slice].set_dcw_for_rhs(calculate_density_compensation_for_rhs(set, slice));

      auto t = calculate_trajectory_for_rhs(0, set, slice);
      auto tmp = this->acc_buffer_sense_cg_[set*this->slices_+slice];
      tmp.preprocess(t.get());
      this->acc_buffer_sense_cg_[set*this->slices_+slice].preprocess(t.get()); //calculate_trajectory_for_rhs(0, set, slice).get());
    }    
  }

  GADGET_FACTORY_DECLARE(gpuRadialSensePrepGadget)
}
