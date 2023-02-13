#include "gpuRadialSpiritPrepGadget.h"
#include "spirit_calibration.h"
#include "cuSpiritBuffer.h"
#include "cuNDFFT.h"
#include "cuSpiritOperator.h"
#include "hoNDArray_fileio.h"

namespace Gadgetron{

  gpuRadialSpiritPrepGadget::gpuRadialSpiritPrepGadget() : gpuRadialPrepGadget() {}

  int 
  gpuRadialSpiritPrepGadget::process_config(ACE_Message_Block* mb)
  {
    return gpuRadialPrepGadget::process_config(mb);
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuRadialSpiritPrepGadget::compute_csm( unsigned int idx )
  {    
    // Estimate and update csm related data structures
    //
  
    cuSpiritBuffer<float,2> *acc_buffer = &this->acc_buffer_spirit_[idx];
  
    boost::shared_ptr< cuNDArray<float_complext> > csm_data = 
      acc_buffer->get_accumulated_coil_images();

    std::vector<size_t> dims_to_xform;
    dims_to_xform.push_back(0); dims_to_xform.push_back(1);    
    cuNDFFT<float>::instance()->fft( csm_data.get(), &dims_to_xform );
    
    boost::shared_ptr< cuNDArray<float_complext> > csm =       
      estimate_spirit_kernels( csm_data.get(), 7 ); // TODO: let the kernel size be user defined

    return csm->to_host(); 
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuRadialSpiritPrepGadget::compute_reg( unsigned int set, unsigned int slice, bool new_frame )
  {    
    // Estimate and update regularization image related data structures
    //
    
    cuSpiritBuffer<float,2> *acc_buffer = &this->acc_buffer_spirit_[set*this->slices_+slice];
    boost::shared_ptr< cuNDArray<float_complext> > reg_image = acc_buffer->get_combined_coil_image();
    
    if( !reg_image.get() ){
      GDEBUG("Error computing regularization image\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }            
    
    return reg_image->to_host();
  }

  void 
  gpuRadialSpiritPrepGadget::allocate_accumulation_buffer( unsigned int size )
  {    
    this->acc_buffer_spirit_ = boost::shared_array< cuSpiritBuffer<float,2> >(new cuSpiritBuffer<float,2>[size]);
  }

  void gpuRadialSpiritPrepGadget::reconfigure(unsigned int set, unsigned int slice, bool use_dcw)
  {    
    gpuRadialPrepGadget::reconfigure(set, slice, use_dcw);
    //gpuRadialPrepGadget::reconfigure(set, slice, false);
    
    cuSpiritBuffer<float,2> *acc_buffer = &this->acc_buffer_spirit_[set*this->slices_+slice];

    if( use_dcw ) 
      acc_buffer->set_dcw_for_rhs(calculate_density_compensation_for_rhs(set, slice));

    acc_buffer->preprocess(calculate_trajectory_for_rhs(0, set, slice).get());
  }

  GADGET_FACTORY_DECLARE(gpuRadialSpiritPrepGadget)
}
