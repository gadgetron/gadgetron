#include "gpuRadialSpiritPrepGadget.h"
#include "cuSpiritBuffer.h"
#include "cuSpiritCalibrationOperator.h"
#include "cuCgSolver.h"
#include "cuIdentityOperator.h"
#include "hoNDArray_fileio.h"

namespace Gadgetron{

  gpuRadialSpiritPrepGadget::gpuRadialSpiritPrepGadget() : gpuRadialPrepGadget(){
    set_parameter(std::string("number_of_iterations").c_str(), "25");
    set_parameter(std::string("cg_limit").c_str(), "1e-6");
    set_parameter(std::string("tikhonov_weight").c_str(), "0.1");
  }

  int 
  gpuRadialSpiritPrepGadget::process_config(ACE_Message_Block* mb)
  {
    gpuRadialPrepGadget::process_config(mb);
    number_of_iterations_ = get_int_value(std::string("number_of_iterations").c_str());
    cg_limit_ = get_double_value(std::string("cg_limit").c_str());
    reg_weight_ = get_double_value(std::string("tikhonov_weight").c_str());
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuRadialSpiritPrepGadget::compute_csm( unsigned int idx )
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

    std::vector<size_t> domain_dims = *csm_data->get_dimensions();
    size_t num_coils = domain_dims.back();
    domain_dims.pop_back();
    domain_dims.push_back(num_coils*num_coils);

    boost::shared_ptr< cuSpirit2DCalibrationOperator<float> > C( new cuSpirit2DCalibrationOperator<float>() );
    C->set_accumulated_kspace(csm_data);
    C->set_domain_dimensions(&domain_dims);
    C->set_codomain_dimensions(csm_data->get_dimensions().get());

    boost::shared_ptr< cuIdentityOperator<float_complext> > R( new cuIdentityOperator<float_complext>() );
    R->set_weight(reg_weight_);

    cuCgSolver<float_complext> cg;
    cg.set_encoding_operator(C);
    cg.add_regularization_operator(R);
    cg.set_max_iterations(number_of_iterations_);
    cg.set_tc_tolerance( cg_limit_ );
    cg.set_output_mode(Gadgetron::cuCgSolver<float_complext>::OUTPUT_VERBOSE);

    boost::shared_ptr< cuNDArray<float_complext> > csm = 
      cg.solve(csm_data.get());
      
    if( !csm.get() ){
      GADGET_DEBUG1("Error during spirit calibration estimation\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }            

    /*static int counter = 0;
    char filename[256];
    sprintf((char*)filename, "_before_%d.real", counter);
    write_nd_array<float>( abs(csm_data.get())->to_host().get(), filename );
    cuNDArray<float_complext> after(csm_data.get()); C->mult_M(csm.get(),&after);
    sprintf((char*)filename, "_after_%d.real", counter);
    write_nd_array<float>( abs(&after)->to_host().get(), filename );
    sprintf((char*)filename, "_spirit_calibration_%d.real", counter);
    write_nd_array<float>( abs(csm.get())->to_host().get(), filename );    
    counter++;*/
      
    return csm->to_host(); 
  }
  
  boost::shared_ptr< hoNDArray<float_complext> > 
  gpuRadialSpiritPrepGadget::compute_reg( unsigned int set, unsigned int slice, bool new_frame )
  {    
    // Estimate and update regularization image related data structures
    //
    
    cuBuffer<float,2> *acc_buffer = &this->acc_buffer_[set*this->slices_+slice];
    boost::shared_ptr< cuNDArray<float_complext> > reg_image = acc_buffer->get_combined_coil_image();
    
    if( !reg_image.get() ){
      GADGET_DEBUG1("Error computing regularization image\n");
      return boost::shared_ptr< hoNDArray<float_complext> >();
    }            
    
    return reg_image->to_host();
  }

  void 
  gpuRadialSpiritPrepGadget::allocate_accumulation_buffer( unsigned int size )
  {    
    this->acc_buffer_ = boost::shared_array< cuBuffer<float,2> >(new cuSpiritBuffer<float,2>[size]);
  }

  void gpuRadialSpiritPrepGadget::reconfigure(unsigned int set, unsigned int slice, bool use_dcw)
  {    
    gpuRadialPrepGadget::reconfigure(set, slice, use_dcw);
    //gpuRadialPrepGadget::reconfigure(set, slice, false);
    
    cuBuffer<float,2> *acc_buffer = &this->acc_buffer_[set*this->slices_+slice];

    if( use_dcw ) static_cast<cuSpiritBuffer<float,2>*>( acc_buffer )->set_dcw_for_rhs(calculate_density_compensation_for_rhs(set, slice));
    static_cast<cuSpiritBuffer<float,2>*>( acc_buffer )->preprocess(calculate_trajectory_for_rhs(0, set, slice).get());
  }
}
