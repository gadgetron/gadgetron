#include "gpuRegistrationAveragingGadget.h"
#include "cuLinearResampleOperator.h"
#include "cuCKOpticalFlowSolver.h"

namespace Gadgetron{

  int gpuRegistrationAveragingGadget2D::setup_solver()
  {
    // Allocate solver
    cuCKOpticalFlowSolver<float,2> *solver = new cuCKOpticalFlowSolver<float,2>();
    this->of_solver_ = solver;

    // Use bilinear resampling for interpolation
    solver->set_interpolator( boost::shared_ptr< cuLinearResampleOperator<float,2> >(new cuLinearResampleOperator<float,2>()) );
    
    // Configurable settings from the xml propoerties
    //
    
    if( this->output_convergence_ )
      solver->set_output_mode( cuCKOpticalFlowSolver<float,2>::OUTPUT_VERBOSE );
    else
      solver->set_output_mode( cuCKOpticalFlowSolver<float,2>::OUTPUT_SILENT );
    
    solver->set_num_multires_levels(this->num_multires_levels_);
    solver->set_max_num_iterations_per_level(this->max_iterations_per_level_);
    solver->set_alpha(this->alpha_);
    solver->set_beta(this->beta_);
    solver->set_limit(this->limit_);

    return GADGET_OK;
  }

  int gpuRegistrationAveragingGadget2D::set_continuation
  ( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, cuNDArray<float> *continuation )
  {
    GadgetContainerMessage< hoNDArray<float> > *m2 = new GadgetContainerMessage< hoNDArray<float> >();      
    m2->getObjectPtr()->create(continuation->dimensions());
    
    if( cudaMemcpy( m2->getObjectPtr()->get_data_ptr(), continuation->get_data_ptr(), 
		    continuation->get_number_of_elements()*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
      throw cuda_error("gpuRegistrationAveragingGadget::set_continuation(): failed to copy memory from device");
    }

    m1->cont(m2);

    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(gpuRegistrationAveragingGadget2D)
}
