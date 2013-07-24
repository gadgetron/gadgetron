#include "cpuRegistrationAveragingGadget.h"
#include "hoLinearResampleOperator.h"
#include "hoCKOpticalFlowSolver.h"

namespace Gadgetron{

  int cpuRegistrationAveragingGadget2D::setup_solver()
  {
    // Allocate solver
    hoCKOpticalFlowSolver<float,2> *solver = new hoCKOpticalFlowSolver<float,2>();
    this->of_solver_ = solver;

    // Use bilinear resampling for interpolation
    solver->set_interpolator( boost::shared_ptr< hoLinearResampleOperator<float,2> >(new hoLinearResampleOperator<float,2>()) );
    
    // Configurable settings from the xml propoerties
    //
    
    if( this->output_convergence_ )
      solver->set_output_mode( hoCKOpticalFlowSolver<float,2>::OUTPUT_VERBOSE );
    else
      solver->set_output_mode( hoCKOpticalFlowSolver<float,2>::OUTPUT_SILENT );
    
    solver->set_num_multires_levels(this->num_multires_levels_);
    solver->set_max_num_iterations_per_level(this->max_iterations_per_level_);
    solver->set_alpha(this->alpha_);
    solver->set_beta(this->beta_);
    solver->set_limit(this->limit_);

    return GADGET_OK;
  }

  int cpuRegistrationAveragingGadget2D::set_continuation
  ( GadgetContainerMessage<ISMRMRD::ImageHeader>* m1, hoNDArray<float> *continuation )
  {
    GadgetContainerMessage< hoNDArray<float> > *m2 = new GadgetContainerMessage< hoNDArray<float> >();      
    *m2->getObjectPtr() = *continuation;
    m1->cont(m2);

    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(cpuRegistrationAveragingGadget2D)
}
