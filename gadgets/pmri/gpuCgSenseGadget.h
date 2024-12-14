#ifndef gpuCgSenseGadget_H
#define gpuCgSenseGadget_H
#pragma once

#include "Gadget.h"
#include "GenericReconJob.h"
#include "cuCgSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuNFFT.h"
#include "cuImageOperator.h"

#include <complex>
#include "gpuSenseGadget.h"

namespace Gadgetron{

  class gpuCgSenseGadget : public gpuSenseGadget
   {

  public:
    gpuCgSenseGadget();
    virtual ~gpuCgSenseGadget();

  protected:
    GADGET_PROPERTY(kappa, float, "Regularization factor kappa", 0.3);
    GADGET_PROPERTY(number_of_iterations, int, "Max number of iterations in CG solver", 5);
    GADGET_PROPERTY(cg_limit, float, "Residual limit for CG convergence", 1e-6);

    virtual int process( GadgetContainerMessage< mrd::ImageHeader > *m1, GadgetContainerMessage< GenericReconJob > *m2 );
    virtual int process_config(const mrd::Header& header);

    unsigned int number_of_iterations_;
    double cg_limit_;
    double kappa_;

    bool output_timing_;
    bool matrix_size_reported_;
    bool is_configured_;

    // Define conjugate gradient solver
    cuCgSolver<float_complext> cg_;

    // Define non-Cartesian Sense Encoding operator
    boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_;

    // Define preconditioner
    boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Define regularization image operator
    boost::shared_ptr< cuImageOperator<float_complext> > R_;

  };
}
#endif //gpuCgSenseGadget
