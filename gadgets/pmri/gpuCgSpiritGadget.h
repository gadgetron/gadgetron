#ifndef gpuCgSpiritGadget_H
#define gpuCgSpiritGadget_H
#pragma once

#include "gadgetron_gpupmri_export.h"
#include "Gadget.h"
#include "GenericReconJob.h"
#include "GadgetMRIHeaders.h"
#include "cuCgSolver.h"
#include "cuNFFT.h"
#include "../../toolboxes/nfft/NFFTOperator.h"
#include "cuSpiritOperator.h"
#include "cuCgPreconditioner.h"
#include "cuNFFT.h"
#include "cuImageOperator.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include "gpuSenseGadget.h"

namespace Gadgetron{

  class EXPORTGADGETS_GPUPMRI gpuCgSpiritGadget : public gpuSenseGadget
   {

  public:

    GADGET_DECLARE(gpuCgSpiritGadget);

    gpuCgSpiritGadget();
    virtual ~gpuCgSpiritGadget();
  protected:
    GADGET_PROPERTY(number_of_iterations, int, "Number of iterations", 5);
    GADGET_PROPERTY(kappa, float, "Kappa regularization factor", 0.3);
    GADGET_PROPERTY(cg_limit, float, "Residual limit for CG convergence", 1e-6);
    virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader > *m1, GadgetContainerMessage< GenericReconJob > *m2 );
    virtual int process_config( ACE_Message_Block* mb );


    unsigned int number_of_iterations_;
    double cg_limit_;
   bool matrix_size_reported_;
    bool is_configured_;

    double kappa_;

    // Define conjugate gradient solver
    cuCgSolver<float_complext> cg_;

    // Define Spirit encoding operator (NFFT)
    boost::shared_ptr< NFFTOperator<cuNDArray,float,2> > E_;

    // Define Spirit regularization operator (convolution consistency)
    boost::shared_ptr< cuSpirit2DOperator<float> > S_;

    // Define preconditioner
    //boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Define regularization image operator
    //boost::shared_ptr< cuImageOperator<float_complext> > R_;
    
  };
}
#endif //gpuCgSpiritGadget
