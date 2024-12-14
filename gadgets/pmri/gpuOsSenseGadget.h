#ifndef gpuSbSenseGadget_H
#define gpuSbSenseGadget_H
#pragma once

#include "Gadget.h"
#include "GenericReconJob.h"
#include "cuNlcgSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuPartialDerivativeOperator.h"
#include "cuNFFT.h"
#include "cuImageOperator.h"
#include "cuTvOperator.h"
#include "cuTvPicsOperator.h"
#include "osSenseOperator.h"
#include "../../toolboxes/nfft/NFFTOperator.h"
#include "osMOMSolver.h"
#include "osSPSSolver.h"
#include "gpuSenseGadget.h"

#include <complex>

namespace Gadgetron{

  class gpuOsSenseGadget : public gpuSenseGadget
  {

  public:

    gpuOsSenseGadget();
    virtual ~gpuOsSenseGadget();

  protected:
    GADGET_PROPERTY(lambda, float, "Lambda regularization factor", 1e-6);
    GADGET_PROPERTY(alpha, float, "Alpha regularization factor", 0.5);
    GADGET_PROPERTY(kappa, float, "Kappa regularization factor", 1.0);
    GADGET_PROPERTY(number_of_iterations, int, "Number of solver iterations", 0);
    GADGET_PROPERTY(exclusive_access, bool,"Forces 1 gadget per GPU",false);
    GADGET_PROPERTY(coils_per_subset, int,"Number of coils to use for each subset",1);
    GADGET_PROPERTY(damping,float,"Relative step size. Reduce if solver fails to converge",1);
    
    virtual int process( GadgetContainerMessage< mrd::ImageHeader >* m1, GadgetContainerMessage< GenericReconJob > * m2 );
    virtual int process_config(const mrd::Header& header);

    int coils_per_subset_;

    uint64d2 matrix_size_;
    uint64d2 matrix_size_os_;
    uint64d2 matrix_size_seq_;

    unsigned int number_of_iterations_;


    bool exclusive_access_;
    double lambda_;
    double alpha_;
    double kappa_;
    double damping_;

    bool is_configured_;
    bool prepared_;

    // Define non-linear conjugate gradient solver
    osMOMSolver<cuNDArray<float_complext>> solver_;

    // Define non-Cartesian Sense Encoding operator
    boost::shared_ptr< osSenseOperator<cuNDArray<float_complext>,2,NFFTOperator<cuNDArray,float,2>>> E_;

    // Average image for regularization
    boost::shared_ptr< cuNDArray<float_complext> > reg_image_;

    boost::shared_ptr<cuTvOperator<float_complext,3> > TV_;
    boost::shared_ptr<cuTvPicsOperator<float_complext,3> > PICS_;

  };
}
#endif //gpuSbSenseGadget

