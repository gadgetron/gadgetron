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

#include <complex>

namespace Gadgetron{

  class gpuNlcgSenseGadget : public Gadget2< mrd::ImageHeader, GenericReconJob >
  {

  public:

    gpuNlcgSenseGadget();
    virtual ~gpuNlcgSenseGadget();

  protected:
    GADGET_PROPERTY(deviceno, int, "GPU device number", 0);
    GADGET_PROPERTY(setno, int, "Which set to process", 0);
    GADGET_PROPERTY(sliceno, int, "Which slice to process", 0);
    GADGET_PROPERTY(cg_limit, float, "Convervence limit for CG", 1e-6);
    GADGET_PROPERTY(oversampling_factor, float, "Oversampling factor for NFFT", 1.5);
    GADGET_PROPERTY(kernel_width, float, "Kernel width for NFFT", 5.5);
    GADGET_PROPERTY(lambda, float, "Lambda regularization factor", 1e-6);
    GADGET_PROPERTY(alpha, float, "Alpha regularization factor", 0.5);
    GADGET_PROPERTY(exclusive_access, bool, "Exclusive access to solver", false);
    GADGET_PROPERTY(number_of_cg_iterations, int, "Number of CG iterations", 0);
    GADGET_PROPERTY(rotations_to_discard, int, "Rotations to discard", 0);
    GADGET_PROPERTY(output_convergence, bool, "Output convergence information", false);
    
    virtual int process( GadgetContainerMessage< mrd::ImageHeader >* m1, GadgetContainerMessage< GenericReconJob > * m2 );
    virtual int process_config(const mrd::Header& header);

    int channels_;
    int device_number_;
    int set_number_;
    int slice_number_;

    uint64d2 matrix_size_;
    uint64d2 matrix_size_os_;
    uint64d2 matrix_size_seq_;

    unsigned int number_of_cg_iterations_;

    double cg_limit_;
    double oversampling_factor_;
    double kernel_width_;

    double lambda_;
    double alpha_;
    unsigned int rotations_to_discard_;

    bool output_convergence_;
    bool exclusive_access_;
    bool is_configured_;
    bool prepared_;

    // Define non-linear conjugate gradient solver
    cuNlcgSolver<float_complext> solver_;

    // Define non-Cartesian Sense Encoding operator
    boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_;

    // Define preconditioner
    boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Average image for regularization
    boost::shared_ptr< cuNDArray<float_complext> > reg_image_;

    boost::shared_ptr<cuTvOperator<float_complext,3> > TV_;
    boost::shared_ptr<cuTvPicsOperator<float_complext,3> > PICS_;

    int frame_counter_;
  };
}
#endif //gpuSbSenseGadget

