#ifndef gpuCgKtSenseGadget_H
#define gpuCgKtSenseGadget_H
#pragma once

#include "Gadget.h"
#include "GenericReconJob.h"
#include "cuCgSolver.h"
#include "cuNonCartesianKtSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuNFFT.h"
#include "cuImageOperator.h"

#include <complex>

namespace Gadgetron{

  class gpuCgKtSenseGadget : public Gadget2<mrd::ImageHeader, GenericReconJob>
  {

  public:
    gpuCgKtSenseGadget();
    virtual ~gpuCgKtSenseGadget();

  protected:
    GADGET_PROPERTY(deviceno, int, "GPU device number", 0);
    GADGET_PROPERTY(setno, int, "Set to process", 0);
    GADGET_PROPERTY(sliceno, int, "Slice to process", 0);
    GADGET_PROPERTY(number_of_iterations, int, "Number of iterations", 5);
    GADGET_PROPERTY(cg_limit, float, "Convergence limit for CG solver", 1e-6);
    GADGET_PROPERTY(oversampling_factor, float, "Recon oversampling factor for NFFT", 1.25);
    GADGET_PROPERTY(kernel_width, float, "Kernel width for NFFT", 5.5);
    GADGET_PROPERTY(kappa, float, "Kappa regularization factor", 0.3);
    GADGET_PROPERTY(training_data_shutter_radius, float, "Shutter radius for training data", 0.0);
    GADGET_PROPERTY(rotations_to_discard, int, "Number of rotations to dump", 0);
    GADGET_PROPERTY(output_convergence, bool, "Print convergence information", false);

    virtual int process( GadgetContainerMessage< mrd::ImageHeader > *m1, GadgetContainerMessage< GenericReconJob > *m2 );
    virtual int process_config(const mrd::Header& header);

    boost::shared_ptr< cuNDArray<float_complext> > compute_regularization_image( GenericReconJob *job );

    int channels_;
    int device_number_;
    int set_number_;
    int slice_number_;

    uint64d2 matrix_size_;
    uint64d2 matrix_size_os_;
    uint64d2 matrix_size_seq_;

    unsigned int number_of_iterations_;
    double cg_limit_;
    double oversampling_factor_;
    double kernel_width_;
    double kappa_;
    double shutter_radius_;
    unsigned int rotations_to_discard_;

    bool output_convergence_;
    bool is_configured_;

    // Define conjugate gradient solver
    cuCgSolver<float_complext> cg_;

    // Define non-Cartesian Sense Encoding operator
    boost::shared_ptr< cuNonCartesianKtSenseOperator<float,2> > E_;

    // Define preconditioner
    boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Define regularization image operator
    boost::shared_ptr< cuImageOperator<float_complext> > R_;

    int frame_counter_;
  };
}
#endif //gpuCgKtSenseGadget
