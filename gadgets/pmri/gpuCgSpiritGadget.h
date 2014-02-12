#ifndef gpuCgSpiritGadget_H
#define gpuCgSpiritGadget_H
#pragma once

#include "gadgetron_gpupmri_export.h"
#include "Gadget.h"
#include "GenericReconJob.h"
#include "GadgetMRIHeaders.h"
#include "cuCgSolver.h"
#include "cuNFFTOperator.h"
#include "cuSpiritOperator.h"
#include "cuCgPreconditioner.h"
#include "cuNFFT.h"
#include "cuImageOperator.h"

#include <ismrmrd.h>
#include <complex>

namespace Gadgetron{

  class EXPORTGADGETS_GPUPMRI gpuCgSpiritGadget : public Gadget2<ISMRMRD::ImageHeader, GenericReconJob>
  {

  public:

    GADGET_DECLARE(gpuCgSpiritGadget);

    gpuCgSpiritGadget();
    virtual ~gpuCgSpiritGadget();

  protected:

    virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader > *m1, GadgetContainerMessage< GenericReconJob > *m2 );
    virtual int process_config( ACE_Message_Block* mb );

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
    unsigned int rotations_to_discard_;

    bool output_convergence_;
    bool output_timing_;
    bool matrix_size_reported_;
    bool is_configured_;

    // Define conjugate gradient solver
    cuCgSolver<float_complext> cg_;

    // Define Spirit encoding operator (NFFT)
    boost::shared_ptr< cuNFFTOperator<float,2> > E_;

    // Define Spirit regularization operator (convolution consistency)
    boost::shared_ptr< cuSpirit2DOperator<float> > S_;

    // Define preconditioner
    //boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Define regularization image operator
    //boost::shared_ptr< cuImageOperator<float_complext> > R_;
    
    unsigned int frame_counter_;
  };
}
#endif //gpuCgSpiritGadget
