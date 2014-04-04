#ifndef gpuSbSenseGadget_H
#define gpuSbSenseGadget_H
#pragma once

#include <ace/Synch.h>
#include <ace/Mutex.h>

#include "gadgetron_gpupmri_export.h"
#include "Gadget.h"
#include "GenericReconJob.h"
#include "GadgetMRIHeaders.h"
#include "cuNlcgSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuPartialDerivativeOperator.h"
#include "cuNFFT.h"
#include "cuImageOperator.h"
#include "ismrmrd.h"
#include "cuTvOperator.h"
#include "cuTvPicsOperator.h"

#include <complex>

namespace Gadgetron{

  class EXPORTGADGETS_GPUPMRI gpuNlcgSenseGadget : public Gadget2< ISMRMRD::ImageHeader, GenericReconJob >
  {

  public:
    GADGET_DECLARE(gpuNlcgSenseGadget);

    gpuNlcgSenseGadget();
    virtual ~gpuNlcgSenseGadget();

  protected:

    virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader >* m1, GadgetContainerMessage< GenericReconJob > * m2 );
    virtual int process_config( ACE_Message_Block* mb );

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

