#ifndef gpuCgSenseGadget_H
#define gpuCgSenseGadget_H
#pragma once

#include <ace/Synch.h>
#include <ace/Mutex.h>

#include "gadgetron_gpusense_export.h"
#include "Gadget.h"
#include "SenseJob.h"
#include "GadgetMRIHeaders.h"
#include "cuCgSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuNFFT.h"
#include "cuImageOperator.h"
#include "ismrmrd.h"

#include <complex>

namespace Gadgetron{

  class EXPORTGADGETS_GPUSENSE gpuCgSenseGadget : public Gadget2< ISMRMRD::ImageHeader, SenseJob >
  {

  public:
    GADGET_DECLARE(gpuCgSenseGadget);

    gpuCgSenseGadget();
    virtual ~gpuCgSenseGadget();

  protected:

    virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader >* m1, GadgetContainerMessage< SenseJob > * m2 );
    virtual int process_config( ACE_Message_Block* mb );

    int channels_;
    int device_number_;
    int slice_number_;

    uintd2 matrix_size_;
    uintd2 matrix_size_os_;
    uintd2 matrix_size_seq_;

    unsigned int number_of_iterations_;
    double cg_limit_;
    double oversampling_factor_;
    double kernel_width_;
    double kappa_;

    bool is_configured_;

    // Define conjugate gradient solver
    cuCgSolver<float_complext> cg_;

    // Define non-Cartesian Sense Encofing operator
    boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_;

    // Define preconditioner
    boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Define regularization image operator
    boost::shared_ptr< cuImageOperator<float_complext> > R_;

    int image_series_;
    int image_counter_;
  };
}
#endif //gpuCgSenseGadget
