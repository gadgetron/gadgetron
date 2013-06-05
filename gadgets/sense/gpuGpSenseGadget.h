#ifndef gpuGpSenseGadget_H
#define gpuGpSenseGadget_H
#pragma once

#include "Gadget.h"
#include "SenseJob.h"
#include "GadgetMRIHeaders.h"
#include "cuGpBbSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuTvOperator.h"
#include "cuTvPicsOperator.h"
#include "cuNFFT.h"
#include "cuSenseRHSBuffer.h"
#include "cuImageOperator.h"
#include "gadgetron_gpusense_export.h"

#include <ismrmrd.h>
#include <complex>

#include <ace/Synch.h>
#include <ace/Mutex.h>

namespace Gadgetron{

  class EXPORTGADGETS_GPUSENSE gpuGpSenseGadget : public Gadget2< ISMRMRD::ImageHeader, SenseJob >
  {

  public:
    GADGET_DECLARE(gpuGpSenseGadget);

    gpuGpSenseGadget();
    virtual ~gpuGpSenseGadget();

  protected:

    virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader >* m1, GadgetContainerMessage< SenseJob > * m2 );
    virtual int process_config( ACE_Message_Block* mb );

    int channels_;
    int device_number_;

    uintd2 matrix_size_;
    uintd2 matrix_size_os_;

    unsigned int number_of_iterations_;
    double oversampling_;
    double kernel_width_;
    double lambda_;
    double alpha_;

    bool is_configured_;
    bool prepared_;

    // Define constraint Split Bregman solver
    cuGpBbSolver<float_complext> gp_;

    // Define non-Cartesian Sense Encofing operator
    boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_;

    // Define preconditioner
    boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Average image for regularization
    boost::shared_ptr< cuNDArray<float_complext> > reg_image_;

    // Define regularization operators
    boost::shared_ptr<cuTvOperator<float_complext,3> > TV_;
    boost::shared_ptr<cuTvPicsOperator<float_complext,3> > TVP_;
	
    int image_series_;
    int image_counter_;
  };
}
#endif //gpuGpSenseGadget
