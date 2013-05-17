#ifndef cuSbGadgetGeneric_H
#define cuSbGadgetGeneric_H
#pragma once

#include <ace/Synch.h>
#include <ace/Mutex.h>

#include "gadgetronsbsense_export.h"
#include "Gadget.h"
#include "SenseJob.h"
#include "GadgetMRIHeaders.h"
#include "cuSbcCgSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuPartialDerivativeOperator.h"
#include "cuNFFT.h"
#include "cuSenseRHSBuffer.h"
#include "cuImageOperator.h"
#include "ismrmrd.h"

#include <complex>

namespace Gadgetron{

  class EXPORTGADGETSSBSENSE cuSbGadgetGeneric : public Gadget2< ISMRMRD::ImageHeader, SenseJob >
  {

  public:
    GADGET_DECLARE(cuSbGadgetGeneric);

    cuSbGadgetGeneric();
    virtual ~cuSbGadgetGeneric();

  protected:

    virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader >* m1, GadgetContainerMessage< SenseJob > * m2 );
    virtual int process_config( ACE_Message_Block* mb );

    virtual int configure_channels();

    int channels_;
    int device_number_;

    uintd2 matrix_size_;
    uintd2 matrix_size_os_;

    unsigned int number_of_cg_iterations_;
    unsigned int number_of_sb_iterations_;
    double cg_limit_;
    double oversampling_;
    double kernel_width_;
    double mu_;
    double lambda_;
    double alpha_;

    bool is_configured_;

    // Define constraint Split Bregman solver
    cuSbcCgSolver<float_complext> sb_;

    // Define non-Cartesian Sense Encofing operator
    boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_;

    // Define preconditioner
    boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Define regularization operators
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,2> > Rx1_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,2> > Rx2_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,2> > Ry1_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,2> > Ry2_;
    //boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Rz1_;
    //boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Rz2_;
	
    int image_series_;
    int image_counter_;
  };
}
#endif //cuSbGadgetGeneric
