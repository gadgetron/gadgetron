#ifndef gpuSbSenseGadget_H
#define gpuSbSenseGadget_H
#pragma once

#include <ace/Synch.h>
#include <ace/Mutex.h>

#include "gadgetron_gpupmri_export.h"
#include "Gadget.h"
#include "GenericReconJob.h"
#include "GadgetMRIHeaders.h"
#include "cuSbcCgSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuPartialDerivativeOperator.h"
#include "cuPartialDerivativeOperator2.h"
#include "cuNFFT.h"
#include "cuImageOperator.h"
#include "ismrmrd/ismrmrd.h"
#include "gpuSenseGadget.h"
#include <complex>
#include "cuDWTOperator.h"

namespace Gadgetron{

  class EXPORTGADGETS_GPUPMRI gpuSbSenseGadget : public gpuSenseGadget
  {

  public:
    GADGET_DECLARE(gpuSbSenseGadget);

    gpuSbSenseGadget();
    virtual ~gpuSbSenseGadget();

  protected:

    virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader >* m1, GadgetContainerMessage< GenericReconJob > * m2 );
    virtual int process_config( ACE_Message_Block* mb );


    unsigned int number_of_cg_iterations_;
    unsigned int number_of_sb_iterations_;
    double cg_limit_;
    double mu_;
    double lambda_;
    double alpha_;
    double gamma_;
    unsigned int rotations_to_discard_;

    bool exclusive_access_;
    bool is_configured_;
    bool prepared_;
    bool is_cyclic_; //True if 3rd dimension of the data is cyclic (i.e. cardiac)

    // Define constraint Split Bregman solver
    cuSbcCgSolver<float_complext> sb_;

    // Define non-Cartesian Sense Encoding operator
    boost::shared_ptr< cuNonCartesianSenseOperator<float,2> > E_;

    // Define preconditioner
    boost::shared_ptr< cuCgPreconditioner<float_complext> > D_;

    // Average image for regularization
    boost::shared_ptr< cuNDArray<float_complext> > reg_image_;

    // Define regularization operators
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Rx1_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Rx2_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Ry1_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Ry2_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Rz1_;
    boost::shared_ptr< cuPartialDerivativeOperator<float_complext,3> > Rz2_;
    boost::shared_ptr< cuPartialDerivativeOperator2<float_complext,3> > Rt1_;
    boost::shared_ptr< cuPartialDerivativeOperator2<float_complext,3> > Rt2_;
    boost::shared_ptr< cuDWTOperator<float_complext,3> > W_;
    boost::shared_ptr< cuDWTOperator<float_complext,3> > W2_;
	
  };
}
#endif //gpuSbSenseGadget
