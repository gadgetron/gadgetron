#ifndef GPUCGGadgetGeneric_H
#define GPUCGGadgetGeneric_H
#pragma once

#include <ace/Synch.h>
#include <ace/Mutex.h>

#include "gadgetroncgsense_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "cuCgSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "NFFT.h"
#include "cuSenseRHSBuffer.h"
#include "cuImageOperator.h"
#include "ismrmrd.h"

#include <complex>

namespace Gadgetron{

class CGSenseJob
{
public:

	CGSenseJob() {}

	~CGSenseJob() {

	}

	boost::shared_ptr< hoNDArray<float_complext> > csm_host_;
	boost::shared_ptr< hoNDArray<float_complext> > reg_host_;
	boost::shared_ptr< hoNDArray<floatd2>        > tra_host_;
	boost::shared_ptr< hoNDArray<float>          > dcw_host_;
	boost::shared_ptr< hoNDArray<float_complext> > dat_host_;
};


class EXPORTGADGETSCGSENSE GPUCGGadgetGeneric : public Gadget2< ISMRMRD::ImageHeader, CGSenseJob >
{

public:
    GADGET_DECLARE(GPUCGGadgetGeneric);

	GPUCGGadgetGeneric();
	virtual ~GPUCGGadgetGeneric();

protected:

	virtual int process( GadgetContainerMessage< ISMRMRD::ImageHeader >* m1, GadgetContainerMessage< CGSenseJob > * m2 );
	virtual int process_config( ACE_Message_Block* mb );

	virtual int configure_channels();

	int channels_;
	int device_number_;

	uintd2 matrix_size_;
	uintd2 matrix_size_os_;

	unsigned int number_of_iterations_;
	double cg_limit_;
	double oversampling_;
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
#endif //GPUCGGadgetGeneric
