/*
 * CSIGadget.cpp
 *
 *  Created on: Nov 11, 2014
 *      Author: dch
 */

#include "CSIGadget.h"
#include <ismrmrd/xml.h>
#include "cudaDeviceManager.h"

namespace Gadgetron {

CSIGadget::CSIGadget() {
	// TODO Auto-generated constructor stub

}

CSIGadget::~CSIGadget() {
	// TODO Auto-generated destructor stub
}


int CSIGadget::process_config(ACE_Message_Block *mb){
	//GADGET_DEBUG1("gpuCgSenseGadget::process_config\n");

	device_number_ = get_int_value(std::string("deviceno").c_str());

	int number_of_devices = 0;
	if (cudaGetDeviceCount(&number_of_devices)!= cudaSuccess) {
		GADGET_DEBUG1( "Error: unable to query number of CUDA devices.\n" );
		return GADGET_FAIL;
	}

	if (number_of_devices == 0) {
		GADGET_DEBUG1( "Error: No available CUDA devices.\n" );
		return GADGET_FAIL;
	}

	if (device_number_ >= number_of_devices) {
		GADGET_DEBUG2("Adjusting device number from %d to %d\n", device_number_,  (device_number_%number_of_devices));
		device_number_ = (device_number_%number_of_devices);
	}

	if (cudaSetDevice(device_number_)!= cudaSuccess) {
		GADGET_DEBUG1( "Error: unable to set CUDA device.\n" );
		return GADGET_FAIL;
	}

	pass_on_undesired_data_ = get_bool_value(std::string("pass_on_undesired_data").c_str());
	cg_limit_ = get_double_value(std::string("cg_limit").c_str());
	oversampling_factor_ = get_double_value(std::string("oversampling_factor").c_str());
	kernel_width_ = get_double_value(std::string("kernel_width").c_str());
	kappa_ = get_double_value(std::string("kappa").c_str());
	output_convergence_ = get_bool_value(std::string("output_convergence").c_str());
	number_of_sb_iterations_ = get_int_value(std::string("number_of_sb_iterations").c_str());
	number_of_cg_iterations_ = get_int_value(std::string("number_of_cg_iterations").c_str());

	mu_ = get_double_value(std::string("mu").c_str());

	// Get the Ismrmrd header
	//
	ISMRMRD::IsmrmrdHeader h;
	ISMRMRD::deserialize(mb->rd_ptr(),h);


	if (h.encoding.size() != 1) {
		GADGET_DEBUG1("This Gadget only supports one encoding space\n");
		return GADGET_FAIL;
	}

	// Get the encoding space and trajectory description
	ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
	ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
	ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

	img_dims_ = {r_space.matrixSize.x,r_space.matrixSize.y,r_space.matrixSize.z};

	matrix_size_ = {r_space.matrixSize.x,r_space.matrixSize.y};

	unsigned int warp_size = cudaDeviceManager::Instance()->warp_size(device_number_);

	matrix_size_os_ =
			uint64d2(((static_cast<unsigned int>(std::ceil(matrix_size_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
					((static_cast<unsigned int>(std::ceil(matrix_size_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);

	if (!is_configured_) {

		if (h.acquisitionSystemInformation) {
			channels_ = h.acquisitionSystemInformation->receiverChannels ? *h.acquisitionSystemInformation->receiverChannels : 1;
		} else {
			channels_ = 1;
		}

		// Allocate encoding operator for non-Cartesian Sense
		E_ = boost::make_shared< CSIOperator<float> >();
		E_->set_weight(mu_);
		S_ = boost::make_shared<cuNonCartesianSenseOperator<float,2>>();
		E_->set_senseOp(S_);


		auto idOp = boost::make_shared<identityOperator<cuNDArray<float_complext>>>();
		idOp->set_domain_dimensions(&img_dims_);
		idOp->set_codomain_dimensions(&img_dims_);
		idOp->set_weight(2*mu_);
		solver_.add_regularization_operator(idOp);
		// Setup solver
		solver_.set_encoding_operator( E_ );        // encoding matrix
		solver_.set_max_outer_iterations( number_of_sb_iterations_ );
		solver_.set_max_inner_iterations(1);
		solver_.get_inner_solver()->set_max_iterations(number_of_cg_iterations_);
		solver_.get_inner_solver()->set_tc_tolerance( cg_limit_ );
		solver_.set_output_mode( (output_convergence_) ? cuCgSolver<float_complext>::OUTPUT_VERBOSE : cuCgSolver<float_complext>::OUTPUT_SILENT);

		is_configured_ = true;
	}

	return GADGET_OK;

}

int CSIGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> * m1, GadgetContainerMessage<GenericReconJob> *m2){


	if (!is_configured_) {
		GADGET_DEBUG1("\nData received before configuration complete\n");
		return GADGET_FAIL;
	}

	GenericReconJob* j = m2->getObjectPtr();

	// Let's first check that this job has the required data...
	if (!j->csm_host_.get() || !j->dat_host_.get() || !j->tra_host_.get() || !j->dcw_host_.get()) {
		GADGET_DEBUG1("Received an incomplete Sense job\n");
		return GADGET_FAIL;
	}


	boost::shared_ptr< cuNDArray<floatd2> > traj(new cuNDArray<floatd2> (j->tra_host_.get()));
	boost::shared_ptr< cuNDArray<float> > dcw(new cuNDArray<float> (j->dcw_host_.get()));
	sqrt_inplace(dcw.get());
	boost::shared_ptr< cuNDArray<float_complext> > csm(new cuNDArray<float_complext> (j->csm_host_.get()));
	boost::shared_ptr< cuNDArray<float_complext> > device_samples(new cuNDArray<float_complext> (j->dat_host_.get()));

	E_->set_domain_dimensions(&img_dims_);
	E_->set_codomain_dimensions(device_samples->get_dimensions().get());

	S_->set_domain_dimensions(&img_dims_);
	S_->set_codomain_dimensions(device_samples->get_dimensions().get());

	S_->set_csm(csm);
	S_->setup( matrix_size_, matrix_size_os_, kernel_width_ );
	S_->preprocess(traj.get());

	auto result = solver_.solve(device_samples.get());


	GadgetContainerMessage< hoNDArray< std::complex<float> > > *cm =
			new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

	GadgetContainerMessage<ISMRMRD::ImageHeader> *m =
			new GadgetContainerMessage<ISMRMRD::ImageHeader>();

	*m->getObjectPtr() = j->image_headers_[0];
	m->getObjectPtr()->matrix_size[0] = matrix_size_[0];
	m->getObjectPtr()->matrix_size[1] = matrix_size_[1];
	m->getObjectPtr()->matrix_size[2] = matrix_size_[2];

	m->cont(cm);



	result->to_host((hoNDArray<float_complext>*)cm->getObjectPtr());

	m->getObjectPtr()->matrix_size[0] = img_dims_[0];
	m->getObjectPtr()->matrix_size[1] = img_dims_[1];
	m->getObjectPtr()->matrix_size[2] = img_dims_[2];
	m->getObjectPtr()->channels       = 1;
	m->getObjectPtr()->image_index    = 1;




	return GADGET_OK;
}

  GADGET_FACTORY_DECLARE(CSIGadget)

} /* namespace Gadgetron */
