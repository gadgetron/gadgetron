/*
 * CSIGadget.cpp
 *
 *  Created on: Nov 11, 2014
 *      Author: dch
 */

#include "CSIGadget.h"
#include <ismrmrd/xml.h>
#include "cudaDeviceManager.h"
#include "cuNDArray_utils.h"
#include "cuNlcgSolver.h"
#include "eigenTester.h"
#include "CSfreqOperator.h"
#include "cuPartialDerivativeOperator.h"
#include "cuDWTOperator.h"
#include <boost/make_shared.hpp>
namespace Gadgetron {

CSIGadget::CSIGadget() {
	// TODO Auto-generated constructor stub

}

CSIGadget::~CSIGadget() {
	// TODO Auto-generated destructor stub
}


int CSIGadget::process_config(ACE_Message_Block *mb){
	//GDEBUG("gpuCgSenseGadget::process_config\n");

        device_number_ = deviceno.value();

	int number_of_devices = 0;
	if (cudaGetDeviceCount(&number_of_devices)!= cudaSuccess) {
		GDEBUG( "Error: unable to query number of CUDA devices.\n" );
		return GADGET_FAIL;
	}

	if (number_of_devices == 0) {
		GDEBUG( "Error: No available CUDA devices.\n" );
		return GADGET_FAIL;
	}

	if (device_number_ >= number_of_devices) {
		GDEBUG("Adjusting device number from %d to %d\n", device_number_,  (device_number_%number_of_devices));
		device_number_ = (device_number_%number_of_devices);
	}

	if (cudaSetDevice(device_number_)!= cudaSuccess) {
		GDEBUG( "Error: unable to set CUDA device.\n" );
		return GADGET_FAIL;
	}

	cg_limit_ = cg_limit.value();
	oversampling_factor_ = oversampling_factor.value();
	kernel_width_ = kernel_width.value();
	output_convergence_ = output_convergence.value();
	number_of_sb_iterations_ = number_of_sb_iterations.value();
	number_of_cg_iterations_ = number_of_cg_iterations.value();

	use_compressed_sensing_ = compressed_sensing.value();

	mu_ = mu.value();

	// Get the Ismrmrd header
	//
	ISMRMRD::IsmrmrdHeader h;
	ISMRMRD::deserialize(mb->rd_ptr(),h);



	if (h.encoding.size() != 1) {
		GDEBUG("This Gadget only supports one encoding space\n");
		return GADGET_FAIL;
	}

	// Get the encoding space and trajectory description
	ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
	ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
	ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;

	img_dims_ = {r_space.matrixSize.x,r_space.matrixSize.y,r_space.matrixSize.z};

	matrix_size_ = vector_td<uint64_t,2>{r_space.matrixSize.x,r_space.matrixSize.y};

	unsigned int warp_size = cudaDeviceManager::Instance()->warp_size(device_number_);

	matrix_size_os_ =
			uint64d2(((static_cast<unsigned int>(std::ceil(matrix_size_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
					((static_cast<unsigned int>(std::ceil(matrix_size_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);


		if (h.acquisitionSystemInformation) {
			channels_ = h.acquisitionSystemInformation->receiverChannels ? *h.acquisitionSystemInformation->receiverChannels : 1;
		} else {
			channels_ = 1;
		}


		/*if (~h.userParameters.is_present()){
			GDEBUG("CSI gadget requires userparameters to be set to obtain timesteps.");
			return GADGET_FAIL;
		}*/


		auto parameters = h.userParameters->userParameterDouble;
		auto bw = std::find_if(parameters.begin(),parameters.end(), [](ISMRMRD::UserParameterDouble d) { return d.name=="bw";});

		if (bw->name != "bw"){
			GDEBUG("CSI gadget: User parameter bw is missing.");
			return GADGET_FAIL;
		}


		auto dte = std::find_if(parameters.begin(),parameters.end(), [](ISMRMRD::UserParameterDouble d) { return d.name=="dte";});
		if (dte->name != "dte"){
			GDEBUG("CSI gadget: User parameter dte is missing.");
			return GADGET_FAIL;
		}

// Allocate encoding operator for non-Cartesian Sense
		E_ = boost::make_shared< CSIOperator<float> >(1/bw->value, dte->value);
		E_->set_weight(mu_);

		std::vector<float> freqs = frequencies.value();

		if (freqs.empty()){
			for (float f = frequency_min.value(); f <= frequency_max.value(); f+= frequency_step.value())
				freqs.push_back(f);
		}


		if (freqs.size() == 0)
			throw std::runtime_error("CSIGadget: Frequencies not set!");

		std::stringstream ss;
		ss << "Frequencies set: ";
		for (auto f : freqs)
				ss << f << " ";
		GDEBUG(ss.str().c_str());

		E_->set_frequencies(freqs);

		img_dims_[2]=freqs.size();

		S_ = boost::make_shared<cuNonCartesianSenseOperator<float,2>>();

		E_->set_senseOp(S_);


		auto idOp = boost::make_shared<identityOperator<cuNDArray<float_complext>>>();
		idOp->set_domain_dimensions(&img_dims_);
		idOp->set_codomain_dimensions(&img_dims_);
		idOp->set_weight(2*mu_);
		solver_.add_regularization_operator(idOp);

		auto dX = boost::make_shared<cuPartialDerivativeOperator<float_complext,3>>(0);
		dX->set_domain_dimensions(&img_dims_);
		dX->set_codomain_dimensions(&img_dims_);
		dX->set_weight(2*mu_);
		auto dY = boost::make_shared<cuPartialDerivativeOperator<float_complext,3>>(1);
		dY->set_domain_dimensions(&img_dims_);
		dY->set_codomain_dimensions(&img_dims_);
		dY->set_weight(2*mu_);
		auto dZ = boost::make_shared<cuPartialDerivativeOperator<float_complext,3>>(2);
		dZ->set_domain_dimensions(&img_dims_);
		dZ->set_codomain_dimensions(&img_dims_);
		dZ->set_weight(2*mu_);


		solver_.add_regularization_group_operator(dX);
		solver_.add_regularization_group_operator(dY);

		//solver_.add_regularization_group_operator(dZ);
		solver_.add_group();
/*
		auto W = boost::make_shared<cuDWTOperator<float_complext,3>>();
		W->set_domain_dimensions(&img_dims_);
		W->set_codomain_dimensions(&img_dims_);
		W->set_weight(2*mu_);


		auto W2 = boost::make_shared<cuDWTOperator<float_complext,3>>();
		W2->set_shift(2);
		W2->set_domain_dimensions(&img_dims_);
		W2->set_codomain_dimensions(&img_dims_);
		W2->set_weight(2*mu_);
		solver_.add_regularization_operator(W);
		solver_.add_regularization_operator(W2);
		*/
		// Setup solver
		solver_.set_encoding_operator( E_ );        // encoding matrix
		solver_.set_max_outer_iterations( number_of_sb_iterations_ );
		solver_.set_max_inner_iterations(1);
		solver_.get_inner_solver()->set_max_iterations(number_of_cg_iterations_);
		solver_.get_inner_solver()->set_tc_tolerance( cg_limit_ );
		solver_.set_output_mode( (output_convergence_) ? cuCgSolver<float_complext>::OUTPUT_VERBOSE : cuCgSolver<float_complext>::OUTPUT_SILENT);
		is_configured_ = true;
	return GADGET_OK;

}

int CSIGadget::process(GadgetContainerMessage<cuSenseData>* m1){


	if (!is_configured_) {
		GDEBUG("\nData received before configuration complete\n");
		return GADGET_FAIL;
	}


	GDEBUG("CSI is on the job\n");



	auto traj = m1->getObjectPtr()->traj;

	auto trajdims2 = std::vector<size_t>{traj->get_size(0),1};
	//Extract initial trajectory
	cuNDArray<floatd2> traj2(trajdims2,traj->get_data_ptr());

	auto data = m1->getObjectPtr()->data;
	auto csm =m1->getObjectPtr()->csm;
	auto dcw = m1->getObjectPtr()->dcw;
	//dcw.reset();




	if (dcw)
		sqrt_inplace(dcw.get());


	E_->set_domain_dimensions(&img_dims_);
	E_->set_codomain_dimensions(data->get_dimensions().get());

	std::vector<size_t> sense_dims = *data->get_dimensions();
	sense_dims[1] = img_dims_[2];



	S_->set_domain_dimensions(&img_dims_);
	S_->set_codomain_dimensions(&sense_dims);
/*
	{
		GDEBUG("Removing CSM maps");
		auto csm_dims = *csm->get_dimensions();
		csm_dims.pop_back();
		cuNDArray<float_complext> csm_view(csm_dims,csm->get_data_ptr());
		fill(&csm_view,complext<float>(1,0));
		size_t nelements = csm_view.get_number_of_elements();
		for (int  i = 1; i< csm->get_size(2); i++){
			cuNDArray<float_complext> csm_view2(csm_dims,csm->get_data_ptr()+i*nelements);
			clear(&csm_view2);
		}
	}
*/
	S_->set_csm(csm);
	S_->set_dcw(dcw);
	S_->setup( matrix_size_, matrix_size_os_, kernel_width_ );
	S_->preprocess(&traj2);

	GDEBUG("Setup done, solving....\n");
	/*
	eigenTester<cuNDArray<float_complext>> tester;
	std::vector<float> freqs{  -575.1223,-450.1223,-360.1223,  -183.1223,140.8777};
	auto T_ = boost::make_shared<CSfreqOperator>(E_->get_pointtime(),E_->get_echotime());
	T_->set_frequencies(freqs);
	T_->set_codomain_dimensions(data->get_dimensions().get());

	std::vector<size_t> tim_dims = *data->get_dimensions();
	tim_dims[2] = freqs.size();
	T_->set_domain_dimensions(&tim_dims);

	tester.add_encoding_operator(E_);

	float_complext eigVal = tester.get_smallest_eigenvalue();

	GDEBUG("Smallest eigenvalue: %f %f /n",real(eigVal),imag(eigVal));
*/
	/*
	cuNlcgSolver<float_complext> solv;
	//cuCgSolver<float_complext> solv;
	solv.set_output_mode(cuCgSolver<float_complext>::OUTPUT_VERBOSE);
	solv.set_max_iterations(10);
	solv.set_encoding_operator(E_);
	solv.set_tc_tolerance(1e-8f);
	*/
	boost::shared_ptr<cuNDArray<float_complext>> result;
	if (use_compressed_sensing_)
		result = solver_.solve(data.get());
	else {
		cgSolver<cuNDArray<float_complext>> cgsolver;
		cgsolver.set_max_iterations(solver_.get_inner_solver()->get_max_iterations());
		cgsolver.set_encoding_operator(E_);
		result = cgsolver.solve(data.get());
	}
	//auto result = solv.solve(data.get());

	//E_->mult_MH(data.get(),result.get(),false);

	GDEBUG("Image sum: %f \n",asum(result.get()));
	m1->release();

	GDEBUG("Solver done, next patient...");

	GadgetContainerMessage< hoNDArray< std::complex<float> > > *cm =
			new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

	GadgetContainerMessage<ISMRMRD::ImageHeader> *m =
			new GadgetContainerMessage<ISMRMRD::ImageHeader>();


	m->cont(cm);



	result->to_host((hoNDArray<float_complext>*)cm->getObjectPtr());

	GDEBUG("Result size: %i %i %i \n",result->get_size(0),result->get_size(1),result->get_size(2));

	m->getObjectPtr()->matrix_size[0] = img_dims_[0];
	m->getObjectPtr()->matrix_size[1] = img_dims_[1];
	m->getObjectPtr()->matrix_size[2] = img_dims_[2];
	m->getObjectPtr()->channels       = 1;
	m->getObjectPtr()->image_index    = 1;
	m->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;


	if (!this->next()->putq(m)){
		GDEBUG("Failed to put image on que");
		return GADGET_FAIL;
	}


	return GADGET_OK;
}

  GADGET_FACTORY_DECLARE(CSIGadget)

} /* namespace Gadgetron */
