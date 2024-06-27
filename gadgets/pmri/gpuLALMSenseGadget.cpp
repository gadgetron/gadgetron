#include "gpuLALMSenseGadget.h"
#include "cuNDArray_operators.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_blas.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "GadgetMRIHeaders.h"
#include "b1_map.h"
#include "GPUTimer.h"
#include "vector_td_utilities.h"
#include "hoNDArray_fileio.h"
#include "ismrmrd/xml.h"
#include <boost/thread/mutex.hpp>
#include "cuNDArray_fileio.h"
#include "eigenTester.h"
namespace Gadgetron{

#define max_number_of_gpus 10
static boost::mutex _mutex[max_number_of_gpus];

gpuLALMSenseGadget::gpuLALMSenseGadget()
: gpuSenseGadget(),
  is_configured_(false)
, prepared_(false)
{
}

gpuLALMSenseGadget::~gpuLALMSenseGadget() {}

int gpuLALMSenseGadget::process_config( ACE_Message_Block* mb )
{
	gpuSenseGadget::process_config(mb);
	number_of_iterations_ = number_of_iterations.value();

	exclusive_access_ = exclusive_access.value();
	lambda_ = lambda.value();
	coils_per_subset_ = coils_per_subset.value();
	huber_value_ = huber_value.value();
	damping_ = damping.value();
	use_preconditioner_ = use_preconditioner.value();

	tau_ = tau.value();
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

	matrix_size_seq_ = uint64d2( r_space.matrixSize.x, r_space.matrixSize.y );

	if (!is_configured_) {

		if (h.acquisitionSystemInformation) {
			channels_ = h.acquisitionSystemInformation->receiverChannels ? *h.acquisitionSystemInformation->receiverChannels : 1;
		} else {
			channels_ = 1;
		}

		// Allocate encoding operator for non-Cartesian Sense
		E_ = boost::make_shared<osSenseOperator<cuNDArray<float_complext>,2,NFFTOperator<cuNDArray,float,2>>>();
		E_->set_coils_per_subset(coils_per_subset_);



		// Setup NLCG solver
		solver_.set_encoding_operator( E_ );

		solver_.set_output_mode( (output_convergence_) ? osLALMSolver<cuNDArray<float_complext>>::OUTPUT_VERBOSE : osLALMSolver<cuNDArray<float_complext>>::OUTPUT_SILENT );
		solver_.set_max_iterations( number_of_iterations_ );

		solver_.set_alpha(huber_value_);
		solver_.set_damping(damping_);
		solver_.set_tau(tau_);
		// Add "TV" regularization
		//

		if( lambda_ > 0.0 ){
			for (auto i = 0u; i < 3; i++){

				auto Dx = boost::make_shared<cuPartialDerivativeOperator<float_complext,3>>(i);
				Dx->set_weight(lambda_);
				TV_ops.push_back(Dx);
			}
			solver_.add_regularization_group(TV_ops);
		}

		// Add "PICCS" regularization
		//


		is_configured_ = true;
	}

	GDEBUG("gpuLALMSenseGadget::end of process_config\n");

	return GADGET_OK;
}

int gpuLALMSenseGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<GenericReconJob> *m2)
{
	// Is this data for this gadget's set/slice?
	//
	GDEBUG("Starting gpuLALMSenseGadget\n");

	if( m1->getObjectPtr()->set != set_number_ || m1->getObjectPtr()->slice != slice_number_ ) {
		// No, pass it downstream...
		return this->next()->putq(m1);
	}

	//GDEBUG("gpuLALMSenseGadget::process\n");
	//GPUTimer timer("gpuLALMSenseGadget::process");

	if (!is_configured_) {
		GDEBUG("\nData received before configuration complete\n");
		return GADGET_FAIL;
	}

	GenericReconJob* j = m2->getObjectPtr();

	// Let's first check that this job has the required data...
	if (!j->csm_host_.get() || !j->dat_host_.get() || !j->tra_host_.get() || !j->dcw_host_.get()) {
		GDEBUG("Received an incomplete Sense job\n");
		return GADGET_FAIL;
	}

	unsigned int samples = j->dat_host_->get_size(0);
	unsigned int channels = j->dat_host_->get_size(1);
	unsigned int rotations = samples / j->tra_host_->get_number_of_elements();
	unsigned int frames = j->tra_host_->get_size(1)*rotations;

	if( samples%j->tra_host_->get_number_of_elements() ) {
		GDEBUG("Mismatch between number of samples (%d) and number of k-space coordinates (%d).\nThe first should be a multiplum of the latter.\n",
				samples, j->tra_host_->get_number_of_elements());
		return GADGET_FAIL;
	}

	boost::shared_ptr< cuNDArray<floatd2> > traj(new cuNDArray<floatd2> (*j->tra_host_));
	boost::shared_ptr< cuNDArray<float> > dcw(new cuNDArray<float> (*j->dcw_host_));
	sqrt_inplace(dcw.get());
	boost::shared_ptr< cuNDArray<float_complext> > csm(new cuNDArray<float_complext> (*j->csm_host_));
	boost::shared_ptr< cuNDArray<float_complext> > device_samples(new cuNDArray<float_complext> (*j->dat_host_));


	// Take the reconstruction matrix size from the regulariaztion image.
	// It could be oversampled from the sequence specified size...

	matrix_size_ = uint64d2( j->reg_host_->get_size(0), j->reg_host_->get_size(1) );

	cudaDeviceProp deviceProp;
	if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
		GDEBUG( "\nError: unable to query device properties.\n" );
		return GADGET_FAIL;
	}

	unsigned int warp_size = deviceProp.warpSize;

	matrix_size_os_ =
			uint64d2(((static_cast<unsigned int>(std::ceil(matrix_size_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
					((static_cast<unsigned int>(std::ceil(matrix_size_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);

	GDEBUG("Matrix size    : [%d,%d] \n", matrix_size_[0], matrix_size_[1]);
	GDEBUG("Matrix size OS : [%d,%d] \n", matrix_size_os_[0], matrix_size_os_[1]);

	std::vector<size_t> image_dims = to_std_vector(matrix_size_);
	image_dims.push_back(frames);

	E_->set_domain_dimensions(&image_dims);
	E_->set_codomain_dimensions(device_samples->get_dimensions().get());
	E_->set_csm(csm);
	E_->setup( matrix_size_, matrix_size_os_, kernel_width_ );
	E_->preprocess(*traj);

	for (auto op : TV_ops){
		op->set_domain_dimensions(&image_dims);
		op->set_codomain_dimensions(&image_dims);
	}


	reg_image_ = boost::shared_ptr< cuNDArray<float_complext> >(new cuNDArray<float_complext>(image_dims));

	// These operators need their domain/codomain set before being added to the solver
	//

	E_->set_dcw(dcw);


	if (!use_preconditioner_){
		GDEBUG("Not using preconditioner, damping %f \n",damping_);
		auto reg = boost::make_shared<cuNDArray<float_complext>>(image_dims);

		fill(reg.get(),float_complext(1));

		solver_.set_preconditioning_image(reg);
	}

	GDEBUG("Prepared\n");

	// Expand the average image to the number of frames
	//


	// Define preconditioning weights
	//

	//Apply weights
	*device_samples *= *dcw;

	// Invoke solver
	//

	boost::shared_ptr< cuNDArray<float_complext> > result;
	{
		GDEBUG("Running NLCG solver\n");
		GPUTimer timer("Running NLCG solver");

		// Optionally, allow exclusive (per device) access to the solver
		// This may not matter much in terms of speed, but it can in terms of memory consumption
		//

		if( exclusive_access_ )
			_mutex[device_number_].lock();

		result = solver_.solve(device_samples.get());

		if( exclusive_access_ )
			_mutex[device_number_].unlock();
	}

	// Provide some info about the scaling between the regularization and reconstruction.
	// If it is not close to one, PICCS does not work optimally...
	//


	if (!result.get()) {
		GDEBUG("\nNon-linear conjugate gradient solver failed\n");
		return GADGET_FAIL;
	}

	// If the recon matrix size exceeds the sequence matrix size then crop
	if( matrix_size_seq_ != matrix_size_ )
		*result = crop<float_complext,2>( (matrix_size_-matrix_size_seq_)>>1, matrix_size_seq_, *result );


	// Now pass on the reconstructed images
	//
	this->put_frames_on_que(frames,rotations,j,result.get(),channels);

	frame_counter_ += frames;
	m1->release();
	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(gpuLALMSenseGadget)
}

