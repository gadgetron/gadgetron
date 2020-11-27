#include "CMRTGadget.h"
#include "cuNFFT.h"
#include "vector_td_utilities.h"
#include "permutationOperator.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "radial_utilities.h"
#include "vector_td_operators.h"
#include "../../toolboxes/nfft/NFFTOperator.h"
#include "multiplicationOperatorContainer.h"
#include "cuCgSolver.h"
#include "cuTvOperator.h"
#include "lbfgsSolver.h"
#include "cuSbcCgSolver.h"
#include "cuPartialDerivativeOperator.h"
#include "cuPartialDerivativeOperator2.h"
#include <numeric>
#include <functional>
#include "cuNlcgSolver.h"
#include <boost/make_shared.hpp>

#include <ismrmrd/xml.h>
#include <cmath>

namespace Gadgetron{


int CMRTGadget::process_config(ACE_Message_Block* mb)
{
	ISMRMRD::IsmrmrdHeader h;
	ISMRMRD::deserialize(mb->rd_ptr(),h);


	if (h.encoding.size() != 1) {
		GDEBUG("This Gadget only supports one encoding space\n");
		return GADGET_FAIL;
	}

	ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
	ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
	ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;


	// Matrix size x is the oversampled readout size : size FH*2 -- not for the hyperplarization?
	// Matrix size y is the AP/RL size

	image_space_dimensions_3D_.push_back(e_space.matrixSize.y);
	image_space_dimensions_3D_.push_back(e_space.matrixSize.y);
	image_space_dimensions_3D_.push_back(e_space.matrixSize.x/*/2*/);


	GDEBUG("Matrix size: %d, %d, %d\n",
			image_space_dimensions_3D_[0],
			image_space_dimensions_3D_[1],
			image_space_dimensions_3D_[2] );

	GDEBUG("Matrix size: %d, %d\n", e_space.matrixSize.x, e_space.matrixSize.y, e_space.matrixSize.z);
	dimensions_.push_back(r_space.matrixSize.x);
	dimensions_.push_back(r_space.matrixSize.y);

	field_of_view_.push_back(e_space.fieldOfView_mm.x);
	field_of_view_.push_back(e_space.fieldOfView_mm.y);
	GDEBUG("FOV: %f, %f\n", r_space.fieldOfView_mm.x, r_space.fieldOfView_mm.y);

	repetitions_ = e_limits.repetition.is_present() ? e_limits.repetition.get().maximum + 1 : 1;
	GDEBUG("#Repetitions: %d\n", repetitions_);


	// Allocate readout and trajectory/dcw queues
	//

	golden_ratio_ =golden_ratio.value();
	use_TV_ = use_TV.value();
	projections_per_recon_ = projections_per_recon.value();
	iterations_ = iterations.value();
	frame_readout_queue_ = boost::shared_ptr< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>());
	frame_traj_queue_ = boost::shared_ptr< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>());

	size_t bsize = sizeof(GadgetContainerMessage< hoNDArray< std::complex<float> > >)*dimensions_[0]*10;

	frame_readout_queue_->high_water_mark(bsize);
	frame_readout_queue_->low_water_mark(bsize);
	frame_traj_queue_->high_water_mark(bsize);
	frame_traj_queue_->low_water_mark(bsize);

	return GADGET_OK;
}

int CMRTGadget::process(GadgetContainerMessage< ISMRMRD::AcquisitionHeader > *m1,        // header
		GadgetContainerMessage< hoNDArray< std::complex<float> > > *m2,  // data
		GadgetContainerMessage< hoNDArray<float> > *m3 )                 // traj/dcw
{
	// Throw away any noise samples if they have been allowed to pass this far down the chain...
	//

	bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
	if (is_noise) {
		m1->release();
		return GADGET_OK;
	}

	// First pass initialization
	//

	if (frame_readout_queue_->message_count() == 0 ) {
		samples_per_readout_ = m1->getObjectPtr()->number_of_samples;
		num_coils_ = m1->getObjectPtr()->active_channels;
		dimensions_.push_back(m1->getObjectPtr()->active_channels);
		dimensions_.push_back(repetitions_);
		num_trajectory_dims_ = m3->getObjectPtr()->get_size(0); // 2 for trajectories only, 3 for both trajectories + dcw
	}

	int samples = m1->getObjectPtr()->number_of_samples;
	int readout = m1->getObjectPtr()->idx.kspace_encode_step_1;
	int repetition = m1->getObjectPtr()->idx.kspace_encode_step_2;

	// Enqueue incoming readouts and trajectories
	//

	frame_readout_queue_->enqueue_tail(duplicate_array(m2));

	//Only extract trajectories for first frame. Assume next frames are equal
	if (frames.size() == 0 )
		frame_traj_queue_->enqueue_tail(duplicate_array(m3));

	// If the last readout for a slice has arrived then perform a reconstruction
	//

	bool is_last_scan_in_repetition =
			m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_REPETITION);

	if (is_last_scan_in_repetition) {
		num_frames++;
		GDEBUG("FRAME # %d \n",num_frames);
		// Get samples for frame
		//
		GDEBUG("Extracting samples \n");
		frames.push_back(extract_samples_from_queue( frame_readout_queue_.get()));
		// Get trajectories/dcw for frame - Only for first frame
		//
		if (frames.size() == 1 ){
			extract_trajectory_and_dcw_from_queue( frame_traj_queue_.get(), this->traj, this->dcw);
			GDEBUG("Extracting trajectory \n");
		}

		bool is_last_scan_in_slice= m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
		GDEBUG("Last scan in slice %i \n",is_last_scan_in_slice);
		//If we have enough projections, get the show on the road
		if (is_last_scan_in_slice ){
			GDEBUG("Framing %i frames \n",num_frames);

			if (num_frames%projections_per_recon_ != 0) {
				GDEBUG("Number of frames must be divisible by number of projecitons");
				return GADGET_FAIL;
			}
			boost::shared_ptr<cuNDArray<float_complext> > data = get_combined_frames();
			auto data_dims = data->get_dimensions();
			size_t ntimeframes = data_dims->back()/projections_per_recon_;
			data_dims->back() = projections_per_recon_;
			data_dims->push_back(ntimeframes);
			data->reshape(data_dims);
			// Initialize plan
			//
			GDEBUG("Data size: %i %i %i",data->get_size(0),data->get_size(1),data->get_size(2));
			boost::shared_ptr<cuNDArray<floatd2> >cu_traj(new cuNDArray<floatd2>(*traj));

			std::vector<size_t> projection_dims;
			projection_dims.push_back(dimensions_[0]*2);
			projection_dims.push_back(dimensions_[1]);

			projection_dims.push_back(projections_per_recon_);
			projection_dims.push_back(ntimeframes);


			//cuNDArray<float_complext> result(&image_space_dimensions_3D_);
			boost::shared_ptr<CMRTOperator<float> > E(new CMRTOperator<float>);
			E->setup(cu_traj,image_space_dimensions_3D_,projection_dims,0,golden_ratio_);

			auto image_space_dimensions_4D = image_space_dimensions_3D_;
			image_space_dimensions_4D.push_back(ntimeframes);
			E->set_domain_dimensions(&image_space_dimensions_4D);
			E->set_codomain_dimensions(data->get_dimensions().get());


			boost::shared_ptr<cuNDArray<float_complext> > result;
			//cuCgSolver<float_complext> solver;
			//cuNlcgSolver<float_complext> solver;

			if (use_TV_){
				cuSbcCgSolver<float_complext> solver;
				solver.set_encoding_operator(E);
				//solver.set_max_iterations(20);
				solver.set_max_outer_iterations(iterations_);
				solver.get_inner_solver()->set_max_iterations(10);
				solver.set_tc_tolerance(1e-8f);
				auto Rx1_ = boost::make_shared< cuPartialDerivativeOperator<float_complext,4> >(0);

				auto Ry1_ = boost::make_shared< cuPartialDerivativeOperator<float_complext,4> >(1);
				auto Rz1_ = boost::make_shared< cuPartialDerivativeOperator<float_complext,4> >(2);

				auto Rt1_ = boost::make_shared< cuPartialDerivativeOperator2<float_complext,4> >();

				Rx1_->set_domain_dimensions(&image_space_dimensions_4D);
				Rx1_->set_codomain_dimensions(&image_space_dimensions_4D);

				Ry1_->set_domain_dimensions(&image_space_dimensions_4D);
				Ry1_->set_codomain_dimensions(&image_space_dimensions_4D);

				Rz1_->set_domain_dimensions(&image_space_dimensions_4D);
				Rz1_->set_codomain_dimensions(&image_space_dimensions_4D);

				Rt1_->set_domain_dimensions(&image_space_dimensions_4D);
				Rt1_->set_codomain_dimensions(&image_space_dimensions_4D);
				float lambda = 2000;
				float mu = 1000;
				Rx1_ ->set_weight(lambda);
				Ry1_ ->set_weight(lambda);
				Rz1_ ->set_weight(lambda);
				Rt1_->set_weight(lambda);
				E->set_weight(mu);
				solver.add_regularization_group_operator(Rx1_);
				solver.add_regularization_group_operator(Ry1_);
				solver.add_regularization_group_operator(Rz1_);
				solver.add_regularization_group_operator(Rt1_);
				solver.add_group();




				solver.set_output_mode(cuCgSolver<float_complext>::OUTPUT_VERBOSE);

				//*data *= *cu_dcw;

				result = solver.solve(data.get());
			} else
			{
				cuCgSolver<float_complext> solver;
				//cuNlcgSolver<float_complext> solver;
				solver.set_encoding_operator(E);
				solver.set_max_iterations(iterations_);
				solver.set_tc_tolerance(1e-8f);
				solver.set_output_mode(cuCgSolver<float_complext>::OUTPUT_VERBOSE);

				result = solver.solve(data.get());
			}
			//boost::shared_ptr<cuNDArray<float_complext> > result(new cuNDArray<float_complext>(&image_space_dimensions_3D_));
			//E->mult_MH(data.get(),result.get());
			GDEBUG(" Penguins report mission accomplished \n");


			size_t nelements3d = std::accumulate(image_space_dimensions_3D_.begin(),image_space_dimensions_3D_.end(),1,std::multiplies<size_t>());

			for (size_t i = 0; i < ntimeframes; i++){

				// Define the image header
				//

				GadgetContainerMessage<ISMRMRD::ImageHeader> *cm1 =
						new GadgetContainerMessage<ISMRMRD::ImageHeader>();

				GadgetContainerMessage< hoNDArray< std::complex<float> > > *cm2 =
						new GadgetContainerMessage<hoNDArray< std::complex<float> > >();

				cm1->getObjectPtr()->flags = 0;
				cm1->cont(cm2);


				cm1->getObjectPtr()->field_of_view[0]   = field_of_view_[0];
				cm1->getObjectPtr()->field_of_view[1]   = field_of_view_[1];
				cm1->getObjectPtr()->channels           = num_coils_;
				cm1->getObjectPtr()->repetition         = m1->getObjectPtr()->idx.repetition;



				memcpy(cm1->getObjectPtr()->patient_table_position,
						m1->getObjectPtr()->patient_table_position, sizeof(float)*3);

				cm1->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;
				cm1->getObjectPtr()->image_index = 0;
				cm1->getObjectPtr()->image_series_index = 0;

				// std::complex<float> and Gadgetron::complext<float> are binary compatible
				cuNDArray<complext<float> > cuView(image_space_dimensions_3D_,result->get_data_ptr()+i*nelements3d);
				boost::shared_ptr< hoNDArray< complext<float> > > host_result = cuView.to_host();
				*cm2->getObjectPtr() = *((hoNDArray< std::complex<float> >*) host_result.get());

				cm1->getObjectPtr()->matrix_size[0] = image_space_dimensions_3D_[0];
				cm1->getObjectPtr()->matrix_size[1] = image_space_dimensions_3D_[1];
				cm1->getObjectPtr()->matrix_size[2] = image_space_dimensions_3D_[2];
				cm1->getObjectPtr()->channels       = 1;
				cm1->getObjectPtr()->image_index    = i+1;

				if (this->next()->putq(cm1) < 0) {
					GDEBUG("Failed to put result image on to queue\n");
					cm1->release();
					return GADGET_FAIL;
				}
			}

			num_frames = 0;
		}

	}

	m1->release();
	return GADGET_OK;
}

template<class T> GadgetContainerMessage< hoNDArray<T> >*
CMRTGadget::duplicate_array( GadgetContainerMessage< hoNDArray<T> > *array )
{
	GadgetContainerMessage< hoNDArray<T> > *copy = new GadgetContainerMessage< hoNDArray<T> >();
	*(copy->getObjectPtr()) = *(array->getObjectPtr());
	return copy;
}

boost::shared_ptr< hoNDArray<float_complext> >
CMRTGadget::extract_samples_from_queue ( ACE_Message_Queue<ACE_MT_SYNCH> *queue )
{
	if(!queue) {
		GDEBUG("Illegal queue pointer, cannot extract samples\n");
		throw std::runtime_error("CMRTGadget::extract_samples_from_queue: illegal queue pointer");
	}

	unsigned int readouts_buffered = queue->message_count();

	std::vector<size_t> dims;
	dims.push_back(samples_per_readout_);
	dims.push_back(readouts_buffered);
	dims.push_back(num_coils_);

	boost::shared_ptr< hoNDArray<float_complext> > host_samples(new hoNDArray<float_complext>(dims));

	for (unsigned int p=0; p<readouts_buffered; p++) {

		ACE_Message_Block* mbq;
		if (queue->dequeue_head(mbq) < 0) {
			GDEBUG("Message dequeue failed\n");
			throw std::runtime_error("CMRTGadget::extract_samples_from_queue: dequeing failed");
		}

		GadgetContainerMessage< hoNDArray< std::complex<float> > > *daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);

		if (!daq) {
			GDEBUG("Unable to interpret data on message queue\n");
			throw std::runtime_error("CMRTGadget::extract_samples_from_queue: failed to interpret data");
		}

		for (unsigned int c = 0; c < num_coils_; c++) {

			float_complext *data_ptr = host_samples->get_data_ptr();
			data_ptr += c*samples_per_readout_*readouts_buffered+p*samples_per_readout_;

			std::complex<float> *r_ptr = daq->getObjectPtr()->get_data_ptr();
			r_ptr += c*daq->getObjectPtr()->get_size(0);

			memcpy(data_ptr, r_ptr, samples_per_readout_*sizeof(float_complext));
		}

		mbq->release();
	}

	return host_samples;
}

boost::shared_ptr< hoNDArray<float> >
CMRTGadget::extract_trajectory_from_queue ( ACE_Message_Queue<ACE_MT_SYNCH> *queue )
{
	if(!queue) {
		GDEBUG("Illegal queue pointer, cannot extract trajectory\n");
		throw std::runtime_error("CMRTGadget::extract_trajectory_from_queue: illegal queue pointer");
	}

	unsigned int readouts_buffered = queue->message_count();

	std::vector<size_t> dims;
	dims.push_back(num_trajectory_dims_); // 2 for trajectories only, 3 for both trajectories + dcw
	dims.push_back(samples_per_readout_);
	dims.push_back(readouts_buffered);

	boost::shared_ptr< hoNDArray<float> > host_traj(new hoNDArray<float>(&dims));

	for (unsigned int p=0; p<readouts_buffered; p++) {
		ACE_Message_Block* mbq;
		if (queue->dequeue_head(mbq) < 0) {
			GDEBUG("Message dequeue failed\n");
			throw std::runtime_error("CMRTGadget::extract_trajectory_from_queue: dequeing failed");
		}

		GadgetContainerMessage< hoNDArray<float> > *daq = AsContainerMessage<hoNDArray<float> >(mbq);

		if (!daq) {
			GDEBUG("Unable to interpret data on message queue\n");
			throw std::runtime_error("CMRTGadget::extract_trajectory_from_queue: failed to interpret data");
		}

		float *data_ptr = host_traj->get_data_ptr();
		data_ptr += num_trajectory_dims_*samples_per_readout_*p;

		float *r_ptr = daq->getObjectPtr()->get_data_ptr();

		memcpy(data_ptr, r_ptr, num_trajectory_dims_*samples_per_readout_*sizeof(float));

		mbq->release();
	}

	return host_traj;
}

void CMRTGadget::extract_trajectory_and_dcw_from_queue
( ACE_Message_Queue<ACE_MT_SYNCH> *queue, boost::shared_ptr< hoNDArray<floatd2> > & traj, boost::shared_ptr< hoNDArray<float> > & dcw )
{
	// Extract trajectory and (if present) density compensation weights.
	// They are stored as a float array of dimensions: {2,3} x #samples_per_readout x #readouts.
	// We need
	// - a floatd2 trajectory array
	// - a float dcw array
	//

	if( num_trajectory_dims_ == 2 ){
		//This is an evil evil hack to get the trajectories out. Ohh the horror.
		boost::shared_ptr<hoNDArray<float> > tmp_traj = extract_trajectory_from_queue( queue );
		std::vector<size_t> dims_1d; dims_1d.push_back(tmp_traj->get_size(1)*tmp_traj->get_size(2));
		traj = boost::shared_ptr<hoNDArray<floatd2> >(new hoNDArray<floatd2>(&dims_1d));
		memcpy(traj->get_data_ptr(),tmp_traj->get_data_ptr(),tmp_traj->get_number_of_elements()*sizeof(float));


	}
	else{

		boost::shared_ptr< hoNDArray<float> > host_traj_dcw = extract_trajectory_from_queue( queue );

		std::vector<size_t> order;
		order.push_back(1); order.push_back(2); order.push_back(0);

		auto host_traj_dcw_shifted = permute( *host_traj_dcw, order );

		std::vector<size_t> dims_1d;
		dims_1d.push_back(host_traj_dcw_shifted.get_size(0)*host_traj_dcw_shifted.get_size(1));

		dcw = boost::shared_ptr<hoNDArray<float> > (new hoNDArray<float>(&dims_1d, host_traj_dcw_shifted.get_data_ptr()+2*dims_1d[0]));


		std::vector<size_t> dims_2d = dims_1d; dims_2d.push_back(2);
		order.clear(); order.push_back(1); order.push_back(0);


		hoNDArray<float> tmp(&dims_2d, host_traj_dcw_shifted.get_data_ptr());

		auto _traj = permute( tmp, order );

		traj = boost::shared_ptr<hoNDArray<floatd2> > (new hoNDArray<floatd2>(&dims_1d, (floatd2*)_traj.get_data_ptr()));
	}

	std::vector<size_t >dims_2d;
	dims_2d.push_back(traj->get_number_of_elements());
	dims_2d.push_back(1); // Number of frames

	traj->reshape(&dims_2d);
	if( num_trajectory_dims_ == 3 ) dcw->reshape(&dims_2d);
}

boost::shared_ptr<cuNDArray<float_complext> > CMRTGadget::get_combined_frames(){
	if (frames.size() == 0)
		throw std::runtime_error("No frames received. This should not be possible. Your RAM might be replaced with live salmon, or you may have set the expected number of frames to 0");

	for (unsigned int i = 1; i < frames.size(); i++){
		if (!frames[0]->dimensions_equal(frames[i].get()))
			throw std::runtime_error("CMRTGadget: Frames received do not have equal size");
	}
	//Get data dimensions. Assume all frames have the same dimensions
	std::vector<size_t> dims = *frames[0]->get_dimensions();
	dims.push_back(frames.size());

	boost::shared_ptr<cuNDArray<float_complext> > combined(new cuNDArray<float_complext>(dims));

	//Copy data into 1 array on device
	size_t offset = 0;
	for (unsigned int i = 0; i < frames.size(); i++){
		cudaMemcpy(combined->get_data_ptr()+offset,frames[i]->get_data_ptr(),frames[i]->get_number_of_elements()*sizeof(float_complext),cudaMemcpyHostToDevice);
		offset += frames[i]->get_number_of_elements();
	}

	frames.clear();
	return combined;






}



GADGET_FACTORY_DECLARE(CMRTGadget)
}
