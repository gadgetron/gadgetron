#include "CMRTGadget.h"
#include "cuNFFT.h"
#include "vector_td_utilities.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "permutationOperator.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_operators.h"
#include "radial_utilities.h"
#include "vector_td_operators.h"
#include "cuNFFTOperator.h"
#include "multiplicationOperatorContainer.h"
#include "cuCgSolver.h"

#include "cuNlcgSolver.h"

#include <ismrmrd/xml.h>
#include <cmath>

namespace Gadgetron{


int CMRTGadget::process_config(ACE_Message_Block* mb)
{
	ISMRMRD::IsmrmrdHeader h;
	ISMRMRD::deserialize(mb->rd_ptr(),h);


	if (h.encoding.size() != 1) {
		GADGET_DEBUG1("This Gadget only supports one encoding space\n");
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


	GADGET_DEBUG2("Matrix size: %d, %d, %d\n",
			image_space_dimensions_3D_[0],
			image_space_dimensions_3D_[1],
			image_space_dimensions_3D_[2] );

	GADGET_DEBUG2("Matrix size: %d, %d\n", e_space.matrixSize.x, e_space.matrixSize.y, e_space.matrixSize.z);
	dimensions_.push_back(r_space.matrixSize.x);
	dimensions_.push_back(r_space.matrixSize.y);

	field_of_view_.push_back(e_space.fieldOfView_mm.x);
	field_of_view_.push_back(e_space.fieldOfView_mm.y);
	GADGET_DEBUG2("FOV: %f, %f\n", r_space.fieldOfView_mm.x, r_space.fieldOfView_mm.y);

	repetitions_ = e_limits.repetition.is_present() ? e_limits.repetition.get().maximum + 1 : 1;
	GADGET_DEBUG2("#Repetitions: %d\n", repetitions_);


	// Allocate readout and trajectory/dcw queues
	//

	golden_ratio_ = get_bool_value("golden_ratio");
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
		GADGET_DEBUG2("FRAME # %d \n",num_frames);
		// Get samples for frame
		//
		GADGET_DEBUG1("Extracting samples \n");
		frames.push_back(extract_samples_from_queue( frame_readout_queue_.get()));
		// Get trajectories/dcw for frame - Only for first frame
		//
		if (frames.size() == 1 ){
			extract_trajectory_and_dcw_from_queue( frame_traj_queue_.get(), this->traj, this->dcw);
			GADGET_DEBUG1("Extracting trajectory \n");
		}

		bool is_last_scan_in_slice= m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
		GADGET_DEBUG2("Last scan in slice %i \n",is_last_scan_in_slice);
		//If we have enough projections, get the show on the road
		if (is_last_scan_in_slice){
			num_frames = 0;
			GADGET_DEBUG1("Framing \n");

			boost::shared_ptr<cuNDArray<float_complext> > data = get_combined_frames();
			// Initialize plan
			//
			GADGET_DEBUG2("Data size: %i %i %i",data->get_size(0),data->get_size(1),data->get_size(2));
			boost::shared_ptr<cuNDArray<float> >cu_dcw(new cuNDArray<float>(
					compute_radial_dcw_golden_ratio_2d<float>(
							samples_per_readout_,data->get_size(1),1.0,1.0f/samples_per_readout_/dimensions_[1],0,GR_ORIGINAL).get()));

			sqrt_inplace(cu_dcw.get());
			boost::shared_ptr<cuNDArray<floatd2> >cu_traj(new cuNDArray<floatd2>(*traj));

			std::vector<size_t> projection_dims;
			projection_dims.push_back(dimensions_[0]);
			projection_dims.push_back(dimensions_[1]);
			projection_dims.push_back(dimensions_[3]);


			//cuNDArray<float_complext> result(&image_space_dimensions_3D_);
			boost::shared_ptr<CMRTOperator<float> > E(new CMRTOperator<float>);
			E->setup(data,cu_dcw,cu_traj,image_space_dimensions_3D_,projection_dims,golden_ratio_);

			E->set_domain_dimensions(&image_space_dimensions_3D_);
			E->set_codomain_dimensions(data->get_dimensions().get());

			//cuCgSolver<float_complext> solver;
			cuNlcgSolver<float_complext> solver;
			solver.set_encoding_operator(E);
			solver.set_max_iterations(20);
			solver.set_tc_tolerance(1e-8f);

			solver.set_output_mode(cuCgSolver<float_complext>::OUTPUT_VERBOSE);

			//*data *= *cu_dcw;

			boost::shared_ptr<cuNDArray<float_complext> > result = solver.solve(data.get());
			//boost::shared_ptr<cuNDArray<float_complext> > result(new cuNDArray<float_complext>(&image_space_dimensions_3D_));
			//E->mult_MH(data.get(),result.get());
			GADGET_DEBUG1(" Penguins report mission accomplished \n");



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
			boost::shared_ptr< hoNDArray< complext<float> > > host_result = result->to_host();
			*cm2->getObjectPtr() = *((hoNDArray< std::complex<float> >*) host_result.get());

			cm1->getObjectPtr()->matrix_size[0] = image_space_dimensions_3D_[0];
			cm1->getObjectPtr()->matrix_size[1] = image_space_dimensions_3D_[1];
			cm1->getObjectPtr()->matrix_size[2] = image_space_dimensions_3D_[2];
			cm1->getObjectPtr()->channels       = 1;
			cm1->getObjectPtr()->image_index    = 1;

			if (this->next()->putq(cm1) < 0) {
				GADGET_DEBUG1("Failed to put result image on to queue\n");
				cm1->release();
				return GADGET_FAIL;
			}
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
		GADGET_DEBUG1("Illegal queue pointer, cannot extract samples\n");
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
			GADGET_DEBUG1("Message dequeue failed\n");
			throw std::runtime_error("CMRTGadget::extract_samples_from_queue: dequeing failed");
		}

		GadgetContainerMessage< hoNDArray< std::complex<float> > > *daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);

		if (!daq) {
			GADGET_DEBUG1("Unable to interpret data on message queue\n");
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
		GADGET_DEBUG1("Illegal queue pointer, cannot extract trajectory\n");
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
			GADGET_DEBUG1("Message dequeue failed\n");
			throw std::runtime_error("CMRTGadget::extract_trajectory_from_queue: dequeing failed");
		}

		GadgetContainerMessage< hoNDArray<float> > *daq = AsContainerMessage<hoNDArray<float> >(mbq);

		if (!daq) {
			GADGET_DEBUG1("Unable to interpret data on message queue\n");
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

		boost::shared_ptr< hoNDArray<float> > host_traj_dcw_shifted = permute( host_traj_dcw.get(), &order );

		std::vector<size_t> dims_1d;
		dims_1d.push_back(host_traj_dcw_shifted->get_size(0)*host_traj_dcw_shifted->get_size(1));

		dcw = boost::shared_ptr<hoNDArray<float> > (new hoNDArray<float>(&dims_1d, host_traj_dcw_shifted->get_data_ptr()+2*dims_1d[0]));


		std::vector<size_t> dims_2d = dims_1d; dims_2d.push_back(2);
		order.clear(); order.push_back(1); order.push_back(0);


		hoNDArray<float> tmp(&dims_2d, host_traj_dcw_shifted->get_data_ptr());

		boost::shared_ptr< hoNDArray<float> > _traj = permute( &tmp, &order );

		traj = boost::shared_ptr<hoNDArray<floatd2> > (new hoNDArray<floatd2>(&dims_1d, (floatd2*)_traj->get_data_ptr()));
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
