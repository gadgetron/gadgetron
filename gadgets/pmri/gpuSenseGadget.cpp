/*
 * gpuSenseGadget.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: dch
 */

#include "gpuSenseGadget.h"
#include "cuNDArray.h"
#include "vector_td_utilities.h"
namespace Gadgetron {

gpuSenseGadget::gpuSenseGadget() {
	set_parameter(std::string("deviceno").c_str(), "0");
	set_parameter(std::string("setno").c_str(), "0");
	set_parameter(std::string("sliceno").c_str(), "0");
	set_parameter(std::string("cg_limit").c_str(), "1e-6");
	set_parameter(std::string("oversampling_factor").c_str(), "1.5");
	set_parameter(std::string("kernel_width").c_str(), "5.5");
	set_parameter(std::string("save_individual_frames").c_str(),"true");

	matrix_size_ = uint64d2(0u,0u);
	matrix_size_os_ = uint64d2(0u,0u);
	matrix_size_seq_ = uint64d2(0u,0u);
}

gpuSenseGadget::~gpuSenseGadget() {
	// TODO Auto-generated destructor stub
}

int gpuSenseGadget::process_config(ACE_Message_Block* mb) {
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
	set_number_ = get_int_value(std::string("setno").c_str());
	slice_number_ = get_int_value(std::string("sliceno").c_str());
	oversampling_factor_ = get_double_value(std::string("oversampling_factor").c_str());
	kernel_width_ = get_double_value(std::string("kernel_width").c_str());
	output_convergence_ = get_bool_value(std::string("output_convergence").c_str());
	output_timing_ = get_bool_value(std::string("output_timing").c_str());
	rotations_to_discard_ = get_int_value(std::string("rotations_to_discard").c_str());

	if( (rotations_to_discard_%2) == 1 ){
		GADGET_DEBUG1("#rotations to discard must be even.\n");
		return GADGET_FAIL;
	}
	save_individual_frames_ = get_bool_value("save_individual_frames");


}

int gpuSenseGadget::put_frames_on_que(int frames,int rotations, GenericReconJob* j, cuNDArray<float_complext>* cgresult) {

	unsigned int frames_per_rotation = frames/rotations;

	if( rotations == 1 ){ // this is the case for golden ratio
		rotations = frames;
		frames_per_rotation = 1;
	}
	if (save_individual_frames_){
		for( unsigned int frame=0; frame<frames; frame++ ){

			unsigned int rotation_idx = frame/frames_per_rotation;

			// Check if we should discard this frame
			if( rotation_idx < (rotations_to_discard_>>1) || rotation_idx >= rotations-(rotations_to_discard_>>1) )
				continue;

			GadgetContainerMessage<ISMRMRD::ImageHeader> *m =
					new GadgetContainerMessage<ISMRMRD::ImageHeader>();

			GadgetContainerMessage< hoNDArray< std::complex<float> > > *cm =
					new GadgetContainerMessage< hoNDArray< std::complex<float> > >();

			*m->getObjectPtr() = j->image_headers_[frame];
			m->cont(cm);

			std::vector<size_t> img_dims(2);
			img_dims[0] = matrix_size_seq_[0];
			img_dims[1] = matrix_size_seq_[1];

			cm->getObjectPtr()->create(&img_dims);

			size_t data_length = prod(matrix_size_seq_);

			cudaMemcpy(cm->getObjectPtr()->get_data_ptr(),
					cgresult->get_data_ptr()+frame*data_length,
					data_length*sizeof(std::complex<float>),
					cudaMemcpyDeviceToHost);

			cudaError_t err = cudaGetLastError();
			if( err != cudaSuccess ){
				GADGET_DEBUG2("Unable to copy result from device to host: %s\n", cudaGetErrorString(err));
				m->release();
				return GADGET_FAIL;
			}

			m->getObjectPtr()->matrix_size[0] = matrix_size_seq_[0];
			m->getObjectPtr()->matrix_size[1] = matrix_size_seq_[1];
			m->getObjectPtr()->matrix_size[2] = 1;
			m->getObjectPtr()->channels       = 1;
			m->getObjectPtr()->image_index    = frame_counter_ + frame;

			if (this->next()->putq(m) < 0) {
				GADGET_DEBUG1("Failed to put result image on to queue\n");
				m->release();
				return GADGET_FAIL;
			}
		}
	} else{
		std::vector<size_t> img_dims(3);
		img_dims[0] = matrix_size_seq_[0];
		img_dims[1] = matrix_size_seq_[1];
		img_dims[2] = frames;

		auto cm =
				new GadgetContainerMessage< hoNDArray< std::complex<float> > >(img_dims);
		cgresult->to_host(reinterpret_cast<hoNDArray<float_complext>*>(cm->getObjectPtr()));
		auto m =
				new GadgetContainerMessage<ISMRMRD::ImageHeader>();

		*m->getObjectPtr() = j->image_headers_[0]; //Just use the first header
		m->cont(cm);

		m->getObjectPtr()->matrix_size[0] = matrix_size_seq_[0];
		m->getObjectPtr()->matrix_size[1] = matrix_size_seq_[1];
		m->getObjectPtr()->matrix_size[2] = frames;
		m->getObjectPtr()->channels       = 1;
		m->getObjectPtr()->image_index    = frame_counter_;

		if (this->next()->putq(m) < 0) {
			GADGET_DEBUG1("Failed to put result image on to queue\n");
			m->release();
			return GADGET_FAIL;
		}

	}

}

} /* namespace Gadgetron */
