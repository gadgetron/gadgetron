/*
 * gpuSenseGadget.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: dch
 */

#include "gpuSenseGadget.h"
#include "cuNDArray.h"
#include "vector_td_utilities.h"
#include "hoNDArray_math.h"
#include "cuNDArray_math.h"
#include "cudaDeviceManager.h"
namespace Gadgetron {

gpuSenseGadget::gpuSenseGadget() 
	: frame_counter_(0) 
{
	matrix_size_ = uint64d2(0u,0u);
	matrix_size_os_ = uint64d2(0u,0u);
	matrix_size_seq_ = uint64d2(0u,0u);
}

gpuSenseGadget::~gpuSenseGadget() {
	// TODO Auto-generated destructor stub
}

int gpuSenseGadget::process_config(ACE_Message_Block* mb) {
  device_number_ = deviceno.value();

  int number_of_devices = cudaDeviceManager::Instance()->getTotalNumberOfDevice();

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

  set_number_ = setno.value();
  slice_number_ = sliceno.value();
  oversampling_factor_ = oversampling_factor.value();
  kernel_width_ = kernel_width.value();
  output_convergence_ = output_convergence.value();
  output_timing_ = output_timing.value();
  rotations_to_discard_ = rotations_to_discard.value();
  
  if( (rotations_to_discard_%2) == 1 ){
    GDEBUG("#rotations to discard must be even.\n");
    return GADGET_FAIL;
  }
  save_individual_frames_ = save_individual_frames.value();
  return GADGET_OK;
}

int gpuSenseGadget::put_frames_on_que(int frames,int rotations, GenericReconJob* j, cuNDArray<float_complext>* cgresult,int channels) {

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

			std::vector<size_t> img_dims {cgresult->get_size(0),cgresult->get_size(1)};

			cm->getObjectPtr()->create(img_dims);

			size_t data_length = cm->getObjectPtr()->get_number_of_bytes();

			cudaMemcpy(cm->getObjectPtr()->get_data_ptr(),
					cgresult->get_data_ptr()+frame*cm->getObjectPtr()->get_number_of_elements(),
					data_length,
					cudaMemcpyDeviceToHost);

			cudaError_t err = cudaGetLastError();
			if( err != cudaSuccess ){
				GDEBUG("Unable to copy result from device to host: %s\n", cudaGetErrorString(err));
				m->release();
				return GADGET_FAIL;
			}

			m->getObjectPtr()->matrix_size[0] = img_dims[0];
			m->getObjectPtr()->matrix_size[1] = img_dims[1];
			m->getObjectPtr()->matrix_size[2] = 1;
			m->getObjectPtr()->channels       = 1;
			m->getObjectPtr()->image_index    = frame_counter_ + frame;

			if (this->next()->putq(m) < 0) {
				GDEBUG("Failed to put result image on to queue\n");
				m->release();
				return GADGET_FAIL;
			}
		}
	} else{
		std::vector<size_t> img_dims { cgresult->get_size(0),cgresult->get_size(1),(size_t)frames};

		auto cm =
				new GadgetContainerMessage< hoNDArray< std::complex<float> > >(img_dims);
		cgresult->to_host(reinterpret_cast<hoNDArray<float_complext>*>(cm->getObjectPtr()));
		auto m =
				new GadgetContainerMessage<ISMRMRD::ImageHeader>();

		*m->getObjectPtr() = j->image_headers_[0]; //Just use the first header
		m->cont(cm);

		m->getObjectPtr()->matrix_size[0] = img_dims[0];
		m->getObjectPtr()->matrix_size[1] = img_dims[1];
		m->getObjectPtr()->matrix_size[2] = img_dims[2];
		m->getObjectPtr()->channels       = channels;
		m->getObjectPtr()->image_index    = frame_counter_;

		if (this->next()->putq(m) < 0) {
			GDEBUG("Failed to put result image on to queue\n");
			m->release();
			return GADGET_FAIL;
		}

	}
	return GADGET_OK;

}

} /* namespace Gadgetron */
