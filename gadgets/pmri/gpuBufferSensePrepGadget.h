/*
 * gpuBufferSensePrepGadget.h
 *
 *  Created on: Dec 10, 2014
 *      Author: dch
 */

#ifndef GPUBUFFERSENSEPREPGADGET_H_
#define GPUBUFFERSENSEPREPGADGET_H_

#include "Gadget.h"
#include "cuNDArray.h"
#include "vector_td.h"
#include "complext.h"

namespace Gadgetron {

class gpuBufferSensePrepGadget: public Gadgetron::Gadget1<mrd::ReconData> {
public:
	gpuBufferSensePrepGadget();
	virtual ~gpuBufferSensePrepGadget();

	virtual int process_config(const mrd::Header& header);

	virtual int process(GadgetContainerMessage<mrd::ReconData>* data);
protected:
	GADGET_PROPERTY(profiles_per_frame,int,"Number of profiles per frame", 0);
	GADGET_PROPERTY(kernel_width,float,"Kernel width for NFFT", 5.5);
	GADGET_PROPERTY(buffer_convolution_oversampling_factor,float,"Oversampling used in buffer NFFT", 1.5);
	GADGET_PROPERTY(reconstruction_os_factor,float,"Oversampling for recon NFFT", 1.5);

	size_t profiles_per_frame_;
	float kernel_width_;
	float oversampling_factor_;
	int ncoils_;
	std::vector<size_t> image_dims_;
	std::vector<size_t> image_dims_recon_;
	uint64d2 image_dims_recon_os_;
	boost::shared_ptr<cuNDArray<float_complext>> reconstruct_regularization(cuNDArray<float_complext>* data, cuNDArray<floatd2>* traj, cuNDArray<float>* dcw, size_t coils );
	static std::tuple<boost::shared_ptr<hoNDArray<floatd2 > >, boost::shared_ptr<hoNDArray<float >>> separate_traj_and_dcw(hoNDArray<float>*);
	mrd::ImageHeader create_image_header(mrd::AcquisitionHeader& header,const mrd::SamplingDescription& samp,size_t idx, size_t num_frames);

};

} /* namespace Gadgetron */
#endif /* GPUBUFFERSENSEPREPGADGET_H_ */
