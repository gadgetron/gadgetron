/*
 * gpuBufferSensePrepGadget.h
 *
 *  Created on: Dec 10, 2014
 *      Author: dch
 */

#ifndef GPUBUFFERSENSEPREPGADGET_H_
#define GPUBUFFERSENSEPREPGADGET_H_

#include "Gadget.h"
#include "mri_core_data.h"
#include "cuNDArray.h"
#include "vector_td.h"
#include "complext.h"

namespace Gadgetron {

class gpuBufferSensePrepGadget: public Gadgetron::Gadget1<IsmrmrdReconData> {
public:
	gpuBufferSensePrepGadget();
	virtual ~gpuBufferSensePrepGadget();

	virtual int process_config(ACE_Message_Block*mb);

	virtual int process(GadgetContainerMessage<IsmrmrdReconData>* data);
protected:
	size_t profiles_per_frame_;
	float kernel_width;
	float oversampling_factor_;
	int ncoils_;
	std::vector<size_t> image_dims_;
	std::vector<size_t> image_dims_recon_;
	uint64d2 image_dims_recon_os_;
	boost::shared_ptr<cuNDArray<float_complext>> reconstruct_regularization(cuNDArray<float_complext>* data, cuNDArray<floatd2>* traj, cuNDArray<float>* dcw, size_t coils );
	static std::tuple<boost::shared_ptr<hoNDArray<floatd2 > >, boost::shared_ptr<hoNDArray<float >>> separate_traj_and_dcw(hoNDArray<float>*);
	ISMRMRD::ImageHeader create_image_header(ISMRMRD::AcquisitionHeader& header,const SamplingDescription& samp,size_t idx, size_t num_frames);

};

} /* namespace Gadgetron */
#endif /* GPUBUFFERSENSEPREPGADGET_H_ */
