/*
 * gpuCSICoilEstimationGadget.h
 *
 *  Created on: Nov 12, 2014
 *      Author: dch
 */

#ifndef gpuCSICoilESTIMATIONGADGET_H_
#define gpuCSICoilESTIMATIONGADGET_H_

#include "Gadget.h"
#include <tuple>
#include "mri_core_data.h"
#include "cuNDArray.h"

namespace Gadgetron {

class gpuCSICoilEstimationGadget: public Gadgetron::Gadget1<IsmrmrdAcquisitionBucket> {
public:
	gpuCSICoilEstimationGadget();
	virtual ~gpuCSICoilEstimationGadget();
	virtual int process_config(ACE_Message_Block* mb);

    virtual int process(GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1);

protected:

    static  std::tuple<hoNDArray<std::complex<float>>*, hoNDArray<float>*> combine_data(std::vector<IsmrmrdAcquisitionData>& aquisitions);

    boost::shared_ptr<cuNDArray<float_complext>> calculate_CSM(cuNDArray<float_complext>* data,cuNDArray<floatd2>* traj, cuNDArray<float>* dcw );
/**
 * Separates trajectory and dcw
 * @param array containing trajectory and dcw, so that the first dimension has size 3
 * @return tuple containing trajectory and dcw
 */
    static std::tuple<boost::shared_ptr<hoNDArray<floatd2 > >, boost::shared_ptr<hoNDArray<float >>> separate_traj_and_dcw(hoNDArray<float>*);
    std::vector<size_t> img_size;
    size_t coils;
    float kernel_width;

    unsigned int skip_lines;

};


struct cuSenseData {
	boost::shared_ptr<cuNDArray<float_complext>> data;
	boost::shared_ptr<cuNDArray<floatd2>> traj;
	boost::shared_ptr<cuNDArray<float_complext>> csm;
	boost::shared_ptr<cuNDArray<float>> dcw;
	boost::shared_ptr<cuNDArray<float_complext>> freq_calibration;

};

} /* namespace Gadgetron */
#endif /* gpuCSICoilESTIMATIONGADGET_H_ */
