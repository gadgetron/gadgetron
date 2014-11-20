/*
 * gpuSenseGadget.h
 *
 *  Created on: Nov 17, 2014
 *      Author: dch
 */

#ifndef GPUSENSEGADGET_H_
#define GPUSENSEGADGET_H_

#include "Gadget.h"
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include "vector_td.h"
#include "GenericReconJob.h"
#include "cuNDArray.h"
namespace Gadgetron {

class gpuSenseGadget: public Gadget2<ISMRMRD::ImageHeader, GenericReconJob>{
public:
	gpuSenseGadget();
	virtual ~gpuSenseGadget();
	virtual int process_config(ACE_Message_Block* mb);

protected:

	virtual int put_frames_on_que(int frames,int rotations, GenericReconJob* j, cuNDArray<float_complext>* cgresult);
	 int channels_;
    int device_number_;
    int set_number_;
    int slice_number_;

    uint64d2 matrix_size_;
    uint64d2 matrix_size_os_;
    uint64d2 matrix_size_seq_;
    double oversampling_factor_;
    double kernel_width_;
    unsigned int rotations_to_discard_;

    bool output_convergence_;
    bool output_timing_;
    bool save_individual_frames_;

    int frame_counter_;
};

} /* namespace Gadgetron */
#endif /* GPUSENSEGADGET_H_ */
