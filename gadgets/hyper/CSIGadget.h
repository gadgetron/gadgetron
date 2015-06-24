/*
 * CSIGadget.h
 *
 *  Created on: Nov 11, 2014
 *      Author: dch
 */

#ifndef CSIGADGET_H_
#define CSIGADGET_H_

#include "Gadget.h"
#include <ismrmrd/ismrmrd.h>
#include "GenericReconJob.h"
#include "CSIOperator.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuSbcCgSolver.h"
#include "gpuCSICoilEstimationGadget.h"
namespace Gadgetron {

class CSIGadget: public Gadgetron::Gadget1<cuSenseData>{
public:
	CSIGadget();
	virtual ~CSIGadget();

	virtual int process(GadgetContainerMessage<cuSenseData>* m1);
    virtual int process_config( ACE_Message_Block* mb );

protected:
    GADGET_PROPERTY(deviceno,int,"GPU device number", 0);
    GADGET_PROPERTY(number_of_cg_iterations,int,"Number of CG iterations", 10);
    GADGET_PROPERTY(number_of_sb_iterations,int,"Number of SB iterations",20);
    GADGET_PROPERTY(cg_limit,float,"CG limit",1e-5f);
    GADGET_PROPERTY(oversampling_factor,float,"Oversampling factor",1.5);
    GADGET_PROPERTY(kernel_width,float,"Kernel width",5.5);
    GADGET_PROPERTY(mu,float,"Mu",1.0);
    GADGET_PROPERTY(output_convergence,bool,"Output convergence",false);
    GADGET_PROPERTY(frequencies,std::vector<float>,"Frequencies",{});
    GADGET_PROPERTY(frequency_min,float,"Minimum frequency",0);
    GADGET_PROPERTY(frequency_max,float,"Maximum frequency",1300);
    GADGET_PROPERTY(frequency_step,float,"Frequency step",1300);
    GADGET_PROPERTY(compressed_sensing,bool,"Use compressed sensing",true);



    int device_number_;
    unsigned int number_of_cg_iterations_;
    unsigned int number_of_sb_iterations_;
    float cg_limit_;
    float oversampling_factor_;
    float kernel_width_;
    float mu_;
    float lambda_;
    bool output_convergence_;
    bool use_compressed_sensing_;

    std::vector<size_t> img_dims_;

    uint64d2 matrix_size_;
    uint64d2 matrix_size_os_;

    boost::shared_ptr<CSIOperator<float> > E_;
    boost::shared_ptr<cuNonCartesianSenseOperator<float,2> > S_;
    cuSbcCgSolver<float_complext> solver_;
    bool is_configured_;
    unsigned int channels_;
};

} /* namespace Gadgetron */
#endif /* CSIGADGET_H_ */
