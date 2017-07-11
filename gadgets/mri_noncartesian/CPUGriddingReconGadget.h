#pragma once 
#include "GenericReconGadget.h"
#include "gadgetron_mri_noncartesian_export.h"
#include "hoNDArray.h"

namespace Gadgetron{
	class EXPORTGADGETSMRINONCARTESIAN CPUGriddingReconGadget : public GenericReconGadget {
	public:
		GADGET_DECLARE(CPUGriddingReconGadget);
		CPUGriddingReconGadget();
		~CPUGriddingReconGadget();
	
	protected:
		GADGET_PROPERTY(kernelWidthProperty, float, "Kernel width", 7);
		GADGET_PROPERTY(oversamplingFactorProperty, float, "Oversmapling factor", 1.5);
		GADGET_PROPERTY(iterateProperty, bool, "Iterate bool", false);

		float kernelWidth;
		std::vector<size_t> imageDims;
		std::vector<size_t> imageDimsOs;

		virtual int process_config(ACE_Message_Block *mb);
		virtual int process(GadgetContainerMessage <IsmrmrdReconData> *m1);

		boost::shared_ptr<hoNDArray<float_complext>> reconstruct(
			hoNDArray<float_complext> *data,
			hoNDArray<floatd2> *traj,
			hoNDArray<float> *dcw,
			size_t nCoils
		);

		boost::shared_ptr<hoNDArray<float_complext>> reconstructChannel(
			hoNDArray<float_complext> *data,
			hoNDArray<floatd2> *traj,
			hoNDArray<float> *dcw
		);

		std::tuple<boost::shared_ptr<hoNDArray<floatd2>>, boost::shared_ptr<hoNDArray<float>>> 
		separateDcwAndTraj(
			hoNDArray<float> *dcwTraj
		);
	};
}
