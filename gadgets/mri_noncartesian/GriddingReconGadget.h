#pragma once

#include "GenericReconGadget.h"
#include "gadgetron_mri_noncartesian_export.h"
#include "cuNDArray.h"

namespace Gadgetron {

	class EXPORTGADGETSMRINONCARTESIAN GriddingReconGadget : public GenericReconGadget
	{
	public:
		GADGET_DECLARE(GriddingReconGadget);

		typedef GenericReconGadget BaseClass;

		GriddingReconGadget();
		~GriddingReconGadget();

	protected:
		GADGET_PROPERTY(kernel_width,float,"Kernel width for NFFT", 5.5);
		GADGET_PROPERTY(gridding_oversampling_factor,float,"Oversampling used in NFFT", 1.5);
		GADGET_PROPERTY(iterate,bool,"Iterate instead of using weights", false);
		GADGET_PROPERTY(iteration_max,int,"Maximum number of iterations", 5);
		GADGET_PROPERTY(iteration_tol,float,"Iteration tolerance", 1e-5);
		GADGET_PROPERTY(replicas, int,"Number of pseudo replicas", 0);
		GADGET_PROPERTY(snr_frame, int,"Frame number for SNR measurement", 20);
	
		float kernel_width_;
		float oversampling_factor_;
		int ncoils_;

		std::vector<size_t> image_dims_;
		uint64d2 image_dims_os_;

		virtual int process_config(ACE_Message_Block* mb) override;
		virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1) override;
		virtual void compute_image_header(IsmrmrdReconBit& recon_bit, IsmrmrdImageArray& res, size_t e) override;

		boost::shared_ptr<cuNDArray<float_complext> > reconstruct(
			cuNDArray<float_complext>* data,
			cuNDArray<floatd2>* traj,
			cuNDArray<float>* dcw,
			size_t ncoils );

		std::tuple<boost::shared_ptr<hoNDArray<floatd2 > >, boost::shared_ptr<hoNDArray<float >>> separate_traj_and_dcw(hoNDArray<float >* traj_dcw);

	};
}
