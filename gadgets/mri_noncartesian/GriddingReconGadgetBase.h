#pragma once

#include "GenericReconGadget.h"
#include "gadgetron_mri_noncartesian_export.h"

namespace Gadgetron {

	template<template<class> class ARRAY>
	class GriddingReconGadgetBase : public ImageArraySendMixin<GriddingReconGadgetBase<ARRAY>>, public Gadget1<IsmrmrdReconData>
	{
	public:
//		GADGET_DECLARE(GriddingReconGadgetBase);
		GriddingReconGadgetBase();
		~GriddingReconGadgetBase();


		GADGET_PROPERTY(kernel_width,float,"Kernel width for NFFT", 5.5);
		GADGET_PROPERTY(gridding_oversampling_factor,float,"Oversampling used in NFFT", 1.5);
		GADGET_PROPERTY(iterate,bool,"Iterate instead of using weights", false);
		GADGET_PROPERTY(iteration_max,int,"Maximum number of iterations", 5);
		GADGET_PROPERTY(iteration_tol,float,"Iteration tolerance", 1e-5);
		GADGET_PROPERTY(replicas, int,"Number of pseudo replicas", 0);
		GADGET_PROPERTY(snr_frame, int,"Frame number for SNR measurement", 20);
		GADGET_PROPERTY(perform_timing, bool,"Perform timing", false);
		GADGET_PROPERTY(image_series,int,"Image Series",1);
		GADGET_PROPERTY(verbose, bool,"Verbose", false);
	protected:
		float kernel_width_;
		float oversampling_factor_;
		int ncoils_;

		unsigned int process_called_times=0;

		std::vector<size_t> image_dims_;
		uint64d2 image_dims_os_;

		virtual int process_config(ACE_Message_Block* mb) override;
		virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1) override;

		void pseudo_replica(const hoNDArray<std::complex<float>>& data,
                             ARRAY<floatd2>& traj,ARRAY<float>& dcw, const ARRAY<float_complext>& csm,
                             const IsmrmrdReconBit& recon_bit, size_t encoding, size_t ncoils);

		boost::shared_ptr<ARRAY<float_complext> > reconstruct(
			ARRAY<float_complext>* data,
			ARRAY<floatd2>* traj,
			ARRAY<float>* dcw,
			size_t ncoils );

		std::tuple<boost::shared_ptr<hoNDArray<floatd2 > >, boost::shared_ptr<hoNDArray<float >>> separate_traj_and_dcw(hoNDArray<float >* traj_dcw);

	};
}

#include "GriddingReconGadgetBase.hpp"
