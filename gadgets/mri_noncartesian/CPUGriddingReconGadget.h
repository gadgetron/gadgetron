/**
	\brief CPU Gridding reconstruction gadget

	Handles reconstruction of 2D float data with 
	density compensation provided. Iterative reconstruction 
	can be easily integreated
*/

#pragma once 
#include "GenericReconGadget.h"
#include "gadgetron_mri_noncartesian_export.h"
#include "hoNDArray.h"

namespace Gadgetron{

	class EXPORTGADGETSMRINONCARTESIAN CPUGriddingReconGadget : public GenericReconGadget {

	/**
		CPUGriddingReconGadget class declaration 
		----------------------------------------

		interface similar to GriddingReconGadget
		
		Note: can probably be combined and templated
	*/
	
	public:

		GADGET_DECLARE(CPUGriddingReconGadget);

		CPUGriddingReconGadget();

		~CPUGriddingReconGadget();
	
	protected:
		
		/**
			Gadget properties declaration
		*/

		GADGET_PROPERTY(kernelWidthProperty, float, "Kernel width", 5.5);
		GADGET_PROPERTY(oversamplingFactorProperty, float, "Oversmapling factor", 1.5);
		GADGET_PROPERTY(iterateProperty, bool, "Iterate bool", false);

		/**
			Storage for the properties above
		*/

		float kernelWidth;
		float oversamplingFactor;

		/**
			Image dimensions
		*/

		std::vector<size_t> imageDims;
		std::vector<size_t> imageDimsOs;

		virtual int process_config(ACE_Message_Block *mb);
		virtual int process(GadgetContainerMessage <IsmrmrdReconData> *m1);

		/**
			Reconstruct multi channel data

			/param data: multi-channel k-space data
			/param traj: trajectories
			/param dcw: density compensation
			/param nCoils: number of channels
		*/

		boost::shared_ptr<hoNDArray<float_complext>> reconstruct(
			hoNDArray<float_complext> *data,
			hoNDArray<floatd2> *traj,
			hoNDArray<float> *dcw,
			size_t nCoils
		);

		/**
			Reconstruct single channel

			/param data: single channel data
			/param traj: trajectory
			/param dcw: density compensation
		*/

		boost::shared_ptr<hoNDArray<float_complext>> reconstructChannel(
			hoNDArray<float_complext> *data,
			hoNDArray<floatd2> *traj,
			hoNDArray<float> *dcw
		);

		/**
			Helper method for seperating the trajectory
			and the density compesnation
		*/

		std::tuple<boost::shared_ptr<hoNDArray<floatd2>>, boost::shared_ptr<hoNDArray<float>>> 
		separateDcwAndTraj(
			hoNDArray<float> *dcwTraj
		);
	};
}
