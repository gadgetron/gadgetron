/**
	\brief CPU Gridding reconstruction gadget

	Handles reconstruction of 2D float data with 
	density compensation provided. Iterative reconstruction 
	can be easily integrated
*/

#pragma once 
#include "generic_recon_gadgets/GenericReconGadget.h"
#include "gadgetron_mri_noncartesian_export.h"
#include "hoNDArray.h"
#include "GriddingReconGadgetBase.h"

namespace Gadgetron{

	class EXPORTGADGETSMRINONCARTESIAN CPUGriddingReconGadget :public GriddingReconGadgetBase<hoNDArray> {

	public:

		GADGET_DECLARE(CPUGriddingReconGadget);

		CPUGriddingReconGadget();

		~CPUGriddingReconGadget();
	


	};
}
