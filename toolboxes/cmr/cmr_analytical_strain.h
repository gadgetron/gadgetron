/** \file   cmr_radial_thickening.h
	\brief  Implement functionalities to handle cardiac time stamps
	\author Angela Gao
*/

#pragma once

#include "cmr_export.h"

#include "GadgetronTimer.h"

#include "hoNDArray.h"

namespace Gadgetron {

	/// compute radial strain map over all phases with reference to ref_phase
	/// endo_mask, epi_mask: [RO, E1, N], 2D+T array of masks over all phases
	/// ref phase: size_t that is the peak phase
	template <typename T> EXPORTCMR void compute_analytical_strain(const hoNDArray<double>& dx, const hoNDArray<double>& dy, const hoNDArray<T>& mask, hoNDArray<T>& rad_strain, hoNDArray<T>& circ_strain);
}
