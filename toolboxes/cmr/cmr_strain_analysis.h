/** \file   cmr_stain_analysis.h
    \brief  Implement functionalities to handle cardiac time stamps
    \author Hui Xue
*/

#pragma once

#include "cmr_export.h"

#include "GadgetronTimer.h"

#include "hoNDArray.h"

namespace Gadgetron { 

    /// compute radial and circ strain map from deformation fields
	/// dx, dy: [RO, E1, N], 2D+T array of deformation fields
    template <typename T> EXPORTCMR void compute_strain(const hoNDArray<double>& dx, const hoNDArray<double>& dy, const hoNDArray<T>& mask, const bool compare_mask,  hoNDArray<T>& radial, hoNDArray<T>& circ, hoNDArray<T>& thetas);
}
