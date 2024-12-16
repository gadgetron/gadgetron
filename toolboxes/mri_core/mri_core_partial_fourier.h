
/** \file   mri_core_partial_fourier.h
    \brief  Implementation partial fourier handling functionalities for 2D and 3D MRI
    \author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron
{
    /// perform the POCS partial fourier handling
    /// kspace: [RO E1 E2 CHA N S SLC]
    /// startRO/E1/E2: mark the start of sampling region along RO/E1/E2
    /// endRO/E1/E2: mark the end of sampling region along RO/E1/E2
    /// transit_band_RO/E1/E2: a transition band can be created between the sampled kspace and filled kspace region; if set to be 0, no trasit band is applied
    /// iter: number of maximal iterations for POCS
    /// thres: threshold to stop the iteration
    /// res: [RO E1 E2 CHA N S SLC], result of POCS
    template <typename T> void partial_fourier_POCS(const hoNDArray<T>& kspace,
        size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
        size_t transit_band_RO, size_t transit_band_E1, size_t transit_band_E2,
        size_t iter, double thres, hoNDArray<T>& res);

    /// perform the partial fourier handling filter
    /// filter_pf_width_RO/E1/E2: tapered width ratio [0 1] for the partial foureir filter
    /// filter_pf_density_comp: whether to perform density compensation
    /// filter_pf_RO/E1/E2: the partial fourier filter; if not set, filter will be computed insided the function
    template <typename T> void partial_fourier_filter(const hoNDArray<T>& kspace,
        size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
        double filter_pf_width_RO, double filter_pf_width_E1, double filter_pf_width_E2, bool filter_pf_density_comp, 
        hoNDArray<T>& filter_pf_RO, hoNDArray<T>& filter_pf_E1, hoNDArray<T>& filter_pf_E2, hoNDArray<T>& res);
}
