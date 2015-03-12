
/** \file   mri_core_utility.h
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"

namespace Gadgetron
{
    // ------------------------------------------------------------------------
    // zero-padding and corresponding cut functions for MRI recon
    // ------------------------------------------------------------------------
    // compute the start and end index for zero padding
    // dstSize >= srcSize
    EXPORTMRICORE void zpadRange(size_t srcSize, size_t dstSize, size_t& start, size_t& end);

    // cut the center part of array
    template<typename T> EXPORTMRICORE void cutpad(const hoNDArray<T>& data, const std::vector<size_t>& cutSize, hoNDArray<T>& dataCut);

    // utility functions for 2D and 3D usage
    template<typename T> EXPORTMRICORE void cutpad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataCut);
    template<typename T> EXPORTMRICORE void cutpad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataCut);
}
