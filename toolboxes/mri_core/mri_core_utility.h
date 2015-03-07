
/** \file   mri_core_utility.h
    \brief  Implementation useful utility functionalities for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"

namespace Gadgetron
{
    // ------------------------------------------------------------------------
    // zero-padding and corresponding cut functions for MRI recon
    // ------------------------------------------------------------------------
    // compute the start and end index for zero padding
    // dstSize >= srcSize
    EXPORTMRICORE void zpadRange(size_t srcSize, size_t dstSize, size_t& start, size_t& end);

    // pad the array around its center (N/2 is the center)
    // paddedSize: the array size after zero-padding
    // if paddedSize.size() < data.get_number_of_dimensions(), then first paddedSize.size() dimensions will be padded
    // if presetZeros==true, the dataPadded will be preset to fill with zeros before padding
    template<typename T> EXPORTMRICORE void zeropad(const hoNDArray<T>& data, const std::vector<size_t>& paddedSize, hoNDArray<T>& dataPadded, bool presetZeros = true);

    // pad the first two dimensions around its center, other dimensions are kept unchanged
    template<typename T> EXPORTMRICORE void zeropad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataPadded, bool presetZeros=true);

    // pad first three dimensions array around its center, other dimensions are kept unchanged
    template<typename T> EXPORTMRICORE void zeropad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded, bool presetZeros=true);

    // cut the center part of array
    template<typename T> EXPORTMRICORE void cutpad(const hoNDArray<T>& data, const std::vector<size_t>& cutSize, hoNDArray<T>& dataCut);

    // utility functions for 2D and 3D usage
    template<typename T> EXPORTMRICORE void cutpad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataCut);
    template<typename T> EXPORTMRICORE void cutpad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataCut);
}
