
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
    // zero-padding resize functions
    // ------------------------------------------------------------------------
    // compute the start and end index for zero padding
    // dstSize >= srcSize
    EXPORTMRICORE void zpadRange(size_t srcSize, size_t dstSize, size_t& start, size_t& end);

    // pad the first two dimensions around its center, other dimensions are kept unchanged
    // if presetZeros=true, the dataPadded will be preset to zero before padding
    template<typename T> EXPORTMRICORE void zeropad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataPadded, bool presetZeros=true);

    // pad first three dimensions array around its center, other dimensions are kept unchanged
    template<typename T> EXPORTMRICORE void zeropad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded, bool presetZeros=true);

    // cut the center part
    template<typename T> EXPORTMRICORE void cutpad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataCut);
    template<typename T> EXPORTMRICORE void cutpad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataCut);
}
