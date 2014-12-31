
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

    // ------------------------------------------------------------------------
    // array manipulation
    // ------------------------------------------------------------------------
    /**
    * @brief sum over last dimension of an array
             e.g. for a 4D array, sum over the 4th dimension and get a 3D array
    */
    template<typename T> EXPORTMRICORE bool sumOverLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r); 

    /**
    * @brief sum over the 1st dimension of an array
             e.g. for a 2D array, sum over the 1st dimension and get an array of [1 E1]
    */
    template<typename T> EXPORTMRICORE bool sumOver1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief sum over the 2nd dimension of an array
             e.g. for a 3D array, sum over the 2nd dimension and get an array of [RO 1 CHA]
    */
    template<typename T> EXPORTMRICORE bool sumOver2ndDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief sum over the 3rd dimension of an array
             e.g. for a 4D array, sum over the 3rd dimension and get an array of [RO E1 1 N]
    */
    template<typename T> EXPORTMRICORE bool sumOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief sum over the 4th dimension of an array
             e.g. for a 5D array [RO E1 CHA N S], sum over the 4th dimension and get an array of [RO E1 CHA 1 S]
    */
    template<typename T> EXPORTMRICORE bool sumOver4thDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief sum over the 5th dimension of an array
             e.g. for a 6D array, sum over the 5th dimension and get an array [RO E1 CHA N 1 P]
    */
    template<typename T> EXPORTMRICORE bool sumOver5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief r = x multiply y for every part of y
    */
    template<typename T> EXPORTMRICORE bool multipleMultiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);
}
