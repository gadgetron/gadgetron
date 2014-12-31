
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
    template<typename T> EXPORTMRICORE void zeropad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataPadded);

    // pad first three dimensions array around its center, other dimensions are kept unchanged
    template<typename T> EXPORTMRICORE void zeropad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded);

    // the dataPadded is not pre cleared to fill with zeros
    template<typename T> EXPORTMRICORE void zeropad3DNoPresetZeros(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded);

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
    * @brief permute E2 dimension of x : [RO E1 CHA SLC E2 ...] to r: [RO E1 E2 CHA SLC ...]
    */
    template<typename T> EXPORTMRICORE bool permuteE2To3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute E2 dimension of x : [RO E1 E2 CHA SLC ...] to r: [RO E1 CHA SLC E2 ...]
    */
    template<typename T> EXPORTMRICORE bool permuteE2To5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute RO dimension of x to the 3rd dimension
             x : [RO E1 E2 ...], r: [E1 E2 RO ...]
    */
    template<typename T> EXPORTMRICORE bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute RO dimension of x to the 4th dimension
             x : [RO E1 E2 CHA ...], r: [E1 E2 CHA RO ...]
    */
    template<typename T> EXPORTMRICORE bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute RO dimension of x back to the 1st dimension
             x : [E1 E2 CHA RO ...], r: [RO E1 E2 CHA ...]
    */
    template<typename T> EXPORTMRICORE bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute the 3rd dimension of x to the 1st dimension
             x : [RO E1 E2 CHA ...], r: [E2 RO E1 CHA ...]
    */
    template<typename T> EXPORTMRICORE bool permute3rdDimensionTo1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute RO dimension of x to the 5th dimension
             x : [RO E1 E2 srcCHA dstCHA ...], r: [E1 E2 srcCHA dstCHA RO ...]
    */
    template<typename T> EXPORTMRICORE bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief r = x add/multiply y for every part of y
    */
    template<typename T> EXPORTMRICORE bool multipleAdd(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);
    template<typename T> EXPORTMRICORE bool multipleMultiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);
}
