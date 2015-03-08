
/** \file   mri_core_coil_map_estimation.h
    \brief  Implementation MRI coil sensitivity map estimation functions.

    ISMRMRD_SOUHEIL coil map estimation is based on:

        Inati SJ, Hansen MS, Kellman P.
        A solution to the phase problem in adaptive coil combination.
        In: ISMRM proceeding; april; salt lake city, utah, USA. ; 2013. 2672.

        Kellman P, McVeigh ER.
        Image reconstruction in SNR units: A general method for SNR measurement.
        Magnetic Resonance in Medicine 2005;54(6):1439-1447.

    ISMRMRD_SOUHEIL_ITER coil map estimation is based on:

        Inati SJ, Hansen MS, Kellman P. Unpublished algorithm.

    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"

namespace Gadgetron
{
    // the Souheil method
    // data: [RO E1 CHA], only 3D array
    // these functions are using 2D data correlation matrix
    // ks: the kernel size for local covariance estimation
    // power: number of iterations to apply power method
    template<typename T> EXPORTMRICORE void coil_map_2d_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t power);

    // data: [RO E1 E2 CHA], this functions uses true 3D data correlation matrix
    template<typename T> EXPORTMRICORE void coil_map_3d_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t power);

    // the Souheil iteration method
    // data: [RO E1 CHA], only 3D array
    // ks: the kernel size for local covariance estimation
    // iterNum: number of iterations to refine coil map
    // thres: threshold to stop the iterations
    template<typename T> EXPORTMRICORE void coil_map_2d_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t iterNum, typename realType<T>::Type thres);

    // data: [RO E1 E2 CHA], true 3D coil map estimation
    template<typename T> EXPORTMRICORE void coil_map_3d_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t kz, size_t iterNum, typename realType<T>::Type thres);
}
