
/** \file   mri_core_coil_map_estimation.h
    \brief  Implementation MRI coil sensitivity map estimation functions.

    Inati coil map estimation is based on:

        Inati SJ, Hansen MS, Kellman P.
        A solution to the phase problem in adaptive coil combination.
        In: ISMRM proceeding; April; Salt Lake City, Utah, USA; 2013. 2672.

        Kellman P, McVeigh ER.
        Image reconstruction in SNR units: A general method for SNR measurement.
        Magnetic Resonance in Medicine 2005;54(6):1439-1447.

    Inati_Iter coil map estimation is based on:

        Inati SJ, Hansen MS, Kellman P. 
        A Fast Optimal Method for Coil Sensitivity Estimation and Adaptive Coil Combination for Complex Images.
        In: ISMRM proceeding; May; Milan, Italy; 2014. 4407.

    \author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron
{
    // the Souheil method
    // data: [RO E1 CHA], only 3D array
    // these functions are using 2D data correlation matrix
    // ks: the kernel size for local covariance estimation
    // power: number of iterations to apply power method
    template<typename T>  void coil_map_2d_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks = 7, size_t power = 3);

    // data: [RO E1 E2 CHA], this functions uses true 3D data correlation matrix
    template<typename T>  void coil_map_3d_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks = 7, size_t kz = 5, size_t power = 3);

    // data: [RO E1 E2 CHA N S SLC ...], if E2==1, the 2D coil map estimation is assumed
    template<typename T>  void coil_map_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks = 7, size_t kz = 5, size_t power = 3);
    template<typename T>  hoNDArray<T> coil_map_Inati(const hoNDArray<T>& data, size_t ks = 7, size_t kz = 5, size_t power = 3);

    // the Souheil iteration method
    // data: [RO E1 CHA], only 3D array
    // ks: the kernel size for local covariance estimation
    // iterNum: number of iterations to refine coil map
    // thres: threshold to stop the iterations
    template<typename T>  void coil_map_2d_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks=7, size_t iterNum=5, typename realType<T>::Type thres=1e-3);

    // data: [RO E1 E2 CHA], true 3D coil map estimation
    template<typename T>  void coil_map_3d_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks=7, size_t kz=5, size_t iterNum=5, typename realType<T>::Type thres=0.001);

    // data: [RO E1 E2 CHA N S SLC ...], if E2==1, the 2D coil map estimation is assumed
    template<typename T>  void coil_map_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks=7, size_t kz=5, size_t iterNum=5, typename realType<T>::Type thres=0.001);

    template<class REAL,unsigned int D>
     hoNDArray<complext<REAL>>  estimate_b1_map(const hoNDArray<complext<REAL>>& data);

    // coil combination
    // the cha_dim = 2 for 2D case, e.g.
    // data: in image domain, at least 3D [RO E1 CHA ...]
    // coilMap: [RO E1 CHA ... ]
    // combined: [RO E1 ...]
    /// for the 3d, cha_dim=3, and 
    // data: in image domain, [RO E1 E2 CHA ...]
    // coilMap: [RO E1 E2 CHA ... ]
    // combined: [RO E1 E2 ...]
    template<typename T>  void coil_combine(const hoNDArray<T>& data, const hoNDArray<T>& coilMap, size_t cha_dim, hoNDArray<T>& combined);
    template<typename T>  hoNDArray<T> coil_combine(const hoNDArray<T>& data, const hoNDArray<T>& coilMap, size_t cha_dim);
}
