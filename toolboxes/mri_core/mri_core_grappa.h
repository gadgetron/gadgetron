
/** \file   mri_core_grappa.h
    \brief  GRAPPA implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"

namespace Gadgetron {
    /// ---------------------------------------------------------------------
    /// 2D grappa
    /// ---------------------------------------------------------------------

    /// grappa 2d calibration function to compute convolution kernel
    /// acsSrc : calibration data for source channel [RO E1 srcCHA], full kspace
    /// acsDst : calibration data for destination channel [RO E1 dstCHA], full kspace
    /// accelFactor: acceleration factor
    /// thres: the threshold for regularization during kernel estimation
    /// kRO: kernel size along RO
    /// kNE1: kernel size along E1
    /// convKer: computed grappa convolution kernel
    template <typename T> EXPORTMRICORE void grappa2d_calib_convolution_kernel(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, ho4DArray<T>& convKer);
    /// dataMask : [RO E1] array, marking fully rectangular sampled region with 1
    template <typename T> EXPORTMRICORE void grappa2d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask, size_t accelFactor, double thres, size_t kRO, size_t kNE1, ho4DArray<T>& convKer);

    /// compute image domain kernel from 2d grappd convolution kernel
    /// RO, E1: the size of image domain kernel
    /// kIm: image domain kernel [RO E1 srcCHA dstCHA]
    template <typename T> EXPORTMRICORE void grappa2d_image_domain_kernel(const ho4DArray<T>& convKer, size_t RO, size_t E1, hoNDArray<T>& kIm);

    /// compute unmixing coefficient from image domain kernel and coil sensitivity
    /// kerIm: [RO E1 srcCHA dstCHA], image domain kernel
    /// coilMap: [RO E1 dstCHA] coil sensitivity map
    /// unmixCoeff: [RO E1 srcCHA] unmixing coefficient
    /// gFactor: [RO E1], gfactor
    template <typename T> EXPORTMRICORE void grappa2d_unmixing_coeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, size_t acceFactorE1, hoNDArray<T>& unmixCoeff, hoNDArray< typename realType<T>::Type >& gFactor);

    /// apply unmixing coefficient on undersampled kspace
    /// kspace: [RO E1 srcCHA ...]
    /// unmixCoeff : [RO E1 srcCHA]
    /// complexIm : [RO E1 ...] wrapped complex images
    template <typename T> EXPORTMRICORE void apply_unmix_coeff_kspace(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

    /// aliasedIm : [RO E1 srcCHA ...]
    template <typename T> EXPORTMRICORE void apply_unmix_coeff_aliased_image(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

    /// ------------------------
    /// grappa 2d low level functions
    /// ------------------------
    /// get the kernel pattern, given the acceleration factor and kernel size
    /// kE1: kernel pattern along E1
    /// oE1: output pattern along E1
    /// e.g. for R=4 and kNE1=4, the kernel pattern kE1 will be [-4 0 4 8] and the output pattern oE1 will be [0 1 2 3] if fitItself==true
    /// if fitItself==false, the output pattern oE1 will be [1 2 3]
    /// if the acsSrc and acsDst are generated in different ways, often fitItself needs to be true; e.g. acsSrc is in the origin acquired channels
    /// and acsDst is in eigen channel
    /// accelFactor: acceleration factor
    /// kNE1: kernel size along E1
    EXPORTMRICORE void grappa2d_kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, size_t accelFactor, size_t kNE1, bool fitItself);

    /// grappa calibration for 2D case
    /// kE1: the kernel pattern along E1
    /// oE1: the output kernel pattern along E1
    /// ker : kernel array [kRO kE1 srcCHA dstCHA oE1]
    template <typename T> EXPORTMRICORE void grappa2d_calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray<T>& ker);

    /// convert the grappa multiplication kernel computed from grappa2d_calib to convolution kernel
    /// convKer : [convRO convE1 srcCHA dstCHA]
    template <typename T> EXPORTMRICORE void grappa2d_convert_to_convolution_kernel(const ho5DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho4DArray<T>& convKer);

}
