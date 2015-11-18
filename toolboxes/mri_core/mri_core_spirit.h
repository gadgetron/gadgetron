
/** \file   mri_core_spirit.h
    \brief  SPIRIT implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"

namespace Gadgetron {
    /// ---------------------------------------------------------------------
    /// 2D spirit
    /// ---------------------------------------------------------------------

    /// spirit 2d calibration function to compute convolution kernel
    /// acsSrc : calibration data for source channel [RO E1 srcCHA], full kspace
    /// acsDst : calibration data for destination channel [RO E1 dstCHA], full kspace
    /// startRO, endRO, startE1, endE1: define the data region [startRO endRO], [startE1 endE1] which is used for calibration
    /// thres: the threshold for regularization during kernel estimation
    /// kRO: kernel size along RO
    /// kE1: kernel size along E1
    /// oRO: kernel output size along RO
    /// oE1: kernel output size along E1
    /// convKer: computed spirit convolution kernel
    template <typename T> EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, 
                                                                            size_t startRO, size_t endRO, size_t startE1, size_t endE1, 
                                                                            hoNDArray<T>& convKer);
    /// entire data in acsSrc and acsDst is used
    template <typename T> EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, 
                                                                            hoNDArray<T>& convKer);

    /// dataMask : [RO E1] array, marking fully rectangular sampled region with 1
    template <typename T> EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, 
                                                                            hoNDArray<unsigned short>& dataMask, double thres, 
                                                                            size_t kRO, size_t kE1, size_t oRO, size_t oE1, 
                                                                            hoNDArray<T>& convKer);

    /// compute image domain kernel from 2d grappd convolution kernel
    /// RO, E1: the size of image domain kernel
    /// kIm: image domain kernel [RO E1 srcCHA dstCHA]
    /// if minusI==true, compute image domain G-I kernel
    template <typename T> EXPORTMRICORE void spirit2d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, hoNDArray<T>& kIm, bool minusI = false);

    /// ---------------------------------------------------------------------
    /// 3D spirit
    /// ---------------------------------------------------------------------

    /// spirit 3d calibration function to compute convolution kernel
    /// acsSrc : calibration data for source channel [RO E1 E2 srcCHA], full kspace
    /// acsDst : calibration data for destination channel [RO E1 E2 dstCHA], full kspace
    /// startRO, endRO, startE1, endE1: define the data region [startRO endRO], [startE1 endE1], [startE2 endE2] which is used for calibration
    /// thres: the threshold for regularization during kernel estimation
    /// overDetermineRatio: the calibration matrix over determination ratio
    /// kRO: kernel size along RO
    /// kE1: kernel size along E1
    /// kE2: kernel size along E2
    /// oRO: kernel output size along RO
    /// oE1: kernel output size along E1
    /// oE2: kernel output size along E2
    /// convKer: computed spirit convolution kernel [convKRO convKE1 convKE2 srcCHA dstCHA]
    template <typename T> EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            double thres, double overDetermineRatio,
                                                                            size_t kRO, size_t kE1, size_t kE2, 
                                                                            size_t oRO, size_t oE1, size_t oE2,
                                                                            size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, 
                                                                            hoNDArray<T>& convKer);

    /// entire data in acsSrc and acsDst is used
    template <typename T> EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            double thres, double overDetermineRatio,
                                                                            size_t kRO, size_t kE1, size_t kE2, 
                                                                            size_t oRO, size_t oE1, size_t oE2, 
                                                                            hoNDArray<T>& convKer);

    /// dataMask : [RO E1 E2] array, marking fully rectangular sampled region with 1
    template <typename T> EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask, 
                                                                            double thres, double overDetermineRatio, 
                                                                            size_t kRO, size_t kE1, size_t kE2, 
                                                                            size_t oRO, size_t oE1, size_t oE2,
                                                                            hoNDArray<T>& convKer);

    /// compute image domain kernel from 3d grappd convolution kernel
    /// RO, E1, E2: the size of image domain kernel
    /// kIm: image domain kernel [RO E1 E2 srcCHA dstCHA]
    /// if minusI==true, compute image domain G-I kernel
    template <typename T> EXPORTMRICORE void spirit3d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, size_t E2, hoNDArray<T>& kIm, bool minusI = false);
}
