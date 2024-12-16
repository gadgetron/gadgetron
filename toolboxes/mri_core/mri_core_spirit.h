
/** \file   mri_core_spirit.h
    \brief  SPIRIT implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

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
    /// if minusI==true, compute image domain G-I kernel
    template <typename T> void spirit2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, 
                                                                            size_t startRO, size_t endRO, size_t startE1, size_t endE1, 
                                                                            hoNDArray<T>& convKer, bool minusI = false);
    /// entire data in acsSrc and acsDst is used
    template <typename T> void spirit2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, 
                                                                            hoNDArray<T>& convKer, bool minusI = false);

    /// dataMask : [RO E1] array, marking fully rectangular sampled region with 1
    template <typename T> void spirit2d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, 
                                                                            hoNDArray<unsigned short>& dataMask, double thres, 
                                                                            size_t kRO, size_t kE1, size_t oRO, size_t oE1, 
                                                                            hoNDArray<T>& convKer, bool minusI = false);

    /// compute image domain kernel from 2d grappd convolution kernel
    /// RO, E1: the size of image domain kernel
    /// kIm: image domain kernel [RO E1 srcCHA dstCHA]
    template <typename T> void spirit2d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, hoNDArray<T>& kIm);

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
    /// if minusI==true, compute image domain G-I kernel
    template <typename T> void spirit3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            double thres, double overDetermineRatio,
                                                                            size_t kRO, size_t kE1, size_t kE2, 
                                                                            size_t oRO, size_t oE1, size_t oE2,
                                                                            size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, 
                                                                            hoNDArray<T>& convKer, bool minusI = false);

    /// entire data in acsSrc and acsDst is used
    template <typename T> void spirit3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            double thres, double overDetermineRatio,
                                                                            size_t kRO, size_t kE1, size_t kE2, 
                                                                            size_t oRO, size_t oE1, size_t oE2, 
                                                                            hoNDArray<T>& convKer, bool minusI = false);

    /// dataMask : [RO E1 E2] array, marking fully rectangular sampled region with 1
    template <typename T> void spirit3d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask, 
                                                                            double thres, double overDetermineRatio, 
                                                                            size_t kRO, size_t kE1, size_t kE2, 
                                                                            size_t oRO, size_t oE1, size_t oE2,
                                                                            hoNDArray<T>& convKer, bool minusI = false);

    /// compute image domain kernel from 3d spirit convolution kernel
    /// RO, E1, E2: the size of image domain kernel
    /// kIm: image domain kernel [RO E1 E2 srcCHA dstCHA]
    /// if preset_kIm_with_zeros==false, caller should make sure the kIm is cleared with zeros
    template <typename T> void spirit3d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, size_t E2, hoNDArray<T>& kIm, bool preset_kIm_with_zeros=true);

    /// compute image domain kernel from 3d spirit convolution kernel
    // E1 and E2 stays in the kspace domain
    // kImRO: kspace-image hybrid kernel [convE1 convE2 srcCHA dstCHA RO]
    template <typename T> void spirit3d_kspace_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, hoNDArray<T>& kImRO);

    /// convert kspace-image hybrid kernel to image domain
    // kImRO : kspace-image hybrid kernel [convE1 convE2 srcCHA dstCHA RO]
    // kIm: image domain kernel [E1 E2 srcCHA dstCHA RO]
    template <typename T> void spirit3d_image_domain_kernel(const hoNDArray<T>& kImRO, size_t E1, size_t E2, hoNDArray<T>& kIm);

    /// perform spirit calibration, get the mulitplication kernel
    /// ker: [kRO kE1 kE2 srcCHA dstCHA oRO oE1 oE2]
    template <typename T> void spirit3d_calib(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                                                        double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2,
                                                        size_t oRO, size_t oE1, size_t oE2,
                                                        size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
                                                        hoNDArray<T>& ker);

    /// convert the spirit multiplication kernel computed from spirit2d_calib to convolution kernel
    /// kerL [kRO, kE1, srcCHA, dstCHA, oRO, oE1]
    /// convKer : [convRO convE1 srcCHA dstCHA]
    template <typename T> void spirit3d_convert_to_convolution_kernel(const hoNDArray<T>& ker, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray<T>& convKer, bool minusI = false);

    /// ---------------------------------------------------------------------
    /// general functionalities
    /// ---------------------------------------------------------------------

    /// compute the image domain adjoint kernel
    /// kIm: [... srcCHA dstCHA]
    /// adjkIm: [... dstCHA srcCHA]
    template <typename T> void spirit_image_domain_adjoint_kernel(const hoNDArray<T>& kIm, hoNDArray<T>& adjkIm);

    /// compute the (G-I)'*(G-I)
    template <typename T> void spirit_adjoint_forward_kernel(const hoNDArray<T>& kImS2D, const hoNDArray<T>& kImD2S, hoNDArray<T>& kIm);
}
