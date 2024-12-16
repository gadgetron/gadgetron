
/** \file   mri_core_grappa.h
    \brief  GRAPPA implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron {
    /// ---------------------------------------------------------------------
    /// 2D grappa
    /// ---------------------------------------------------------------------

    /// grappa 2d calibration function to compute convolution kernel
    /// acsSrc : calibration data for source channel [RO E1 srcCHA], full kspace
    /// acsDst : calibration data for destination channel [RO E1 dstCHA], full kspace
    /// startRO, endRO, startE1, endE1: define the data region [startRO endRO], [startE1 endE1] which is used for calibration
    /// accelFactor: acceleration factor
    /// thres: the threshold for regularization during kernel estimation
    /// kRO: kernel size along RO
    /// kNE1: kernel size along E1
    /// convKer: computed grappa convolution kernel
    template <typename T> void grappa2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray<T>& convKer);
    /// entire data in acsSrc and acsDst is used
    template <typename T> void grappa2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, hoNDArray<T>& convKer);
    /// dataMask : [RO E1] array, marking fully rectangular sampled region with 1
    template <typename T> void grappa2d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask, size_t accelFactor, double thres, size_t kRO, size_t kNE1, hoNDArray<T>& convKer);

    /// compute image domain kernel from 2d grappd convolution kernel
    /// RO, E1: the size of image domain kernel
    /// kIm: image domain kernel [RO E1 srcCHA dstCHA]
    template <typename T> void grappa2d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, hoNDArray<T>& kIm);

    /// compute unmixing coefficient from image domain kernel and coil sensitivity
    /// kerIm: [RO E1 srcCHA dstCHA], image domain kernel
    /// coilMap: [RO E1 dstCHA] coil sensitivity map
    /// unmixCoeff: [RO E1 srcCHA] unmixing coefficient
    /// gFactor: [RO E1], gfactor
    template <typename T> void grappa2d_unmixing_coeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, size_t acceFactorE1, hoNDArray<T>& unmixCoeff, hoNDArray< typename realType<T>::Type >& gFactor);

    ///  apply kspace domain kernel to unwarp multi-channel images
    /// kspace: [RO E1 srcCHA ... ]
    /// kerIm: [RO E1 srcCHA dstCHA]
    /// complexIm: [RO E1 dstCHA ... ]
    template <typename T> void grappa2d_image_domain_unwrapping(const hoNDArray<T>& kspace, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm);
    /// aliasedIm : [RO E1 srcCHA ...]
    template <typename T> void grappa2d_image_domain_unwrapping_aliased_image(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm);

    /// apply unmixing coefficient on undersampled kspace
    /// kspace: [RO E1 srcCHA ...]
    /// unmixCoeff : [RO E1 srcCHA]
    /// complexIm : [RO E1 ...] wrapped complex images
    template <typename T> void apply_unmix_coeff_kspace(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

    /// aliasedIm : [RO E1 srcCHA ...]
    template <typename T> void apply_unmix_coeff_aliased_image(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

    /// ------------------------
    /// grappa 2d low level functions
    /// ------------------------
    /// get the kernel pattern, given the acceleration factor and kernel size
    /// kE1: kernel pattern along E1
    /// oE1: output pattern along E1
    /// convKRO: convolution kernel size along RO
    /// convKE1: convolution kernel size along E1
    /// e.g. for R=4 and kNE1=4, the kernel pattern kE1 will be [-4 0 4 8] and the output pattern oE1 will be [0 1 2 3] if fitItself==true
    /// if fitItself==false, the output pattern oE1 will be [1 2 3]
    /// if the acsSrc and acsDst are generated in different ways, often fitItself needs to be true; e.g. acsSrc is in the origin acquired channels
    /// and acsDst is in eigen channel
    /// accelFactor: acceleration factor
    /// kRO: kernel size along RO
    /// kNE1: kernel size along E1
    void grappa2d_kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, size_t& convKRO, size_t& convKE1, size_t accelFactor, size_t kRO, size_t kNE1, bool fitItself);

    /// grappa calibration for 2D case
    /// kE1: the kernel pattern along E1
    /// oE1: the output kernel pattern along E1
    /// ker : kernel array [kRO kE1 srcCHA dstCHA oE1]
    /// prepare calibration for A*ker = B
    template <typename T> void grappa2d_prepare_calib(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray<T>& A, hoNDArray<T>& B);

    /// solve for ker
    template <typename T> void grappa2d_perform_calib(const hoNDArray<T>& A, const hoNDArray<T>& B, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, double thres, hoNDArray<T>& ker);

    template <typename T> void grappa2d_calib(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, double thres, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray<T>& ker);

    /// convert the grappa multiplication kernel computed from grappa2d_calib to convolution kernel
    /// convKer : [convRO convE1 srcCHA dstCHA]
    template <typename T> void grappa2d_convert_to_convolution_kernel(const hoNDArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, hoNDArray<T>& convKer);

    /// perform kspace multiplication grappa recon
    /// if periodic_boundary_condition == true, use the periodic boundary condition; other wise, use zero-filling
    /// kspace: [RO E1 srcCHA], every 2D kspace will be reconed to [RO E1 dstCHA]
    /// data matrix A [M kRO*kNE1*srcCHA], using boundary condition for edge cells
    /// data matrix indexes AInd, [M 2], the mth line of data matrix is for oE1 points at (ro, e1)
    template <typename T> void grappa2d_prepare_recon(const hoNDArray<T>& kspace, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, bool periodic_boundary_condition, hoNDArray<T>& A, hoNDArray<unsigned short>& AInd);

    /// perform kspace recon, A*ker
    /// and fill back points to kspace
    template <typename T> void grappa2d_perform_recon(const hoNDArray<T>& A, const hoNDArray<T>& ker, const hoNDArray<unsigned short>& AInd, const std::vector<int>& oE1, size_t RO, size_t E1, hoNDArray<T>& res);

    /// perform ksapce recon
    template <typename T> void grappa2d_recon(hoNDArray<T>& kspace, const hoNDArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, bool periodic_boundary_condition);

    /// in case the recon=A*ker has already been computed, assign them back to res
    template <typename T> void grappa2d_fill_reconed_kspace(const hoNDArray<unsigned short>& AInd, const hoNDArray<T>& recon, const std::vector<int>& oE1, size_t RO, size_t E1, hoNDArray<T>& res);

    /// ---------------------------------------------------------------------
    /// 3D grappa
    /// ---------------------------------------------------------------------

    /// grappa 3d calibration function to compute convolution kernel
    /// acsSrc : calibration data for source channel [RO E1 E2 srcCHA], full kspace
    /// acsDst : calibration data for destination channel [RO E1 E2 dstCHA], full kspace
    /// startRO, endRO, startE1, endE1: define the data region [startRO endRO], [startE1 endE1], [startE2 endE2] which is used for calibration
    /// accelFactorE1, accelFactorE2,: acceleration factor along E1 and E2
    /// thres: the threshold for regularization during kernel estimation
    /// overDetermineRatio: the calibration matrix over determination ratio
    /// kRO: kernel size along RO
    /// kNE1: kernel size along E1
    /// kNE2: kernel size along E2
    /// convKer: computed grappa convolution kernel [convKRO convKE1 convKE2 srcCHA dstCHA]
    template <typename T> void grappa3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            size_t accelFactorE1, size_t accelFactorE2, 
                                                                            double thres, double overDetermineRatio,
                                                                            size_t kRO, size_t kNE1, size_t kNE2, 
                                                                            size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, 
                                                                            hoNDArray<T>& convKer);

    /// entire data in acsSrc and acsDst is used
    template <typename T> void grappa3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                                            size_t accelFactorE1, size_t accelFactorE2,
                                                                            double thres, double overDetermineRatio,
                                                                            size_t kRO, size_t kNE1, size_t kNE2, 
                                                                            hoNDArray<T>& convKer);

    /// dataMask : [RO E1 E2] array, marking fully rectangular sampled region with 1
    template <typename T> void grappa3d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask, 
                                                                            size_t accelFactorE1, size_t accelFactorE2, 
                                                                            double thres, double overDetermineRatio, 
                                                                            size_t kRO, size_t kNE1, size_t kNE2, 
                                                                            hoNDArray<T>& convKer);

    /// compute image domain kernel from 3d grappd convolution kernel
    /// RO, E1, E2: the size of image domain kernel
    /// kIm: image domain kernel [RO E1 E2 srcCHA dstCHA]
    /// if preset_kIm_with_zeros==true, the kIm is cleared to be all zeros inside this function, before the kernel padding
    /// if preset_kIm_with_zeros==false, caller should make sure the kIm is cleared with zeros
    template <typename T> void grappa3d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, size_t E2, hoNDArray<T>& kIm, bool preset_kIm_with_zeros=true);

    /// compute unmixing coefficient from grappa convolution kernel and coil sensitivity
    /// convKer: 3D kspace grappa convolution kernel
    /// coilMap: [RO E1 E2 dstCHA] coil sensitivity map
    /// unmixCoeff: [RO E1 E2 srcCHA] unmixing coefficient
    /// gFactor: [RO E1 E2], gfactor
    template <typename T> void grappa3d_unmixing_coeff(const hoNDArray<T>& convKer, const hoNDArray<T>& coilMap, 
                                                                size_t acceFactorE1, size_t acceFactorE2, 
                                                                hoNDArray<T>& unmixCoeff, 
                                                                hoNDArray< typename realType<T>::Type >& gFactor);

    /// apply grappa convolution kernel to perform per-channel unwrapping
    /// convKer: 3D kspace grappa convolution kernel [convKRO convKE1 convKE2 srcCHA dstCHA]
    /// kspace: undersampled kspace [RO E1 E2 srcCHA] or [RO E1 E2 srcCHA N]
    /// complexIm: [RO E1 E2 dstCHA N] unwrapped complex images
    template <typename T> void grappa3d_image_domain_unwrapping(const hoNDArray<T>& convKer, const hoNDArray<T>& kspace,
                                                                        size_t acceFactorE1, size_t acceFactorE2,
                                                                        hoNDArray<T>& complexIm);

    /// aliasedIm: wrapped complex images [RO E1 E2 srcCHA] or [RO E1 E2 srcCHA N]
    template <typename T> void grappa3d_image_domain_unwrapping_aliasedImage(const hoNDArray<T>& convKer, const hoNDArray<T>& aliasedIm,
                                                                        size_t acceFactorE1, size_t acceFactorE2,
                                                                        hoNDArray<T>& complexIm);

    /// apply unmixing coefficient on undersampled kspace
    /// kspace: [RO E1 E2 srcCHA ...]
    /// unmixCoeff : [RO E1 E2 srcCHA]
    /// complexIm : [RO E1 E2 ...] unwrapped complex images
    template <typename T> void apply_unmix_coeff_kspace_3D(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

    /// aliasedIm : [RO E1 E2 srcCHA ...]
    template <typename T> void apply_unmix_coeff_aliased_image_3D(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

    /// ------------------------
    /// grappa 3d low level functions
    /// ------------------------
    /// get the kernel pattern, given the acceleration factor and kernel size
    /// kE1, kE2: kernel pattern along E1 and E2
    /// oE1, oE2: output pattern along E1 and E2
    /// convKRO: convolution kernel size along RO
    /// convKE1: convolution kernel size along E1
    /// convKE2: convolution kernel size along E2
    /// e.g. for R=4 and kNE1=kNE2=4, the kernel pattern kE1 and kE2 will be [-4 0 4 8] and the output pattern oE1 and oE2 will be [0 1 2 3] if fitItself==true
    /// if fitItself==false, the output pattern oE1 and oE2 will be [1 2 3]
    /// if the acsSrc and acsDst are generated in different ways, often fitItself needs to be true; e.g. acsSrc is in the origin acquired channels
    /// and acsDst is in eigen channel
    /// accelFactor: acceleration factor
    /// kRO: kernel size along RO
    /// kNE1: kernel size along E1
    /// kNE2: kernel size along E2
    void grappa3d_kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, 
                                        std::vector<int>& kE2, std::vector<int>& oE2, 
                                        size_t& convKRO, size_t& convKE1, size_t& convKE2, 
                                        size_t accelFactorE1, size_t accelFactorE2, 
                                        size_t kRO, size_t kNE1, size_t kNE2, bool fitItself);

    /// grappa calibration for 3D case
    /// kE1: the kernel pattern along E1
    /// kE2: the kernel pattern along E2
    /// oE1: the output kernel pattern along E1
    /// oE2: the output kernel pattern along E2
    /// ker : kernel array [kRO kE1 kE2 srcCHA dstCHA oE1 oE2]
    template <typename T> void grappa3d_calib(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, 
                                                    double thres, double overDetermineRatio, size_t kRO, 
                                                    const std::vector<int>& kE1, const std::vector<int>& oE1, 
                                                    const std::vector<int>& kE2, const std::vector<int>& oE2, 
                                                    hoNDArray<T>& ker);

    /// convert the grappa multiplication kernel computed from grappa3d_calib to convolution kernel
    /// convKer : [convRO convE1 convE2 srcCHA dstCHA]
    template <typename T> void grappa3d_convert_to_convolution_kernel(const hoNDArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, const std::vector<int>& kE2, const std::vector<int>& oE2, hoNDArray<T>& convKer);

}
