
/** \file   mri_core_grappa.h
    \brief  GRAPPA implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_parallel_imaging.h"
#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"
#include "ho6DArray.h"
#include "ho7DArray.h"

namespace Gadgetron
{

template <typename T> 
class EXPORTMRICORE grappa : public parallelImaging<T>
{
public:

    typedef typename realType<T>::Type value_type;

    grappa(){}
    virtual ~grappa() {}

    virtual void printInfo(std::ostream& os);

    /// ---------------------------------------------------------------------
    /// 2D grappa
    /// ---------------------------------------------------------------------

    /// get the kernel pattern, given the acceleration factor and kernel size
    /// kE1: kernel pattern along E1
    /// oE1: output pattern along E1
    /// e.g. for R=4, the kernel size is 5*4, the kernel pattern kE1 will be [-4 0 4 8] and the output pattern oE1 will be [0 1 2 3] if fitItself==true
    /// if fitItself==false, the output pattern oE1 will be [1 2 3]
    /// if the acsSrc and acsDst are generated in different ways, often fitItself needs to be true; e.g. acsSrc is in the origin acquired channels
    /// and acsDst is in eigen channel
    /// accelFactor: acceleration factor
    /// kNE1: kernel size along E1
    void kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, size_t accelFactor, size_t kNE1, bool fitItself);

    /// grappa calibration for 2D case
    /// acsSrc : calibration data for source channel [RO E1 srcCHA]
    /// acsDst : calibration data for destination channel [RO E1 dstCHA]
    /// thres: the threshold for regularization during kernel estimation
    /// kRO: kernel size along RO
    /// kE1: the kernel pattern along E1
    /// oE1: the output kernel pattern along E1
    /// ker : kernel array [kRO kE1 srcCHA dstCHA oE1]
    void calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, 
            size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray<T>& ker);

    /// compute image domain kernel for 2D kernel
    /// ro, e1: the size of image domain kernel
    /// kIm: image domain kernel [ro e1 srcCHA dstCHA]
    void imageDomainKernel(const ho5DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t ro, size_t e1, hoNDArray<T>& kIm);

    /// ---------------------------------------------------------------------
    /// 3D grappa
    /// ---------------------------------------------------------------------
    /// grappa calibration for 3D case
    /// acsSrc : [RO E1 E2 srcCHA]
    /// acsDst : [RO E1 E2 dstCHA]
    /// overDetermineRatio: the kernel calibration overDetermineRatio determines how many equations will be used for kernel calculation
    /// ker : [kRO kE1 kE2 srcCHA dstCHA oE1 oE2]
    void calib3D(const ho4DArray<T>& acsSrc, const ho4DArray<T>& acsDst, double thres, double overDetermineRatio, 
            size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, ho7DArray<T>& ker);

    /// image domain kernel for 3D kernel
    /// ro, e1, e2: the size of image domain kernel
    /// kIm: image domain kernel [ro e1 e2 srcCHA dstCHA]
    void imageDomainKernel3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, size_t ro, size_t e1, size_t e2, hoNDArray<T>& kIm);

    /// ---------------------------------------------------------------------
    /// 3D grappa with E1 and E2 in kspace, for readout decoupled reconstruction workflow
    /// ---------------------------------------------------------------------

    /// compute hybrid image domain kernel for 3D kernel, only RO direction is converted to image domain
    /// sometime we want to reduce the peak memory usage of 3D reconstruction, therefore, only RO dimension is converted to image domain
    /// then the data will be unwrapped along every chunk of RO dimenions, e.g. if RO=256, the data is split into 5 chunks, then every time 50 2D slices are unwrapped
    /// E1 and E2 stays in the kspace domain
    /// kImRO: kspace-image hybrid kernel [convE1 convE2 RO srcCHA dstCHA]
    void imageDomainKernelRO3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, size_t ro, hoNDArray<T>& kImRO);

    /// image domain kernel for 3D kernel, E1 and E2 directions are converted to image domain
    /// kImRO : kspace-image hybrid kernel where first two dimensions are E1 and E2 and in kspace
    void imageDomainKernelE1E2RO(const hoNDArray<T>& kImRO, size_t e1, size_t e2, hoNDArray<T>& kImE1E2RO);

protected:

    /// convert the calibrated kernel to the convlution kernel in kspace
    /// if ROis3rdDim == true, the kernel dimension is [E1 E2 RO], otherwise [RO E1 E2]
    void kspaceDomainConvKernel3D(const ho7DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& kE2, const std::vector<int>& oE1, const std::vector<int>& oE2, ho5DArray<T>& convKerFlip, bool ROis3rdDim=true);
};

}
