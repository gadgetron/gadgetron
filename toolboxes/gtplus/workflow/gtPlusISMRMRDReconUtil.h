/** \file   gtPlusISMRMRDReconUtil.h
    \brief  Define the symbols and implement common functionalities for GtPlus toolbox

            The ISMRMRD format is fully supported in this toolbox.

            Other functinalities implemented here include:
            Karhunen-Lo�ve Transform (KLT) or Principle Component Analysis (PCA)
            KSpace filter
            Several MR sensitivity map estimation methods

            Ref to :

            http://en.wikipedia.org/wiki/Karhunen%E2%80%93Lo%C3%A8ve_theorem

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

#include "GtPlusExport.h"

#include "ismrmrd/ismrmrd.h"
#include "ismrmrd/meta.h"

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_comparison.hpp"
#include "boost/tuple/tuple_io.hpp"

#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"
#include "ho6DArray.h"
#include "ho7DArray.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"
#include "hoNDArray_reductions.h"
#include "hoNDArray_linalg.h"
#include "gtPlusIOAnalyze.h"
#include "hoNDArrayMemoryManaged.h"
#include "GadgetronTimer.h"
#include "log.h"

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

#include "GtPlusDefinition.h"

// MACROS for EXPORTING
#define GADGET_EXPORT_ARRAY(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportArray(a, debugFolder+filename); } }
#define GADGET_EXPORT_ARRAY_COMPLEX(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportArrayComplex(a, debugFolder+filename); } }
#define GADGET_EXPORT_ARRAY_COMPLEX_REAL_IMAG(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportArrayComplexRealImag(a, debugFolder+filename); } }

#define GADGET_EXPORT_IMAGE(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportImage(a, debugFolder+filename); } }
#define GADGET_EXPORT_IMAGE_COMPLEX(debugFolder, exporter, a, filename) { if ( !debugFolder.empty() ) { exporter.exportImageComplex(a, debugFolder+filename); } }

namespace Gadgetron {

    /**
    * @brief sum over last dimension of an array
             e.g. for a 4D array, sum over the 4th dimension and get a 3D array
    */
    template<typename T> EXPORTGTPLUS bool sumOverLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r); // 

    /**
    * @brief sum over the second last dimension of an array
             e.g. for a 4D array, sum over the 3rd dimension and get a 3D array
    */
    template<typename T> EXPORTGTPLUS bool sumOverSecondLastDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief multiply over the last dimension of y by x
             e.g. x is 3D and y is 4D array, r(:,:,:,n) = y(:,:,:,n) .* x
    */
    template<typename T> EXPORTGTPLUS bool multiplyOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief divide the last dimension of y by x
             e.g. x is 3D and y is 4D array, r(:,:,:,n) = y(:,:,:,n) ./ x
    */
    template<typename T> EXPORTGTPLUS bool divideOverLastDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief sum over the 1st dimension of an array
             e.g. for a 2D array, sum over the 1st dimension and get an array of [1 E1]
    */
    template<typename T> EXPORTGTPLUS bool sumOver1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief sum over the 2nd dimension of an array
             e.g. for a 3D array, sum over the 2nd dimension and get an array of [RO 1 CHA]
    */
    template<typename T> EXPORTGTPLUS bool sumOver2ndDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief sum over the 3rd dimension of an array
             e.g. for a 4D array, sum over the 3rd dimension and get an array of [RO E1 1 N]
    */
    template<typename T> EXPORTGTPLUS bool sumOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief sum over the 4th dimension of an array
             e.g. for a 5D array [RO E1 CHA N S], sum over the 4th dimension and get an array of [RO E1 CHA 1 S]
    */
    template<typename T> EXPORTGTPLUS bool sumOver4thDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief sum over the 5th dimension of an array
             e.g. for a 6D array, sum over the 5th dimension and get an array [RO E1 CHA N 1 P]
    */
    template<typename T> EXPORTGTPLUS bool sumOver5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief multiply over the 3rd/4th/5th dimension of y by x
             e.g. x is 3D and y is 4D array, r(:,:,n,:) = y(:,:,n,:) .* x
             e.g. x is 4D and y is 5D array, r(:,:,:,n,:) = y(:,:,:,n,:) .* x
             e.g. x is 5D and y is 6D array, r(:,:,:,:, n,:) = y(:,:,:,:,n,:) .* x
    */
    template<typename T> EXPORTGTPLUS bool multiplyOver3rdDimension(const hoNDArray<T>& x3D, const hoNDArray<T>& y4D, hoNDArray<T>& r);
    template<typename T> EXPORTGTPLUS bool multiplyOver4thDimension(const hoNDArray<T>& x4D, const hoNDArray<T>& y5D, hoNDArray<T>& r);
    template<typename T> EXPORTGTPLUS bool multiplyOver5thDimension(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief multiply over the 4th/5th dimension of y by x except for dimension index n
             e.g. x is 4D and y is 5D array, r(:,:,:,t,:) = y(:,:,:,t,:) .* x, except for r(:,:,:,n,:) = y(:,:,:,n,:)
             e.g. x is 5D and y is 6D array, r(:,:,:,:,t,:) = y(:,:,:,:,t,:) .* x, except for r(:,:,:,:,n,:) = y(:,:,:,:,n,:)
    */
    template<typename T> EXPORTGTPLUS bool multiplyOver4thDimensionExcept(const hoNDArray<T>& x4D, const hoNDArray<T>& y5D, size_t n, hoNDArray<T>& r, bool copyY2R=true);
    template<typename T> EXPORTGTPLUS bool multiplyOver5thDimensionExcept(const hoNDArray<T>& x, const hoNDArray<T>& y, size_t n, hoNDArray<T>& r, bool copyY2R=true);

    /**
    * @brief r = x add/multiply/divide y for every part of y
    */
    template<typename T> EXPORTGTPLUS bool multipleAdd(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);
    template<typename T> EXPORTGTPLUS bool multipleMultiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);
    template<typename T> EXPORTGTPLUS bool multipleDivide(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r);

    /**
    * @brief copy the sub-array of x to r
             the sub-array is defined by its starting index and array size
    */
    template<typename T> EXPORTGTPLUS bool cropUpTo11DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const std::vector<size_t>& startND, std::vector<size_t>& size);

    /**
    * @brief set the sub-array of r from x
             the sub-array is defined by its starting index and array size
    */
    template<typename T> EXPORTGTPLUS bool setSubArrayUpTo11DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const std::vector<size_t>& startND, std::vector<size_t>& size);

    /**
    * @brief extract sampled lines from an NDArray
             timeStamp indicates sampled lines; -1 for unsampled lines
             x : [Ro E1 Cha Slice E2 Con Phase Rep Set Seg AVE]
             timeStamp: [1 E1 1 Slice E2 Con Phase Rep Set Seg AVE]
    */
    template<typename T> EXPORTGTPLUS bool extractSampledLinesUpTo11DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const hoNDArray<float>& timeStamp, double acceFactorE1, double acceFactorE2);

    /**
    * @brief fill sampled lines to an NDArray
             timeStamp indicates sampled lines; -1 for unsampled lines
    */
    template<typename T> EXPORTGTPLUS bool fillSampledLinesUpTo11DArray(const hoNDArray<T>& x, hoNDArray<T>& r, const hoNDArray<float>& timeStamp);

    /**
    * @brief copy the sub-array of x to r only along the 3rd dimensions
             e.g. x is [RO E1 D3 ...], r will be [RO E1 end-start+1 ... ]
    */
    template<typename T> EXPORTGTPLUS bool cropOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end);

    /**
    * @brief set the sub-array of r from x only along the 3rd dimensions
             e.g. r(:, :, start:end, :, ...) will be replaced by x
    */
    template<typename T> EXPORTGTPLUS bool setSubArrayOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t start, size_t end);

    /**
    * @brief compute the standard deviation along the 3rd dimension, if NMinusOne == true, divided by N-1; otherwise, divided by N
    */
    template<typename T> EXPORTGTPLUS bool stdOver3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& std, bool NMinusOne);

    /**
    * @brief permute E2 dimension of x : [RO E1 CHA SLC E2 ...] to r: [RO E1 E2 CHA SLC ...]
    */
    template<typename T> EXPORTGTPLUS bool permuteE2To3rdDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute E2 dimension of x : [RO E1 E2 CHA SLC ...] to r: [RO E1 CHA SLC E2 ...]
    */
    template<typename T> EXPORTGTPLUS bool permuteE2To5thDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute RO dimension of x to the 3rd dimension
             x : [RO E1 E2 ...], r: [E1 E2 RO ...]
    */
    template<typename T> EXPORTGTPLUS bool permuteROTo3rdDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute RO dimension of x to the 4th dimension
             x : [RO E1 E2 CHA ...], r: [E1 E2 CHA RO ...]
    */
    template<typename T> EXPORTGTPLUS bool permuteROTo4thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute RO dimension of x back to the 1st dimension
             x : [E1 E2 CHA RO ...], r: [RO E1 E2 CHA ...]
    */
    template<typename T> EXPORTGTPLUS bool permuteROTo1stDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute the 3rd dimension of x to the 1st dimension
             x : [RO E1 E2 CHA ...], r: [E2 RO E1 CHA ...]
    */
    template<typename T> EXPORTGTPLUS bool permute3rdDimensionTo1stDimension(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief permute RO dimension of x to the 5th dimension
             x : [RO E1 E2 srcCHA dstCHA ...], r: [E1 E2 srcCHA dstCHA RO ...]
    */
    template<typename T> EXPORTGTPLUS bool permuteROTo5thDimensionFor3DRecon(const hoNDArray<T>& x, hoNDArray<T>& r);

    /**
    * @brief Image domain unwrapping for 2D
             x : [RO E1 srcCHA], ker [RO E1 srcCHA dstCHA]
             buf is a buffer for computer, need to be pre-allocated [RO E1 srcCHA], y [RO E1 dstCHA]
             for the sake of speed, no check is made in this function
    */
    template<typename T> EXPORTGTPLUS bool imageDomainUnwrapping2D(const hoNDArray<T>& x, const hoNDArray<T>& ker, hoNDArray<T>& buf, hoNDArray<T>& y);

    /**
    * @brief Image domain unwrapping for 2D
             x : [RO E1 srcCHA N], ker [RO E1 srcCHA dstCHA 1 or N], 
             buf is a buffer for computer, need to be pre-allocated [RO E1 srcCHA], y [RO E1 dstCHA N]
             for the sake of speed, no check is made in this function
    */
    template<typename T> EXPORTGTPLUS bool imageDomainUnwrapping2DT(const hoNDArray<T>& x, const hoNDArray<T>& ker, hoNDArray<T>& buf, hoNDArray<T>& y);

    /**
    * @brief compute periodic boundary values for an array
             x : [N 1] the data point location, y[N M] data point values at x
             r : [N+2 M], the data point values with computed boundaries
    */
    template<typename CoordType, typename T> EXPORTGTPLUS bool computePeriodicBoundaryValues(const hoNDArray<CoordType>& x, const hoNDArray<T>& y, CoordType start, CoordType end, hoNDArray<CoordType>& vx, hoNDArray<T>& vy);
}

namespace Gadgetron { namespace gtPlus {

// ================================================================================================== //

template <typename T> 
class gtPlusISMRMRDReconUtil
{
public:

    typedef typename realType<T>::Type value_type;

    gtPlusISMRMRDReconUtil();
    virtual ~gtPlusISMRMRDReconUtil();

    void printInfo(std::ostream& os);

    typedef std::pair<ISMRMRDDIM, size_t> DimensionRecordType;

    // ------------------------------------------------------------------------
    // coil compression and KarhunenLoeverTransform
    // ------------------------------------------------------------------------
    // data: M rows and N cols matrix
    // the KLT direction is along the N
    // eigenVectors: N*N eigen vectors, every column is a eigen vector
    // eigenValues: N*1 eigen values, descending order
    bool KLT_eigenAnalysis(const hoMatrix<T>& data, hoMatrix<T>& eigenVectors, hoMatrix<T>& eigenValues);

    // apply the eigen transform
    // data: M*N data matrix
    // eigenVectors: N*K eigen vector matrix, every column is a eigen vector
    // dataEigen: M*K eigen data matrix
    bool KLT_applyEigen(const hoMatrix<T>& data, hoMatrix<T>& dataEigen, const hoMatrix<T>& eigenVectors);
    bool KLT_applyEigen(const hoNDArray<T>& data, hoNDArray<T>& dataEigen, const hoMatrix<T>& eigenVectors);

    // number of kept eigen modes
    // all modes with eigen values greater than thres*max(eigenValues) are kept
    bool KLT_numberOfKeptModes(const hoMatrix<T>& eigenValues, double thres, long long& numOfModesKept);

    // prune the eigen vector matrixes to keep the last numOfModesKept columns
    bool pruneEigenVectorMatrix(const hoMatrix<T>& eigenVectors, long long numOfModesKept, hoMatrix<T>& eigenVectorsPruned);

    // KLT based coil compression
    // data: at least 3D [RO E1 CHA ...]
    // the KL transform is applied along CHA
    // coeff: CHA*numOfModesKept eigen vector matrix
    // eigenValues: CHA*1 eigen values
    // thres <0 or numOfModesKept==-1, keep all modes
    // if isChaLastDim==true, the CHA is the last dimension
    bool computeKLCoilCompressionCoeff(const hoNDArray<T>& data, double thres, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, bool isChaLastDim=false);
    bool computeKLCoilCompressionCoeff(const hoNDArray<T>& data, int numOfModesKept, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, bool isChaLastDim=false);
    // coeff: CHA*CHA eigen vector matrix
    bool computeKLTCoeff(const hoNDArray<T>& data, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, bool isChaLastDim=false);

    // dataEigen: [RO E1 numOfModesKept ...] 
    bool computeKLCoilCompression(const hoNDArray<T>& data, double thres, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, hoNDArray<T>& dataEigen, bool isChaLastDim=false);
    bool computeKLCoilCompression(const hoNDArray<T>& data, int numOfModesKept, hoMatrix<T>& coeff, hoMatrix<T>& eigenValues, hoNDArray<T>& dataEigen, bool isChaLastDim=false);

    // apply coil compression coefficients
    bool appyKLCoilCompressionCoeff(const hoNDArray<T>& data, const hoMatrix<T>& coeff, hoNDArray<T>& dataEigen, bool isChaLastDim=false);

    // apply coil compression coefficients on array [RO E1 srcCHA ...]
    // dataEigen: [RO E1 dstCHA ...]
    // coeff: [srcCHA dstCHA] matrixes for every last dimension
    // every last dimension has different compression coefficients
    bool applyKLCoilCompressionCoeff(const hoNDArray<T>& data, const std::vector<hoMatrix<T> >& coeff, hoNDArray<T>& dataEigen, bool isChaLastDim=false);

    // compute KL transform and perform filtering
    // the KL dimension is the last dimension
    bool computeKLFilter(const hoNDArray<T>& data, size_t numOfModesKept, hoNDArray<T>& dataKLF);

    // ------------------------------------------------------------------------
    // zero-padding resize
    // ------------------------------------------------------------------------
    // compute the start and end index for zero padding
    // dstSize >= srcSize
    bool zpadRange(size_t srcSize, size_t dstSize, size_t& start, size_t& end);

    // pad the first two dimensions around its center, other dimensions are kept unchanged
    bool zeropad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataPadded);

    // pad first three dimensions array around its center, other dimensions are kept unchanged
    bool zeropad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded);

    // the dataPadded is not pre cleared to fill with zeros
    bool zeropad3DNoPresetZeros(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataPadded);

    // cut the center part
    bool cutpad2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataCut);
    bool cutpad3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataCut);

    // ------------------------------------------------------------------------
    // kspace filter
    // ------------------------------------------------------------------------
    bool compute2DFilterFromTwo1D(const hoNDArray<T>& fx, const hoNDArray<T>& fy, hoNDArray<T>& fxy);
    bool compute2DFilterFromTwo1D(const hoNDArray<float>& fx, const hoNDArray<float>& fy, hoNDArray< std::complex<float> >& fxy);
    bool compute2DFilterFromTwo1D(const hoNDArray<double>& fx, const hoNDArray<double>& fy, hoNDArray< std::complex<double> >& fxy);

    bool compute3DFilterFromThree1D(const hoNDArray<T>& fx, const hoNDArray<T>& fy, const hoNDArray<T>& fz, hoNDArray<T>& fxyz);
    bool compute3DFilterFromThree1D(const hoNDArray<float>& fx, const hoNDArray<float>& fy, const hoNDArray<float>& fz, hoNDArray< std::complex<float> >& fxyz);
    bool compute3DFilterFromThree1D(const hoNDArray<double>& fx, const hoNDArray<double>& fy, const hoNDArray<double>& fz, hoNDArray< std::complex<double> >& fxyz);

    // data: in kspace, [RO E1 E2 CHA SLC CON PHS REP SET]
    bool kspacefilterRO(hoNDArray<T>& data, const hoNDArray<T>& fRO);
    bool kspacefilterRO(const hoNDArray<T>& data, const hoNDArray<T>& fRO, hoNDArray<T>& dataFiltered);
    bool kspacefilterROE1(const hoNDArray<T>& data, const hoNDArray<T>& fROE1, hoNDArray<T>& dataFiltered);
    bool kspacefilterROE1(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered);
    bool kspacefilterE1(const hoNDArray<T>& data, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered);

    // kspace fitler for ISMRMRD dimension order
    bool kspacefilterE2(const hoNDArray<T>& data, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    bool kspacefilterROE2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    bool kspacefilterE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    bool kspacefilterROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fROE1E2, hoNDArray<T>& dataFiltered);
    bool kspacefilterROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);

    // kspace filter for the array whose first three dimensions are RO, E1 and E2; 
    bool kspace3DfilterE2(const hoNDArray<T>& data, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    bool kspace3DfilterROE2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    bool kspace3DfilterE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    bool kspace3DfilterROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fROE1E2, hoNDArray<T>& dataFiltered);
    bool kspace3DfilterROE1E2(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);

    // ------------------------------------------------------------------------
    // generate kspace filters
    // ------------------------------------------------------------------------
    // symmetric filter, used for image filtering
    // sigma: for Gaussian, in the unit of pixel
    // width: for Tukey filter etc., the length of transition band
    bool generateSymmetricFilter(size_t len, size_t start, size_t end, hoNDArray<T>& filter, ISMRMRDKSPACEFILTER filterType, double sigma, size_t width);

    // asymmetric filter, used for partial fourier/asymmetric echo filtering
    // start, end: the data range
    // tapered hanning filer is implemented for this
    // if filterType==ISMRMRD_FILTER_NONE and densityComp==true, the 0-1-2 filter will be generated
    // if filterType==ISMRMRD_FILTER_TAPERED_HANNING and the densityComp is true, the density compensation version of tapered filter will be generated
    // where unacquired region has filter values 0 and symmetric region 1 and nonsymmetric region 2
    // if densityComp==false, the one side tapered filter will be generated
    bool generateAsymmetricFilter(size_t len, size_t start, size_t end, hoNDArray<T>& filter, ISMRMRDKSPACEFILTER filterType, size_t width, bool densityComp=false);

    // generate ref data filter
    bool generateSymmetricFilterForRef(size_t len, size_t start, size_t end, hoNDArray<T>& filter, ISMRMRDKSPACEFILTER filterType, double sigma, size_t width);

    // find the symmetric sampled region
    bool findSymmetricSampledRegion(size_t start, size_t end, size_t center, size_t& startSym, size_t& endSym);

    // compute the filter SNR unit scale factor
    bool computeFilterSNRUnitScaleFactor(const hoNDArray<T>& filter, T& scalFactor);

    // ------------------------------------------------------------------------
    // detect sampled region
    // ------------------------------------------------------------------------
    // data : [RO E1 SLC E2 CON PHS REP SET] array
    bool detectSampledRegion2D(const hoNDArray<T>& data, size_t& startRO, size_t& endRO, size_t& startE1, size_t& endE1);
    bool detectSampledRegion3D(const hoNDArray<T>& data, size_t& startRO, size_t& endRO, size_t& startE1, size_t& endE1, size_t& startE2, size_t& endE2);

    // ------------------------------------------------------------------------
    // coil sensitivity
    // ------------------------------------------------------------------------
    // average kspace along the 4th dimension
    // data: [RO E1 CHA N S ...]
    // ave: [RO E1 CHA 1 S ... ]
    // simple average
    bool averageKSpace4D(const hoNDArray<T>& data, hoNDArray<T>& ave);
    // the sampled times are considered for averaging for E1 dimension
    // sampledTimes: [E1 1], recording the number of sampled times for each lines
    bool averageKSpace4D(const hoNDArray<T>& data, hoNDArray<T>& ave, std::vector<size_t>& sampledTimes);

    // average kspace along the 5th dimension
    // data: [RO E1 E2 CHA N ...]
    // ave: [RO E1 E2 CHA 1 ... ]
    // simple average
    bool averageKSpace5D(const hoNDArray<T>& data, hoNDArray<T>& ave);
    // the sampled times are considered for averaging for E1 and E2 dimension
    // sampledTimes: [E1 E2], recording the number of sampled times for each lines
    bool averageKSpace5D(const hoNDArray<T>& data, hoNDArray<T>& ave, hoNDArray<size_t>& sampledTimes);

    // sampled region along E1
    // data: [RO E1 CHA N]
    bool detectSampledRegionE1(const hoNDArray<T>& data, size_t& startE1, size_t& endE1);

    // sampled times along E1, if not sampled, sampledTimes[e1] == 0
    bool detectSampledTimesE1(const hoNDArray<T>& data4D, std::vector<size_t>& sampledTimes);

    // sampled region along E1 and E2
    // data: [RO E1 E2 CHA N]
    bool detectSampledRegionE1E2(const hoNDArray<T>& data, size_t& startE1, size_t& endE1, size_t& startE2, size_t& endE2);

    // sampled times along E1 and E2, if not sampled, sampledTimes(e1, e2) == 0
    // data5D: [RO E1 E2 CHA N]
    bool detectSampledTimesE1E2(const hoNDArray<T>& data5D, hoNDArray<size_t>& sampledTimes);

    // copy along E1
    bool copyAlongE1(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startE1, size_t endE1);

    // copy along RO and E1
    bool copyAlongROE1(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1);

    // copy along RO, E1 and E2
    bool copyAlongROE1E2(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2);

    // copy along RO and E1, but a transition band is used to make sure the smoothing transition on the dst kspace
    // the transition band is achieved via the tapered hanning filter
    bool copyAlongROE1TransitionBand(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, 
            size_t startE1, size_t endE1, size_t transBandRO, size_t transBandE1);

    // copy along RO, E1 and E2, but a transition band is used to make sure the smoothing transition on the dst kspace
    // the transition band is achieved via the tapered hanning filter
    // src, dst: [RO E1 E2 ...]
    bool copyAlongROE1E2TransitionBand(const hoNDArray<T>& src, hoNDArray<T>& dst, 
                                    size_t startRO, size_t endRO, 
                                    size_t startE1, size_t endE1, 
                                    size_t startE2, size_t endE2, 
                                    size_t transBandRO, size_t transBandE1, 
                                    size_t transBandE2);

    // ------------------------------------------------------------------------
    // ISMRMRDDIM related functions
    // ------------------------------------------------------------------------
    // get the dimension name
    std::string getISMRMRDDimName(const ISMRMRDDIM& dim);
    ISMRMRDDIM getISMRMRDDimFromName(const std::string& name);

    // get the dimension order index in the ISMRMRD format
    // this function is for the kspace
    //  [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
    //   0  1  2   3     4  5   6     7   8   9
    bool getISMRMRDDimIndex(const ISMRMRDDIM& dim, long long& ind);

    // find the dimension indexes
    bool findDimIndex(const std::vector<DimensionRecordType>& dimStartingIndexes, ISMRMRDDIM dim, size_t ind);

    // get recon algorithm from name
    ISMRMRDALGO getISMRMRDReconAlgoFromName(const std::string& name);

    // get coil map algorithm from name
    ISMRMRDCOILMAPALGO getISMRMRDCoilMapAlgoFromName(const std::string& name);

    // get the partial fourier/asymmetric echo handling algorithm from name
    ISMRMRDPFALGO getISMRMRDPartialFourierReconAlgoFromName(const std::string& name);

    // get the partial fourier/asymmetric echo handling algorithm name from algorithm
    std::string getNameFromISMRMRDPartialFourierReconAlgo(ISMRMRDPFALGO algo);

    // get the kspace filter algorithm from name
    ISMRMRDKSPACEFILTER getISMRMRDKSpaceFilterFromName(const std::string& name);

    // get retro-gating interpolation method from name
    ISMRMRDINTERPRETROGATING getISMRMRDRetroGatingInterpFromName(const std::string& name);

    // extract sub array for a dimension
    // if lessEqual ==  true, [0:value] are extracted for dim
    bool extractSubArrayForDim(const hoNDArray<T>& x, hoNDArray<T>& r, ISMRMRDDIM& dim, size_t value, bool lessEqual);
    // if lessEqual ==  true, [0:value1 0:value2] are extracted for dim
    bool extractSubArrayForDim(const hoNDArray<T>& x, hoNDArray<T>& r, ISMRMRDDIM& dim1, size_t value1, ISMRMRDDIM& dim2, size_t value2, bool lessEqual);

    // extract sub array for two dimensions
    // [0:value1] and [value2]
    bool extractSubArrayForDim1LessEqualDim2Equal(const hoNDArray<T>& x, hoNDArray<T>& r, ISMRMRDDIM& dim1, size_t value1, ISMRMRDDIM& dim2, size_t value2);

    // extract sub array limited by the max encoding counters
    bool extractSubArrayForMaxEncodingCounters(const hoNDArray<T>& x, hoNDArray<T>& r, const ISMRMRD::EncodingCounters& maxIdx);

    // ------------------------------------------------------------------------
    // ISMRMRD acquisition header
    // ------------------------------------------------------------------------
    void clearAcquisitionHeaderISMRMRD(ISMRMRD::AcquisitionHeader& acqHeader);

    // compute the image geometry for two acquisition header
    bool hasIdenticalGeometryISMRMRD(const ISMRMRD::AcquisitionHeader& acqHeader1, const ISMRMRD::AcquisitionHeader& acqHeader2);

    // add zeros pre/post data array
    // 1 : pre zeros
    // 2 : post zeros
    // 0 : no zeros
    long long addPrePostZeros(size_t centre_column, size_t samples);

    // find RO ranges from centre_column and number of samples
    void findStartEndRO(size_t centre_column, size_t samples, long long& startRO, long long& endRO);

    // find RO ranges from centre_column and number of samples after zero-filling
    void findStartEndROAfterZeroFilling(size_t centre_column, size_t samples_zerofilled, int& startRO, int& endRO);

    // ------------------------------------------------------------------------
    // ISMRMRD image header
    // ------------------------------------------------------------------------
    // set the meta attributes from the ISMRMRD image header
    bool setMetaAttributesFromImageHeaderISMRMRD(const ISMRMRD::ImageHeader& imgHeader, ISMRMRD::MetaContainer& attrib);

    // compute the image geometry for two acquisition header
    bool setImageHeaderISMRMRDFromMetaAttributes(const ISMRMRD::MetaContainer& attrib, ISMRMRD::ImageHeader& imgHeader);

    // ------------------------------------------------------------------------
    // utility functions for various things
    // ------------------------------------------------------------------------
    // jobSchedule : for every valid device, it records the job allocated to it
    // what is stored are valid device id and job packages allocated to it
    // for one valid device, multiple job packages can be given to it

    // load two hoNDArray and compute differences
    void compareAgainstGroundTruthArray(const std::string& gt_filename, const hoNDArray<T>& x, typename realType<T>::Type& normDiff, typename realType<T>::Type& maxNormDiff);
    void compareAgainstGroundTruthArray(const hoNDArray<T>& gt, const hoNDArray<T>& x, typename realType<T>::Type& normDiff, typename realType<T>::Type& maxNormDiff);

    template <typename T2, unsigned int D> void compareAgainstGroundTruthImage(const std::string& gt_filename, const hoNDImage<T2, D>& x, typename realType<T2>::Type& normDiff, typename realType<T2>::Type& maxNormDiff)
    {
        hoNDImage<T2, D> gt;

        gtPlusIOAnalyze gt_io;
        gt_io.importImage(gt, gt_filename);

        compareAgainstGroundTruthImage(gt, x, normDiff, maxNormDiff);
    }

    template <typename T2, unsigned int D> void compareAgainstGroundTruthImage(const hoNDImage<T2, D>& gt, const hoNDImage<T2, D>& x, typename realType<T2>::Type& normDiff, typename realType<T2>::Type& maxNormDiff)
    {
        hoNDImage<T2, D> diff(x);
        Gadgetron::subtract(gt, x, diff);

        hoNDImage<T2, D> gtEps(gt);
        Gadgetron::addEpsilon(gtEps);

        Gadgetron::norm2(diff, normDiff);

        Gadgetron::divide(diff, gtEps, diff);

        T2 maxV;
        size_t ind;
        Gadgetron::maxAbsolute(diff, maxV, ind);
        maxNormDiff = std::abs(maxV);
    }

    void getCurrentMoment(std::string& procTime)
    {
        char timestamp[100];
        time_t mytime;
        struct tm *mytm;
        mytime=time(NULL);
        mytm=localtime(&mytime);
        strftime(timestamp, sizeof(timestamp),"%a, %b %d %Y, %H:%M:%S",mytm);
        procTime = timestamp;
    }

    void getCurrentMomentForFileName(std::string& procTime)
    {
        char timestamp[100];
        time_t mytime;
        struct tm *mytm;
        mytime=time(NULL);
        mytm=localtime(&mytime);
        strftime(timestamp, sizeof(timestamp),"%a_%b_%d_%Y_%H_%M_%S",mytm);
        procTime = timestamp;
    }
};

// utility functions only meaningful for complex data type
template <typename T> 
class gtPlusISMRMRDReconUtilComplex : public gtPlusISMRMRDReconUtil<T>
{
public:

    typedef typename realType<T>::Type value_type;

    gtPlusISMRMRDReconUtilComplex();
    virtual ~gtPlusISMRMRDReconUtilComplex();

    void printInfo(std::ostream& os);

    // ------------------------------------------------------------------------
    // noise prewhitening
    // ------------------------------------------------------------------------
    // compute the noise prewhitening matrix
    // noise: the noise scan [RO E1 CHA]
    // noiseBandWidth: the noise bandwidth, Hz/pixel
    // receiverBWRatio: system receiver noise equivaluent bandwidth ratio
    // ADCSamplingTimeinSecond: ADC sampling time in second
    // prewhiteningMatrix: the computed noise prewhitening matrix [CHA CHA]
    bool computeNoisePrewhiteningMatrix(const hoNDArray<T>& noise, double noiseBandWidth, double receiverBWRatio, double ADCSamplingTimeinSecond, hoMatrix<T>& prewhiteningMatrix);

    // perform the noise prewhitening matrix on the image/ref data
    // result = prewhiteningMatrix * data
    // data should at least have three dimensions [R0 E1 CHA], up to 10D
    bool performNoisePrewhitening(hoNDArray<T>& data, const hoMatrix<T>& prewhiteningMatrix);

    // ------------------------------------------------------------------------
    // zero-padding resize
    // ------------------------------------------------------------------------
    // zero padding resize for kspace and complex image
    // data is first fft to kspace and then zero padding then ifft to image domain
    // the scaling is handled to keep the noise variance
    // data: the 1st and 2rd dimensions are resized
    bool zpadResize2D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, hoNDArray<T>& dataResized);
    // data: [RO E1 SLC E2 CON PHS REP SET], 1, 2 and 4th dimensions are resized
    bool zpadResize3D(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataResized);

    // zero-padding resize with kspace as input
    bool zpadResize2DOnKSpace(const hoNDArray<T>& kspace, size_t sizeX, size_t sizeY, hoNDArray<T>& dataResized);
    bool zpadResize3DOnKSpace(const hoNDArray<T>& kspace, size_t sizeX, size_t sizeY, size_t sizeZ, hoNDArray<T>& dataResized);

    // zero padding resize with filter
    // data is first fft to kspace, then zero padding, then filtered and ifft to image domain
    // filter2D: 2D array for kspace filter
    // data: the 1st and 2rd dimensions are resized
    bool zpadResize2DFilter(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, const hoNDArray<T>& filter2D, hoNDArray<T>& dataResized);
    // filter3D: 3D array for kspace filter
    // data: [RO E1 SLC E2 CON PHS REP SET], 1, 2 and 4th dimensions are resized
    bool zpadResize3DFilter(const hoNDArray<T>& data, size_t sizeX, size_t sizeY, size_t sizeZ, const hoNDArray<T>& filter3D, hoNDArray<T>& dataResized);

    // ------------------------------------------------------------------------
    // kspace filter in image domain
    // ------------------------------------------------------------------------
    // in image domain for ISMRMRD dimension order
    bool kspacefilterROImage(hoNDArray<T>& data, const hoNDArray<T>& fRO);
    bool kspacefilterROImage(const hoNDArray<T>& data, const hoNDArray<T>& fRO, hoNDArray<T>& dataFiltered);
    bool kspacefilterE1Image(const hoNDArray<T>& data, const hoNDArray<T>& fE1, hoNDArray<T>& dataFiltered);
    bool kspacefilterE2Image(const hoNDArray<T>& data, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    bool kspacefilterE1E2Image(const hoNDArray<T>& data, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);
    bool kspacefilterROE1E2Image(const hoNDArray<T>& data, const hoNDArray<T>& fRO, const hoNDArray<T>& fE1, const hoNDArray<T>& fE2, hoNDArray<T>& dataFiltered);

    // ------------------------------------------------------------------------
    // coil sensitivity
    // ------------------------------------------------------------------------
    // coil estimation using NIH method
    // data: in image domain, at least 3D [RO E1 CHA], the coil map will be estimated for every 2D kspace
    bool coilMap2DNIH(const hoNDArray<T>& data, hoNDArray<T>& coilMap, ISMRMRDCOILMAPALGO algo, size_t ks=11, size_t power=3, size_t iterNum=5, typename realType<T>::Type thres=1e-3, bool useGPU=true);

    // data: in image domain, at least 4D [RO E1 E2 CHA], the coil map will be estimated for every 2D kspace [RO E1 CHA] across E2
    bool coilMap3DNIH(const hoNDArray<T>& data, hoNDArray<T>& coilMap, ISMRMRDCOILMAPALGO algo, size_t ks=7, size_t power=3, size_t iterNum=5, typename realType<T>::Type thres=1e-3, bool true3D=false);

    // the Souheil method
    // data: [RO E1 CHA], only 3D array
    // these functions are using 2D data correlation matrix
    bool coilMap2DNIHInner(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t power);

    // data: [RO E1 E2 CHA], this functions uses true 3D data correlation matrix
    bool coilMap3DNIHInner(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t power);

    // the Souheil iteration method
    // data: [RO E1 CHA], only 3D array
    bool coilMap2DNIH2Inner(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t iterNum, typename realType<T>::Type thres);

    // data: [RO E1 E2 CHA], true 3D coil map estimation
    bool coilMap3DNIH2Inner(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t kz, size_t iterNum, typename realType<T>::Type thres);

    // sum of square coil combination
    // data: in image domain, at least 3D [RO E1 CHA]
    bool sumOfSquare(const hoNDArray<T>& data, hoNDArray<T>& sos);

    // coil map weighted coil combination
    // data: in image domain, at least 3D [RO E1 CHA ...]
    // coilMap: [RO E1 CHA ... ]
    bool coilCombine(const hoNDArray<T>& data, const hoNDArray<T>& coilMap, hoNDArray<T>& combined);

    // data: in image domain, at least 4D [RO E1 E2 CHA ...]
    // coilMap: [RO E1 E2 CHA ... ]
    bool coilCombine3D(const hoNDArray<T>& data, const hoNDArray<T>& coilMap, hoNDArray<T>& combined);

    // ------------------------------------------------------------------------
    // kspace utility functions
    // ------------------------------------------------------------------------
    // get the conjugate symmetric kspace for 2D case
    // kspace : [RO E1 ...]
    // kspaceConj: [RO E1 ...]
    bool conjugateSymmetry2D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj);

    // kspace : [RO E1 E2 ...]
    // kspaceConj: [RO E1 E2 ...]
    bool conjugateSymmetry3D(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj);
};

}}

#include "gtPlusISMRMRDReconUtil.hxx"
