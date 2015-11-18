
/** \file   mri_core_spirit.cpp
    \brief  SPIRIT implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue

    References to the implementation can be found in:

    Lustig M, Pauly JM. 
    SPIRiT: Iterative self-consistent parallel imaging reconstruction from arbitrary k-space. 
    Magnetic Resonance in Medicine 2010;64(2):457-471.

    ISMRM 2013 sunrise course on Parallel Imaging
    Michael S. Hansen, Philip Beatty
    http://gadgetron.sourceforge.net/sunrise/
    http://cds.ismrm.org/protected/13MPresentations/abstracts/7059.pdf 
*/

#include "mri_core_spirit.h"
#include "mri_core_utility.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

template <typename T> 
void spirit2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                                    double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1,
                                    size_t startRO, size_t endRO, size_t startE1, size_t endE1,
                                    hoNDArray<T>& convKer)
{
    try
    {
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit2d_calib_convolution_kernel(start, end) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void spirit2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                                    double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray<T>& convKer)
{
    try
    {
        size_t startRO = 0;
        size_t endRO = acsSrc.get_size(0) - 1;
        size_t startE1 = 0;
        size_t endE1 = acsSrc.get_size(1) - 1;

        GADGET_CATCH_THROW( Gadgetron::spirit2d_calib_convolution_kernel(acsSrc, acsDst, thres, kRO, kE1, oRO, oE1, startRO, endRO, startE1, endE1, convKer) );
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit2d_calib_convolution_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void spirit2d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst,
                                    hoNDArray<unsigned short>& dataMask, double thres,
                                    size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray<T>& convKer)
{
    try
    {
        GADGET_CHECK_THROW(dataSrc.dimensions_equal(&dataMask));
        GADGET_CHECK_THROW(dataDst.dimensions_equal(&dataMask));

        // find the fully sampled region
        size_t RO = dataMask.get_size(0);
        size_t E1 = dataMask.get_size(1);
        size_t srcCHA = dataSrc.get_size(2);
        size_t dstCHA = dataDst.get_size(2);

        size_t startRO(0), endRO(0), startE1(0), endE1(0);

        size_t ro, e1, scha, dcha;

        for (e1 = 0; e1 < E1; e1++)
        {
            for (ro = 0; ro < RO; ro++)
            {
                if (dataMask(ro, e1)>0)
                {
                    if (ro < startRO) startRO = ro;
                    if (ro > endRO) endRO = ro;

                    if (e1 < startE1) startE1 = e1;
                    if (e1 > endE1) endE1 = e1;
                }
            }
        }

        GADGET_CHECK_THROW(endRO>startRO + kRO);
        GADGET_CHECK_THROW(endE1>startE1 + kE1);

        GADGET_CATCH_THROW(Gadgetron::spirit2d_calib_convolution_kernel(dataSrc, dataDst, thres, kRO, kE1, oRO, oE1, startRO, endRO, startE1, endE1, convKer));
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit2d_calib_convolution_kernel(dataMask) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& dataSrc, const hoNDArray< std::complex<float> >& dataDst, hoNDArray<unsigned short>& dataMask, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& dataSrc, const hoNDArray< std::complex<double> >& dataDst, hoNDArray<unsigned short>& dataMask, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void spirit2d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, hoNDArray<T>& kIm, bool minusI)
{
    try
    {
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit2d_image_domain_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit2d_image_domain_kernel(const hoNDArray< std::complex<float> >& convKer, size_t RO, size_t E1, hoNDArray< std::complex<float> >& kIm, bool minusI);
template EXPORTMRICORE void spirit2d_image_domain_kernel(const hoNDArray< std::complex<double> >& convKer, size_t RO, size_t E1, hoNDArray< std::complex<double> >& kIm, bool minusI);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                                    double thres, double overDetermineRatio,
                                    size_t kRO, size_t kE1, size_t kE2,
                                    size_t oRO, size_t oE1, size_t oE2,
                                    size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
                                    hoNDArray<T>& convKer)
{
    try
    {
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_calib_convolution_kernel(start, end) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                                    double thres, double overDetermineRatio,
                                    size_t kRO, size_t kE1, size_t kE2,
                                    size_t oRO, size_t oE1, size_t oE2,
                                    hoNDArray<T>& convKer)
{
    try
    {
        size_t startRO = 0;
        size_t endRO = acsSrc.get_size(0) - 1;
        size_t startE1 = 0;
        size_t endE1 = acsSrc.get_size(1) - 1;
        size_t startE2 = 0;
        size_t endE2 = acsSrc.get_size(2) - 1;

        GADGET_CATCH_THROW(Gadgetron::spirit3d_calib_convolution_kernel(acsSrc, acsDst, thres, overDetermineRatio, kRO, kE1, kE2, oRO, oE1, oE2, startRO, endRO, startE1, endE1, startE2, endE2, convKer));
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_calib_convolution_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask,
                                    double thres, double overDetermineRatio,
                                    size_t kRO, size_t kE1, size_t kE2,
                                    size_t oRO, size_t oE1, size_t oE2,
                                    hoNDArray<T>& convKer)
{
    try
    {
        GADGET_CHECK_THROW(dataSrc.dimensions_equal(&dataMask));
        GADGET_CHECK_THROW(dataDst.dimensions_equal(&dataMask));

        // find the fully sampled region
        size_t RO = dataMask.get_size(0);
        size_t E1 = dataMask.get_size(1);
        size_t E2 = dataMask.get_size(2);
        size_t srcCHA = dataSrc.get_size(3);
        size_t dstCHA = dataDst.get_size(3);

        size_t startRO(0), endRO(0), startE1(0), endE1(0), startE2(0), endE2(0);

        size_t ro, e1, e2, scha, dcha;

        for (e2 = 0; e2 < E2; e2++)
        {
            for (e1 = 0; e1 < E1; e1++)
            {
                for (ro = 0; ro < RO; ro++)
                {
                    if (dataMask(ro, e1, e2) > 0)
                    {
                        if (ro < startRO) startRO = ro;
                        if (ro > endRO) endRO = ro;

                        if (e1 < startE1) startE1 = e1;
                        if (e1 > endE1) endE1 = e1;

                        if (e2 < startE2) startE2 = e2;
                        if (e2 > endE2) endE2 = e2;
                    }
                }
            }
        }

        GADGET_CHECK_THROW(endRO>startRO + kRO);
        GADGET_CHECK_THROW(endE1>startE1 + kE1);
        GADGET_CHECK_THROW(endE2>startE2 + kE2);

        GADGET_CATCH_THROW(Gadgetron::spirit3d_calib_convolution_kernel(dataSrc, dataDst, thres, overDetermineRatio, kRO, kE1, kE2, oRO, oE1, oE2, startRO, endRO, startE1, endE1, startE2, endE2, convKer));
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_calib_convolution_kernel(dataMask) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, hoNDArray<unsigned short>& dataMask, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, hoNDArray<unsigned short>& dataMask, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, size_t E2, hoNDArray<T>& kIm, bool minusI)
{
    try
    {
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_image_domain_kernel(...) ... ");
    }

    return;
}

// ------------------------------------------------------------------------

}
