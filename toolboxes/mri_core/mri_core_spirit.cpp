
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
#include "mri_core_grappa.h"
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
                                    hoNDArray<T>& convKer, bool minusI)
{
    try
    {
        //hoNDArray<T> ker;
        //GADGET_CATCH_THROW(Gadgetron::spirit2d_calib(acsSrc, acsDst, thres, kRO, kE1, oRO, oE1, startRO, endRO, startE1, endE1, ker));
        //GADGET_CATCH_THROW(Gadgetron::spirit2d_convert_to_convolution_kernel(ker, kRO, kE1, oRO, oE1, convKer, minusI));

        size_t RO = acsSrc.get_size(0);
        size_t E1 = acsSrc.get_size(1);
        size_t srcCHA = acsSrc.get_size(2);
        size_t dstCHA = acsDst.get_size(2);

        hoNDArray<T> acsSrc3D, acsDst3D;
        acsSrc3D.create(RO, E1, 1, srcCHA, const_cast<T*>(acsSrc.begin()));
        acsDst3D.create(RO, E1, 1, dstCHA, const_cast<T*>(acsDst.begin()));

        hoNDArray<T> ker, convKer3D;
        GADGET_CATCH_THROW(Gadgetron::spirit3d_calib(acsSrc3D, acsDst3D, thres, 0, kRO, kE1, 1, oRO, oE1, 1, startRO, endRO, startE1, endE1, 0, 0, ker));
        GADGET_CATCH_THROW(Gadgetron::spirit3d_convert_to_convolution_kernel(ker, kRO, kE1, 1, oRO, oE1, 1, convKer3D, minusI));

        convKer.create(convKer3D.get_size(0), convKer3D.get_size(1), convKer3D.get_size(3), convKer3D.get_size(4));
        memcpy(convKer.begin(), convKer3D.begin(), convKer.get_number_of_bytes());
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit2d_calib_convolution_kernel(start, end) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<float> >& convKer, bool minusI);
template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<double> >& convKer, bool minusI);

// ------------------------------------------------------------------------

template <typename T> 
void spirit2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
    double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray<T>& convKer, bool minusI)
{
    try
    {
        size_t startRO = 0;
        size_t endRO = acsSrc.get_size(0) - 1;
        size_t startE1 = 0;
        size_t endE1 = acsSrc.get_size(1) - 1;

        GADGET_CATCH_THROW(Gadgetron::spirit2d_calib_convolution_kernel(acsSrc, acsDst, thres, kRO, kE1, oRO, oE1, startRO, endRO, startE1, endE1, convKer, minusI));
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit2d_calib_convolution_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray< std::complex<float> >& convKer, bool minusI);
template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray< std::complex<double> >& convKer, bool minusI);

// ------------------------------------------------------------------------

template <typename T> 
void spirit2d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst,
                                    hoNDArray<unsigned short>& dataMask, double thres,
                                    size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray<T>& convKer, bool minusI)
{
    try
    {
        GADGET_CHECK_THROW(dataSrc.dimensions_equal(dataMask));
        GADGET_CHECK_THROW(dataDst.dimensions_equal(dataMask));

        // find the fully sampled region
        size_t RO = dataMask.get_size(0);
        size_t E1 = dataMask.get_size(1);
        size_t srcCHA = dataSrc.get_size(2);
        size_t dstCHA = dataDst.get_size(2);

        size_t startRO(0), endRO(0), startE1(0), endE1(0);

        size_t ro, e1;

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

        GADGET_CATCH_THROW(Gadgetron::spirit2d_calib_convolution_kernel(dataSrc, dataDst, thres, kRO, kE1, oRO, oE1, startRO, endRO, startE1, endE1, convKer, minusI));
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit2d_calib_convolution_kernel(dataMask) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& dataSrc, const hoNDArray< std::complex<float> >& dataDst, hoNDArray<unsigned short>& dataMask, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray< std::complex<float> >& convKer, bool minusI);
template EXPORTMRICORE void spirit2d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& dataSrc, const hoNDArray< std::complex<double> >& dataDst, hoNDArray<unsigned short>& dataMask, double thres, size_t kRO, size_t kE1, size_t oRO, size_t oE1, hoNDArray< std::complex<double> >& convKer, bool minusI);

// ------------------------------------------------------------------------

template <typename T> 
void spirit2d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, hoNDArray<T>& kIm)
{
    try
    {
        Gadgetron::grappa2d_image_domain_kernel(convKer, RO, E1, kIm);
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit2d_image_domain_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit2d_image_domain_kernel(const hoNDArray< std::complex<float> >& convKer, size_t RO, size_t E1, hoNDArray< std::complex<float> >& kIm);
template EXPORTMRICORE void spirit2d_image_domain_kernel(const hoNDArray< std::complex<double> >& convKer, size_t RO, size_t E1, hoNDArray< std::complex<double> >& kIm);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_calib(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                double thres, double overDetermineRatio, 
                size_t kRO, size_t kE1, size_t kE2,
                size_t oRO, size_t oE1, size_t oE2,
                size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
                hoNDArray<T>& ker)
{
    try
    {
        GADGET_CHECK_THROW(acsSrc.get_size(0) == acsDst.get_size(0));
        GADGET_CHECK_THROW(acsSrc.get_size(1) == acsDst.get_size(1));
        GADGET_CHECK_THROW(acsSrc.get_size(2) == acsDst.get_size(2));
        GADGET_CHECK_THROW(acsSrc.get_size(3) >= acsDst.get_size(3));

        size_t RO = acsSrc.get_size(0);
        size_t E1 = acsSrc.get_size(1);
        size_t E2 = acsSrc.get_size(2);
        size_t srcCHA = acsSrc.get_size(3);
        size_t dstCHA = acsDst.get_size(3);

        long long kROhalf = kRO / 2;
        if (2 * kROhalf == kRO)
        {
            GWARN_STREAM("spirit3d_calib - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2 * kROhalf + 1;

        long long kE1half = kE1 / 2;
        if (2 * kE1half == kE1)
        {
            GWARN_STREAM("spirit3d_calib - 2*kE1half == kE1 " << kE1);
        }
        kE1 = 2 * kE1half + 1;

        long long kE2half = kE2 / 2;
        if (2 * kE2half == kE2)
        {
            GWARN_STREAM("spirit3d_calib - 2*kE2half == kE2 " << kE2);
        }
        kE2 = 2 * kE2half + 1;

        if (oRO > kRO) oRO = kRO;
        if (oE1 > kE1) oE1 = kE1;
        if (oE2 > kE2) oE2 = kE2;

        long long oROhalf = oRO / 2;
        if (2 * oROhalf == oRO)
        {
            GWARN_STREAM("spirit3d_calib - 2*oROhalf == oRO " << oRO);
        }
        oRO = 2 * oROhalf + 1;

        long long oE1half = oE1 / 2;
        if (2 * oE1half == oE1)
        {
            GWARN_STREAM("spirit3d_calib - 2*oE1half == oE1 " << oE1);
        }
        oE1 = 2 * oE1half + 1;

        long long oE2half = oE2 / 2;
        if (2 * oE2half == oE2)
        {
            GWARN_STREAM("spirit3d_calib - 2*oE2half == oE2 " << oE2);
        }
        oE2 = 2 * oE2half + 1;

        // allocate kernel
        ker.create(kRO, kE1, kE2, srcCHA, dstCHA, oRO, oE1, oE2);

        // loop over the calibration region and assemble the equation
        // Ax = b

        size_t sRO = kROhalf;
        size_t eRO = RO - kROhalf - 1;
        size_t lenRO = eRO - sRO + 1;

        size_t sE1 = kE1half;
        size_t eE1 = E1 - kE1half - 1;
        size_t lenE1 = eE1 - sE1 + 1;

        size_t sE2 = kE2half;
        size_t eE2 = E2 - kE2half - 1;
        size_t lenE2 = eE2 - sE2 + 1;

        size_t colA = (kRO*kE1*kE2 - 1)*srcCHA;
        if (overDetermineRatio > 1.0)
        {
            size_t maxRowA = (size_t)std::ceil(overDetermineRatio*colA);
            size_t maxROUsed = maxRowA / (lenE1*lenE2);
            if (maxROUsed < lenRO)
            {
                // find the peak signal of acsSrc
                hoNDArray<T> acsSrc1stCha(RO, E1, E2, const_cast<T*>(acsSrc.begin()));
                hoNDArray<T> acsSrc1stChaSumE2(RO, E1, 1), acsSrc1stChaSumE2E1(RO, 1, 1);

                try
                {
                    Gadgetron::sum_over_dimension(acsSrc1stCha, acsSrc1stChaSumE2, 2);
                    Gadgetron::sum_over_dimension(acsSrc1stChaSumE2, acsSrc1stChaSumE2E1, 1);

                    T maxSignal;
                    size_t roInd(0);
                    try
                    {
                        Gadgetron::maxAbsolute(acsSrc1stChaSumE2E1, maxSignal, roInd);

                        if (roInd > maxROUsed / 2 + kROhalf)
                        {
                            sRO = roInd - maxROUsed / 2;
                        }
                        else
                        {
                            sRO = kROhalf;
                        }

                        if (sRO + maxROUsed - 1 <= RO - kROhalf - 1)
                        {
                            eRO = sRO + maxROUsed - 1;
                        }
                        else
                        {
                            eRO = RO - kROhalf - 1;
                        }

                        lenRO = eRO - sRO + 1;
                        GDEBUG_STREAM("spirit3d_calib - overDetermineRatio = " << overDetermineRatio << " ; RO data range used : [" << sRO << " " << eRO << "] ...");
                    }
                    catch (...)
                    {
                        GWARN_STREAM("spirit3d_calib - overDetermineRatio is ignored ... ");
                        throw;
                    }
                }
                catch (...)
                {
                    GWARN_STREAM("spirit3d_calib - overDetermineRatio is ignored ... ");
                    throw;
                }
            }
        }

        size_t rowA = lenRO*lenE1*lenE2;
        size_t colB = dstCHA;

#pragma omp parallel default(none) shared(sRO, eRO, sE1, eE1, sE2, eE2, oRO, oE1, oE2, lenRO, lenE1, lenE2, rowA, colA, colB, kROhalf, kE1half, kE2half, oROhalf, oE1half, oE2half, acsSrc, acsDst, srcCHA, dstCHA, thres, ker, std::cerr) num_threads( (int)(oRO*oE1*oE2) ) if (oRO*oE1*oE2>=3 && oRO*oE1*oE2<9)
        {
            hoMatrix<T> A(rowA, colA);
            hoMatrix<T> B(rowA, colB);
            hoMatrix<T> x(A.cols(), B.cols());

            T* pA = A.begin();
            T* pB = B.begin();

            long long kInd = 0;
#pragma omp for
            for (kInd = 0; kInd<(long long)(oRO*oE1*oE2); kInd++)
            {
                long long oe2 = kInd / (oRO*oE1);
                long long oe1 = kInd - oe2*oRO*oE1;
                oe1 /= oRO;
                long long oro = kInd - oe2*oRO*oE1 - oe1*oRO;

                oe2 -= oE2half;
                oe1 -= oE1half;
                oro -= oROhalf;

                long long dRO, dE1, dE2;

                for (long long e2 = (long long)sE2; e2 <= (long long)eE2; e2++)
                {
                    dE2 = e2 + oe2;

                    for (long long e1 = (long long)sE1; e1 <= (long long)eE1; e1++)
                    {
                        dE1 = e1 + oe1;

                        for (long long ro = sRO; ro <= (long long)eRO; ro++)
                        {
                            dRO = ro + oro;

                            long long rInd = (e2 - sE2)*lenRO*lenE1 + (e1 - sE1)*lenRO + ro - sRO;

                            // fill matrix A
                            size_t col = 0;
                            for (size_t src = 0; src<srcCHA; src++)
                            {
                                for (long long ke2 = -kE2half; ke2 <= kE2half; ke2++)
                                {
                                    for (long long ke1 = -kE1half; ke1 <= kE1half; ke1++)
                                    {
                                        for (long long kro = -kROhalf; kro <= kROhalf; kro++)
                                        {
                                            if (kro != oro || ke1 != oe1 || ke2 != oe2)
                                            {
                                                //A(rInd, col++) = acsSrc(ro+kro, e1+ke1, e2+ke2, src);
                                                pA[rInd + col*rowA] = acsSrc(ro + kro, e1 + ke1, e2 + ke2, src);
                                                col++;
                                            }
                                        }
                                    }
                                }
                            }

                            // fill matrix B
                            for (size_t dst = 0; dst<dstCHA; dst++)
                            {
                                //B(rInd, dst) = acsDst(dRO, dE1, dE2, dst);
                                pB[rInd + dst*rowA] = acsDst(dRO, dE1, dE2, dst);
                            }
                        }
                    }
                }

                SolveLinearSystem_Tikhonov(A, B, x, thres);

                long long ind(0);

                std::vector<size_t> kerInd(8);
                kerInd[7] = oe2 + oE2half;
                kerInd[6] = oe1 + oE1half;
                kerInd[5] = oro + oROhalf;

                for (size_t src = 0; src<srcCHA; src++)
                {
                    kerInd[3] = src;
                    for (long long ke2 = -kE2half; ke2 <= kE2half; ke2++)
                    {
                        kerInd[2] = ke2 + kE2half;
                        for (long long ke1 = -kE1half; ke1 <= kE1half; ke1++)
                        {
                            kerInd[1] = ke1 + kE1half;
                            for (long long kro = -kROhalf; kro <= kROhalf; kro++)
                            {
                                kerInd[0] = kro + kROhalf;

                                if (kro != 0 || ke1 != 0 || ke2 != 0)
                                {
                                    for (size_t dst = 0; dst<dstCHA; dst++)
                                    {
                                        kerInd[4] = dst;
                                        size_t offset = ker.calculate_offset(kerInd);
                                        ker(offset) = x(ind, dst);
                                    }
                                    ind++;
                                }
                                else
                                {
                                    for (size_t dst = 0; dst<dstCHA; dst++)
                                    {
                                        kerInd[4] = dst;
                                        size_t offset = ker.calculate_offset(kerInd);
                                        ker(offset) = 0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for(size_t kk=0; kk>ker.get_number_of_elements(); kk++)
        {
            if(std::isnan(ker(kk).real()) || std::isnan(ker(kk).imag()))
            {
                GERROR_STREAM("nan detected in spirit3d_calib ker ... ");
                throw;
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_calib(...) ... ");
    }

    return;
}


template EXPORTMRICORE void spirit3d_calib(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, hoNDArray< std::complex<float> >& ker);
template EXPORTMRICORE void spirit3d_calib(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, hoNDArray< std::complex<double> >& ker);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_convert_to_convolution_kernel(const hoNDArray<T>& ker, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray<T>& convKer, bool minusI)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(3));
        long long dstCHA = (long long)(ker.get_size(4));

        GADGET_CHECK_THROW(kRO == ker.get_size(0));
        GADGET_CHECK_THROW(kE1 == ker.get_size(1));
        GADGET_CHECK_THROW(kE2 == ker.get_size(2));
        GADGET_CHECK_THROW(oRO == ker.get_size(5));
        GADGET_CHECK_THROW(oE1 == ker.get_size(6));
        GADGET_CHECK_THROW(oE2 == ker.get_size(7));

        long long kROhalf = kRO / 2;
        long long kE1half = kE1 / 2;
        long long kE2half = kE2 / 2;
        long long oROhalf = oRO / 2;
        long long oE1half = oE1 / 2;
        long long oE2half = oE2 / 2;

        /// fill the convolution kernels
        long long convKRO = 2 * kRO - 1;
        long long convKE1 = 2 * kE1 - 1;
        long long convKE2 = 2 * kE2 - 1;

        /// fill in convolution kernel
        hoNDArray<T> convKerBuf;
        convKerBuf.create(convKRO, convKE1, convKE2, srcCHA, dstCHA, oRO, oE1, oE2);
        Gadgetron::clear(&convKerBuf);

        long long oro, oe1, oe2, kro, ke1, ke2, src, dst;
        std::vector<size_t> kerInd(8), convKerInd(8);
        for (oe2 = -oE2half; oe2 <= oE2half; oe2++)
        {
            kerInd[7] = oe2 + oE2half;
            convKerInd[7] = oe2 + oE2half;

            for (oe1 = -oE1half; oe1 <= oE1half; oe1++)
            {
                kerInd[6] = oe1 + oE1half;
                convKerInd[6] = oe1 + oE1half;

                for (oro = -oROhalf; oro <= oROhalf; oro++)
                {
                    kerInd[5] = oro + oROhalf;
                    convKerInd[5] = oro + oROhalf;

                    for (ke2 = -kE2half; ke2 <= kE2half; ke2++)
                    {
                        kerInd[2] = ke2 + kE2half;

                        for (ke1 = -kE1half; ke1 <= kE1half; ke1++)
                        {
                            kerInd[1] = ke1 + kE1half;

                            for (kro = -kROhalf; kro <= kROhalf; kro++)
                            {
                                long long iro = kro - oro + kRO - 1;
                                long long ie1 = ke1 - oe1 + kE1 - 1;
                                long long ie2 = ke2 - oe2 + kE2 - 1;

                                kerInd[0] = kro + kROhalf;

                                convKerInd[0] = iro;
                                convKerInd[1] = ie1;
                                convKerInd[2] = ie2;

                                for (dst = 0; dst<dstCHA; dst++)
                                {
                                    kerInd[4] = dst;
                                    convKerInd[4] = dst;

                                    for (src = 0; src<srcCHA; src++)
                                    {
                                        kerInd[3] = src;
                                        convKerInd[3] = src;

                                        size_t offsetKer = ker.calculate_offset(kerInd);
                                        size_t offsetConvKer = convKerBuf.calculate_offset(convKerInd);

                                        convKerBuf(offsetConvKer) = ker(offsetKer);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        hoNDArray<T> convKer2, convKer3;
        hoNDArray<T> convKernMean(convKRO, convKE1, convKE2, srcCHA, dstCHA, 1, 1, 1);
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(convKerBuf, convKer2, 7));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(convKer2, convKer3, 6));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(convKer3, convKernMean, 5));

        GADGET_CATCH_THROW(Gadgetron::scal((typename realType<T>::Type)(1.0 / (oRO*oE1*oE2)), convKernMean));

        // flip the kernel
        convKer.create(convKRO, convKE1, convKE2, srcCHA, dstCHA);
        Gadgetron::clear(&convKer);

        for (ke2 = 0; ke2<convKE2; ke2++)
        {
            for (ke1 = 0; ke1<convKE1; ke1++)
            {
                for (kro = 0; kro<convKRO; kro++)
                {
                    for (dst = 0; dst<dstCHA; dst++)
                    {
                        for (src = 0; src<srcCHA; src++)
                        {
                            T value = convKernMean(convKRO - 1 - kro, convKE1 - 1 - ke1, convKE2 - 1 - ke2, src, dst);
                            convKer(kro, ke1, ke2, src, dst) = value;
                        }
                    }
                }
            }
        }

        if (minusI)
        {
            for (dst = 0; dst<dstCHA; dst++)
            {
                T value = convKer(kRO - 1, kE1 - 1, kE2 - 1, dst, dst);
                convKer(kRO - 1, kE1 - 1, kE2 - 1, dst, dst) = value - T(1.0);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_convert_to_convolution_kernel(...) ... ");
    }

    return;
}

template <typename T> void spirit3d_convert_to_convolution_kernel(const hoNDArray< std::complex<float> >& ker, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<float> >& convKer, bool minusI);
template <typename T> void spirit3d_convert_to_convolution_kernel(const hoNDArray< std::complex<double> >& ker, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<double> >& convKer, bool minusI);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                                    double thres, double overDetermineRatio,
                                    size_t kRO, size_t kE1, size_t kE2,
                                    size_t oRO, size_t oE1, size_t oE2,
                                    size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
                                    hoNDArray<T>& convKer, bool minusI)
{
    try
    {
        hoNDArray<T> ker;
        GADGET_CATCH_THROW(Gadgetron::spirit3d_calib(acsSrc, acsDst, thres, overDetermineRatio, kRO, kE1, kE2, oRO, oE1, oE2, startRO, endRO, startE1, endE1, startE2, endE2, ker));
        GADGET_CATCH_THROW(Gadgetron::spirit3d_convert_to_convolution_kernel(ker, kRO, kE1, kE2, oRO, oE1, oE2, convKer, minusI));
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_calib_convolution_kernel(start, end) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, hoNDArray< std::complex<float> >& convKer, bool minusI);
template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, hoNDArray< std::complex<double> >& convKer, bool minusI);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                                    double thres, double overDetermineRatio,
                                    size_t kRO, size_t kE1, size_t kE2,
                                    size_t oRO, size_t oE1, size_t oE2,
                                    hoNDArray<T>& convKer, bool minusI)
{
    try
    {
        size_t startRO = 0;
        size_t endRO = acsSrc.get_size(0) - 1;
        size_t startE1 = 0;
        size_t endE1 = acsSrc.get_size(1) - 1;
        size_t startE2 = 0;
        size_t endE2 = acsSrc.get_size(2) - 1;

        GADGET_CATCH_THROW(Gadgetron::spirit3d_calib_convolution_kernel(acsSrc, acsDst, thres, overDetermineRatio, kRO, kE1, kE2, oRO, oE1, oE2, startRO, endRO, startE1, endE1, startE2, endE2, convKer, minusI));
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_calib_convolution_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<float> >& convKer, bool minusI);
template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<double> >& convKer, bool minusI);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask,
                                    double thres, double overDetermineRatio,
                                    size_t kRO, size_t kE1, size_t kE2,
                                    size_t oRO, size_t oE1, size_t oE2,
                                    hoNDArray<T>& convKer, bool minusI)
{
    try
    {
        GADGET_CHECK_THROW(dataSrc.dimensions_equal(dataMask));
        GADGET_CHECK_THROW(dataDst.dimensions_equal(dataMask));

        // find the fully sampled region
        size_t RO = dataMask.get_size(0);
        size_t E1 = dataMask.get_size(1);
        size_t E2 = dataMask.get_size(2);
        size_t srcCHA = dataSrc.get_size(3);
        size_t dstCHA = dataDst.get_size(3);

        size_t startRO(0), endRO(0), startE1(0), endE1(0), startE2(0), endE2(0);

        size_t ro, e1, e2;

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

        GADGET_CATCH_THROW(Gadgetron::spirit3d_calib_convolution_kernel(dataSrc, dataDst, thres, overDetermineRatio, kRO, kE1, kE2, oRO, oE1, oE2, startRO, endRO, startE1, endE1, startE2, endE2, convKer, minusI));
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_calib_convolution_kernel(dataMask) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, hoNDArray<unsigned short>& dataMask, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<float> >& convKer, bool minusI);
template EXPORTMRICORE void spirit3d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, hoNDArray<unsigned short>& dataMask, double thres, double overDetermineRatio, size_t kRO, size_t kE1, size_t kE2, size_t oRO, size_t oE1, size_t oE2, hoNDArray< std::complex<double> >& convKer, bool minusI);

// ------------------------------------------------------------------------

template <typename T> 
void spirit3d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, size_t E2, hoNDArray<T>& kIm, bool preset_kIm_with_zeros)
{
    try
    {
        Gadgetron::grappa3d_image_domain_kernel(convKer, RO, E1, E2, kIm, preset_kIm_with_zeros);
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_image_domain_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit3d_image_domain_kernel(const hoNDArray< std::complex<float> >& convKer, size_t RO, size_t E1, size_t E2, hoNDArray< std::complex<float> >& kIm, bool preset_kIm_with_zeros);
template EXPORTMRICORE void spirit3d_image_domain_kernel(const hoNDArray< std::complex<double> >& convKer, size_t RO, size_t E1, size_t E2, hoNDArray< std::complex<double> >& kIm, bool preset_kIm_with_zeros);

// ------------------------------------------------------------------------

template <typename T>
void spirit3d_kspace_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, hoNDArray<T>& kIm)
{
    try
    {
        size_t kRO = convKer.get_size(0);
        size_t kE1 = convKer.get_size(1);
        size_t kE2 = convKer.get_size(2);

        size_t srcCHA = convKer.get_size(3);
        size_t dstCHA = convKer.get_size(4);

        if (kIm.get_size(0) != kE1 || kIm.get_size(1) != kE2 || kIm.get_size(2) != srcCHA || kIm.get_size(3) != dstCHA || kIm.get_size(4) != RO)
        {
            kIm.create(kE1, kE2, srcCHA, dstCHA, RO);
        }

        hoNDArray<T> convKerScaled;
        convKerScaled = convKer;

        // since only RO is converted to image domain
        Gadgetron::scal((typename realType<T>::Type)(std::sqrt((double)(RO))), convKerScaled);

        // pad the kernel and go to image domain
        hoNDArray<T> kImRO(RO, kE1, kE2, srcCHA, dstCHA);
        Gadgetron::clear(kImRO);
        Gadgetron::pad(RO, kE1, kE2, convKer, kImRO, false);

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(kImRO);

        // go to the required dimension order
        std::vector<size_t> dim_order(5);
        dim_order[0] = 1;
        dim_order[1] = 2;
        dim_order[2] = 3;
        dim_order[3] = 4;
        dim_order[4] = 0;

        Gadgetron::permute(kImRO, kIm, dim_order);
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_image_domain_kernel_E1E2RO(...) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit3d_kspace_image_domain_kernel(const hoNDArray< std::complex<float> >& convKer, size_t RO, hoNDArray< std::complex<float> >& kIm);
template EXPORTMRICORE void spirit3d_kspace_image_domain_kernel(const hoNDArray< std::complex<double> >& convKer, size_t RO, hoNDArray< std::complex<double> >& kIm);

// ------------------------------------------------------------------------

template<typename T> 
void spirit3d_image_domain_kernel(const hoNDArray<T>& kImRO, size_t E1, size_t E2, hoNDArray<T>& kIm)
{
    try
    {
        std::vector<size_t> dim;
        kImRO.get_dimensions(dim);

        std::vector<size_t> dimR(dim);
        dimR[0] = E1;
        dimR[1] = E2;

        kIm.create(dimR);
        Gadgetron::clear(kIm);

        hoNDArray<T> kImROScaled(kImRO);
        Gadgetron::scal((typename realType<T>::Type)(std::sqrt((double)(E1*E2))), kImROScaled);
        Gadgetron::pad(E1, E2, dimR[2], kImROScaled, kIm, false);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kIm);
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit3d_image_domain_kernel(const hoNDArray<T>& kImRO, size_t E1, size_t E2, hoNDArray<T>& kIm) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit3d_image_domain_kernel(const hoNDArray< std::complex<float> >& kImRO, size_t E1, size_t E2, hoNDArray< std::complex<float> >& kIm);
template EXPORTMRICORE void spirit3d_image_domain_kernel(const hoNDArray< std::complex<double> >& kImRO, size_t E1, size_t E2, hoNDArray< std::complex<double> >& kIm);

// ------------------------------------------------------------------------

template <typename T> 
void spirit_image_domain_adjoint_kernel(const hoNDArray<T>& kIm, hoNDArray<T>& adjkIm)
{
    try
    {
        std::vector<size_t> dim, dimAdj, dimOrder;
        kIm.get_dimensions(dim);

        size_t NDim = dim.size();

        dimAdj = dim;
        dimAdj[NDim - 1] = dim[NDim - 2];
        dimAdj[NDim - 2] = dim[NDim - 1];

        if (!adjkIm.dimensions_equal(dimAdj))
        {
            adjkIm.create(dimAdj);
        }

        dimOrder.resize(NDim);

        for (size_t d = 0; d < NDim; d++)
        {
            dimOrder[d] = d;
        }
        dimOrder[NDim - 2] = NDim - 1;
        dimOrder[NDim - 1] = NDim - 2;

        Gadgetron::permute(kIm, adjkIm, dimOrder);
        Gadgetron::conjugate(adjkIm, adjkIm);
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit_image_domain_adjoint_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void spirit_image_domain_adjoint_kernel(const hoNDArray< std::complex<float> >& kIm, hoNDArray< std::complex<float> >& adjkIm);
template EXPORTMRICORE void spirit_image_domain_adjoint_kernel(const hoNDArray< std::complex<double> >& kIm, hoNDArray< std::complex<double> >& adjkIm);

// ------------------------------------------------------------------------

template <typename T> 
void spirit_adjoint_forward_kernel(const hoNDArray<T>& kImS2D, const hoNDArray<T>& kImD2S, hoNDArray<T>& kIm)
{
    try
    {
        std::vector<size_t> dimS2D;
        kImS2D.get_dimensions(dimS2D);

        size_t NDim = kImS2D.get_number_of_dimensions();

        long long srcCHA = dimS2D[NDim - 2];
        long long dstCHA = dimS2D[NDim - 1];

        GADGET_CHECK_THROW(kImD2S.get_number_of_dimensions() == NDim);
        GADGET_CHECK_THROW(kImD2S.get_number_of_elements() == kImS2D.get_number_of_elements());

        std::vector<size_t> dimRes(dimS2D);
        dimRes[NDim - 2] = dstCHA;

        kIm.create(dimRes);
        Gadgetron::clear(&kIm);

        size_t N = kImS2D.get_number_of_elements() / srcCHA / dstCHA;

        std::vector<size_t> dim(1);
        dim[0] = N;

        long long d;
    #pragma omp parallel default(none) private(d) shared(N, dim, dstCHA, srcCHA, kIm, kImS2D, kImD2S) num_threads( (int)dstCHA ) if (dstCHA > 4)
        {
            hoNDArray<T> ker(N);
            hoNDArray<T> dKer, kerS2D, kerD2S;

    #pragma omp for
            for (d = 0; d<dstCHA; d++)
            {
                for (long long dprime = 0; dprime<dstCHA; dprime++)
                {
                    dKer.create(dim, kIm.begin() + d*N + dprime*N*dstCHA);

                    for (long long s = 0; s<srcCHA; s++)
                    {
                        kerS2D.create(dim, const_cast<T*>(kImS2D.begin()) + s*N + dprime*N*srcCHA);
                        kerD2S.create(dim, const_cast<T*>(kImD2S.begin()) + d*N + s*N*dstCHA);

                        Gadgetron::multiply(kerS2D, kerD2S, ker);
                        Gadgetron::add(dKer, ker, dKer);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in spirit_adjoint_forward_kernel(...) ... ");
    }
}

template EXPORTMRICORE void spirit_adjoint_forward_kernel(const hoNDArray< std::complex<float> >& kImS2D, const hoNDArray< std::complex<float> >& kImD2S, hoNDArray< std::complex<float> >& kIm);
template EXPORTMRICORE void spirit_adjoint_forward_kernel(const hoNDArray< std::complex<double> >& kImS2D, const hoNDArray< std::complex<double> >& kImD2S, hoNDArray< std::complex<double> >& kIm);

// ------------------------------------------------------------------------

}
