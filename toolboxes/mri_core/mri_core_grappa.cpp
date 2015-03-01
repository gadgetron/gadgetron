
/** \file   mri_core_grappa.cpp
    \brief  GRAPPA implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue

    References to the implementation can be found in:

    Griswold MA, Jakob PM, Heidemann RM, Nittka M, Jellus V, Wang J, Kiefer B, Haase A. 
    Generalized autocalibrating partially parallel acquisitions (GRAPPA). 
    Magnetic Resonance in Medicine 2002;47(6):1202-1210.

    Kellman P, Epstein FH, McVeigh ER. 
    Adaptive sensitivity encoding incorporating temporal filtering (TSENSE). 
    Magnetic Resonance in Medicine 2001;45(5):846-852.

    Breuer FA, Kellman P, Griswold MA, Jakob PM. .
    Dynamic autocalibrated parallel imaging using temporal GRAPPA (TGRAPPA). 
    Magnetic Resonance in Medicine 2005;53(4):981-985.

    Saybasili H., Kellman P., Griswold MA., Derbyshire JA. Guttman, MA. 
    HTGRAPPA: Real-time B1-weighted image domain TGRAPPA reconstruction. 
    Magnetic Resonance in Medicine 2009;61(6): 1425-1433. 
*/

#include "mri_core_grappa.h"
#include "ismrmrd/ismrmrd.h"
#include "GadgetronTimer.h"
#include "mri_core_utility.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

void grappa2d_kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, size_t accelFactor, size_t kNE1, bool fitItself)
{
    if ( accelFactor == 1 )
    {
        kE1.resize(1, 0);
        oE1.resize(1, 0);
        return;
    }

    kE1.resize(kNE1, 0);
    if ( kNE1%2 == 0 )
    {
        long long k;
        for ( k=-((long long)kNE1/2-1); k<=(long long)kNE1/2; k++ )
        {
            kE1[k+kNE1/2-1] = (int)(k*accelFactor);
        }
    }
    else
    {
        long long k;
        for ( k=-(long long)kNE1/2; k<=(long long)kNE1/2; k++ )
        {
            kE1[k+kNE1/2] = (int)(k*accelFactor);
        }
    }

    if ( fitItself )
    {
        oE1.resize(accelFactor);
        for ( long long a=0; a<(long long)accelFactor; a++ )
        {
            oE1[a] = (int)a;
        }
    }
    else
    {
        oE1.resize(accelFactor-1);
        for ( long long a=1; a<(long long)accelFactor; a++ )
        {
            oE1[a-1] = (int)a;
        }
    }

    return;
}

template <typename T> 
void grappa2d_calib(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, double thres, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray<T>& ker)
{
    try
    {
        GADGET_CHECK_THROW(acsSrc.get_size(0)==acsDst.get_size(0));
        GADGET_CHECK_THROW(acsSrc.get_size(1)==acsDst.get_size(1));
        GADGET_CHECK_THROW(acsSrc.get_size(2)>=acsDst.get_size(2));

        size_t RO = acsSrc.get_size(0);
        size_t E1 = acsSrc.get_size(1);
        size_t srcCHA = acsSrc.get_size(2);
        size_t dstCHA = acsDst.get_size(2);

        const T* pSrc = acsSrc.begin();
        const T* pDst = acsDst.begin();

        long long kROhalf = kRO/2;
        if ( 2*kROhalf == kRO )
        {
            GWARN_STREAM("grappa<T>::calib(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2*kROhalf + 1;

        size_t kNE1 = kE1.size();
        size_t oNE1 = oE1.size();

        /// allocate kernel
        GADGET_CHECK_THROW(ker.createArray(kRO, kNE1, srcCHA, dstCHA, oNE1));

        /// loop over the calibration region and assemble the equation
        /// Ax = b

        size_t eRO = RO - kROhalf -1;
        size_t sE1 = std::abs(kE1[0]);
        size_t eE1 = E1 -1 - kE1[kNE1-1];

        size_t lenRO = eRO-kROhalf+1;

        size_t rowA = (eE1-sE1+1)*lenRO;
        size_t colA = kRO*kNE1*srcCHA;
        size_t colB = dstCHA*oNE1;

        hoMatrix<T> A;
        hoMatrix<T> B;
        hoMatrix<T> x( colA, colB );

        hoNDArray<T> A_mem(rowA, colA);
        A.createMatrix( rowA, colA, A_mem.begin() );
        T* pA = A.begin();

        hoNDArray<T> B_mem(rowA, colB);
        B.createMatrix( A.rows(), colB, B_mem.begin() );
        T* pB = B.begin();

        long long e1;
        for ( e1=(long long)sE1; e1<=(long long)eE1; e1++ )
        {
            for ( long long ro=kROhalf; ro<=(long long)eRO; ro++ )
            {
                long long rInd = (e1-sE1)*lenRO+ro-kROhalf;

                size_t src, dst, ke1, oe1;
                long long kro;

                /// fill matrix A
                size_t col = 0;
                size_t offset = 0;
                for ( src=0; src<srcCHA; src++ )
                {
                    for ( ke1=0; ke1<kNE1; ke1++ )
                    {
                        offset = src*RO*E1 + (e1+kE1[ke1])*RO;
                        for ( kro=-kROhalf; kro<=kROhalf; kro++ )
                        {
                            /// A(rInd, col++) = acsSrc(ro+kro, e1+kE1[ke1], src);
                            pA[rInd + col*rowA] = pSrc[ro+kro+offset];
                            col++;
                        }
                    }
                }

                /// fill matrix B
                col = 0;
                for ( oe1=0; oe1<oNE1; oe1++ )
                {
                    for ( dst=0; dst<dstCHA; dst++ )
                    {
                        B(rInd, col++) = acsDst(ro, e1+oE1[oe1], dst);
                    }
                }
            }
        }

        SolveLinearSystem_Tikhonov(A, B, x, thres);
        memcpy(ker.begin(), x.begin(), ker.get_number_of_bytes());
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa2d_calib(...) ... ");
    }

    return;
}

template EXPORTMRICORE void grappa2d_calib(const ho3DArray< std::complex<float> >& acsSrc, const ho3DArray< std::complex<float> >& acsDst, double thres, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray< std::complex<float> >& ker);
template EXPORTMRICORE void grappa2d_calib(const ho3DArray< std::complex<double> >& acsSrc, const ho3DArray< std::complex<double> >& acsDst, double thres, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho5DArray< std::complex<double> >& ker);

// ------------------------------------------------------------------------

template <typename T>
void grappa2d_convert_to_convolution_kernel(const ho5DArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho4DArray<T>& convKer)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(2));
        long long dstCHA = (long long)(ker.get_size(3));
        long long kNE1 = (long long)(kE1.size());
        long long oNE1 = (long long)(oE1.size());

        long long kROhalf = kRO / 2;
        if (2 * kROhalf == kRO)
        {
            GWARN_STREAM("grappa2d_convert_to_convolution_kernel - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2 * kROhalf + 1;

        //// fill the convolution kernels
        long long convKRO = 2 * kRO + 3;

        long long maxKE1 = std::abs(kE1[0]);
        if (std::abs(kE1[kNE1 - 1]) > maxKE1)
        {
            maxKE1 = std::abs(kE1[kNE1 - 1]);
        }
        long long convKE1 = 2 * maxKE1 + 1;

        //// allocate the convolution kernel
        convKer.createArray(convKRO, convKE1, srcCHA, dstCHA);
        Gadgetron::clear(&convKer);

        //// index
        long long oe1, kro, ke1, src, dst;

        //// fill the convolution kernel and sum up multiple kernels
        for (oe1 = 0; oe1<oNE1; oe1++)
        {
            for (ke1 = 0; ke1<kNE1; ke1++)
            {
                for (kro = -kROhalf; kro <= kROhalf; kro++)
                {
                    for (dst = 0; dst<dstCHA; dst++)
                    {
                        for (src = 0; src<srcCHA; src++)
                        {
                            convKer(-kro + kRO + 1, oE1[oe1] - kE1[ke1] + maxKE1, src, dst) = ker(kro + kROhalf, ke1, src, dst, oe1);
                        }
                    }

                }
            }
        }

        if (oE1[0] != 0)
        {
            for (dst = 0; dst<dstCHA; dst++)
            {
                convKer(kRO + 1, maxKE1, dst, dst) = 1.0;
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_convert_to_convolution_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void grappa2d_convert_to_convolution_kernel(const ho5DArray< std::complex<float> >& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho4DArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa2d_convert_to_convolution_kernel(const ho5DArray< std::complex<double> >& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, ho4DArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T>
void grappa2d_calib_convolution_kernel(const ho3DArray<T>& acsSrc, const ho3DArray<T>& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, ho4DArray<T>& convKer)
{
    try
    {
        std::vector<int> kE1, oE1;

        bool fitItself = false;
        if (&acsSrc != &acsDst) fitItself = true;

        grappa2d_kerPattern(kE1, oE1, accelFactor, kNE1, fitItself);

        ho5DArray<T> ker;
        grappa2d_calib(acsSrc, acsDst, thres, kRO, kE1, oE1, ker);

        grappa2d_convert_to_convolution_kernel(ker, kRO, kE1, oE1, convKer);

    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_calib_convolution_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void grappa2d_calib_convolution_kernel(const ho3DArray< std::complex<float> >& acsSrc, const ho3DArray< std::complex<float> >& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, ho4DArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa2d_calib_convolution_kernel(const ho3DArray< std::complex<double> >& acsSrc, const ho3DArray< std::complex<double> >& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, ho4DArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask, size_t accelFactor, double thres, size_t kRO, size_t kNE1, ho4DArray<T>& convKer)
{
    try
    {
        bool fitItself = false;
        if (&dataSrc != &dataDst) fitItself = true;

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

        GADGET_CHECK_THROW(endRO>startRO);
        GADGET_CHECK_THROW(endE1>startE1 + accelFactor);

        size_t acsRO = endRO - startRO + 1;
        size_t acsE1 = endE1 - startE1 + 1;

        ho3DArray<T> acsSrc(acsRO, acsE1, srcCHA), acsDst(acsRO, acsE1, dstCHA);

        for (scha = 0; scha < srcCHA; scha++)
        {
            for (e1 = startE1; e1 <= endE1; e1++)
            {
                for (ro = startRO; ro <= endRO; ro++)
                {
                    acsSrc(ro, e1, scha) = dataSrc(ro, e1, scha);
                }
            }
        }

        for (dcha = 0; dcha < dstCHA; dcha++)
        {
            for (e1 = startE1; e1 <= endE1; e1++)
            {
                for (ro = startRO; ro <= endRO; ro++)
                {
                    acsDst(ro, e1, dcha) = dataDst(ro, e1, dcha);
                }
            }
        }

        GADGET_CATCH_THROW(grappa2d_calib_convolution_kernel(acsSrc, acsDst, accelFactor, thres, kRO, kNE1, convKer));
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_calib_convolution_kernel(dataMask) ... ");
    }
}

template <typename T> EXPORTMRICORE void grappa2d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask, size_t accelFactor, double thres, size_t kRO, size_t kNE1, ho4DArray<T>& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_image_domain_kernel(const ho4DArray<T>& convKer, size_t RO, size_t E1, hoNDArray<T>& kIm)
{
    try
    {
        ho4DArray<T> convKerScaled(convKer);
        Gadgetron::scal((typename realType<T>::Type)(std::sqrt((double)(RO*E1))), convKerScaled);
        Gadgetron::zeropad2D(convKerScaled, RO, E1, kIm);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kIm);
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa2d_image_domain_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void grappa2d_image_domain_kernel(const ho4DArray< std::complex<float> >& convKer, size_t RO, size_t E1, hoNDArray< std::complex<float> >& kIm);
template EXPORTMRICORE void grappa2d_image_domain_kernel(const ho4DArray< std::complex<double> >& convKer, size_t RO, size_t E1, hoNDArray< std::complex<double> >& kIm);

// ------------------------------------------------------------------------

template <typename T>
void grappa2d_unmixing_coeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, size_t acceFactorE1, hoNDArray<T>& unmixCoeff, hoNDArray< typename realType<T>::Type >& gFactor)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t srcCHA = kerIm.get_size(2);
        size_t dstCHA = kerIm.get_size(3);

        GADGET_CHECK_THROW(acceFactorE1 >= 1);

        GADGET_CHECK_THROW(coilMap.get_size(0) == RO);
        GADGET_CHECK_THROW(coilMap.get_size(1) == E1);
        GADGET_CHECK_THROW(coilMap.get_size(2) == dstCHA);

        unmixCoeff.create(RO, E1, srcCHA);
        Gadgetron::clear(&unmixCoeff);

        gFactor.create(RO, E1);
        Gadgetron::clear(&gFactor);

        int src;

        T* pKerIm = const_cast<T*>(kerIm.begin());
        T* pCoilMap = const_cast<T*>(coilMap.begin());
        T* pCoeff = unmixCoeff.begin();

        std::vector<size_t> dim(2);
        dim[0] = RO;
        dim[1] = E1;

#pragma omp parallel default(none) private(src) shared(RO, E1, srcCHA, dstCHA, pKerIm, pCoilMap, pCoeff, dim)
        {
            hoNDArray<T> coeff2D, coeffTmp(&dim);
            hoNDArray<T> coilMap2D;
            hoNDArray<T> kerIm2D;

#pragma omp for
            for (src = 0; src<(int)srcCHA; src++)
            {
                coeff2D.create(&dim, pCoeff + src*RO*E1);

                for (size_t dst = 0; dst<dstCHA; dst++)
                {
                    kerIm2D.create(&dim, pKerIm + src*RO*E1 + dst*RO*E1*srcCHA);
                    coilMap2D.create(&dim, pCoilMap + dst*RO*E1);
                    Gadgetron::multiplyConj(kerIm2D, coilMap2D, coeffTmp);
                    Gadgetron::add(coeff2D, coeffTmp, coeff2D);
                }
            }
        }

        hoNDArray<T> conjUnmixCoeff(unmixCoeff);
        Gadgetron::multiplyConj(unmixCoeff, conjUnmixCoeff, conjUnmixCoeff);
        // Gadgetron::sumOverLastDimension(conjUnmixCoeff, gFactor);

        hoNDArray<T> gFactorBuf(RO, E1, 1);
        Gadgetron::sum_over_dimension(conjUnmixCoeff, gFactorBuf, 2);
        Gadgetron::sqrt(gFactorBuf, gFactorBuf);
        Gadgetron::scal((value_type)(1.0 / acceFactorE1), gFactorBuf);

        Gadgetron::complex_to_real(gFactorBuf, gFactor);
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_unmixing_coeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor) ... ");
    }
}

template EXPORTMRICORE void grappa2d_unmixing_coeff(const hoNDArray< std::complex<float> >& kerIm, const hoNDArray< std::complex<float> >& coilMap, size_t acceFactorE1, hoNDArray< std::complex<float> >& unmixCoeff, hoNDArray<float>& gFactor);
template EXPORTMRICORE void grappa2d_unmixing_coeff(const hoNDArray< std::complex<double> >& kerIm, const hoNDArray< std::complex<double> >& coilMap, size_t acceFactorE1, hoNDArray< std::complex<double> >& unmixCoeff, hoNDArray<double>& gFactor);

// ------------------------------------------------------------------------

template <typename T>
void apply_unmix_coeff_kspace(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_THROW(kspace.get_size(0) == unmixCoeff.get_size(0));
        GADGET_CHECK_THROW(kspace.get_size(1) == unmixCoeff.get_size(1));
        GADGET_CHECK_THROW(kspace.get_size(2) == unmixCoeff.get_size(2));

        hoNDArray<T> buffer2DT(kspace);

        GADGET_CHECK_THROW(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspace, buffer2DT));
        Gadgetron::apply_unmix_coeff_aliased_image(buffer2DT, unmixCoeff, complexIm);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_unmix_coeff_kspace(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
    }
}

// ------------------------------------------------------------------------

template <typename T>
void apply_unmix_coeff_aliased_image(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_THROW(aliasedIm.get_size(0) == unmixCoeff.get_size(0));
        GADGET_CHECK_THROW(aliasedIm.get_size(1) == unmixCoeff.get_size(1));
        GADGET_CHECK_THROW(aliasedIm.get_size(2) == unmixCoeff.get_size(2));

        boost::shared_ptr< std::vector<size_t> > dim = aliasedIm.get_dimensions();

        std::vector<size_t> dimIm(*dim);
        dimIm[2] = 1;

        if (!complexIm.dimensions_equal(&dimIm))
        {
            complexIm.create(&dimIm);
        }
        Gadgetron::clear(&complexIm);

        hoNDArray<T> buffer2DT(aliasedIm);

        Gadgetron::multipleMultiply(aliasedIm, unmixCoeff, buffer2DT);
        Gadgetron::sum_over_dimension(buffer2DT, complexIm, 2);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_unmix_coeff_aliased_image(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
    }
}
}
