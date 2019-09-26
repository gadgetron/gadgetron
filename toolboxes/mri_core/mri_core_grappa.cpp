
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
#include "mri_core_utility.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "ImageIOAnalyze.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

void grappa2d_kerPattern(std::vector<int>& kE1, std::vector<int>& oE1, size_t& convKRO, size_t& convKE1, size_t accelFactor, size_t kRO, size_t kNE1, bool fitItself)
{
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

    convKRO = 2 * kRO + 3;

    long long maxKE1 = std::abs(kE1[0]);
    if (std::abs(kE1[kNE1 - 1]) > maxKE1)
    {
        maxKE1 = std::abs(kE1[kNE1 - 1]);
    }
    convKE1 = 2 * maxKE1 + 1;

    return;
}

// ------------------------------------------------------------------------

template <typename T> EXPORTMRICORE void grappa2d_prepare_calib(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray<T>& A_mem, hoNDArray<T>& B_mem)
{
    try
    {
        GADGET_CHECK_THROW(acsSrc.get_size(0) == acsDst.get_size(0));
        GADGET_CHECK_THROW(acsSrc.get_size(1) == acsDst.get_size(1));
        GADGET_CHECK_THROW(acsSrc.get_size(2) >= acsDst.get_size(2));

        size_t RO = acsSrc.get_size(0);
        size_t E1 = acsSrc.get_size(1);
        size_t srcCHA = acsSrc.get_size(2);
        size_t dstCHA = acsDst.get_size(2);

        const T* pSrc = acsSrc.begin();
        const T* pDst = acsDst.begin();

        long long kROhalf = kRO / 2;
        if (2 * kROhalf == kRO)
        {
            GWARN_STREAM("grappa2d_prepare_calib(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2 * kROhalf + 1;

        size_t kNE1 = kE1.size();
        size_t oNE1 = oE1.size();

        /// loop over the calibration region and assemble the equation
        /// Ax = b

        size_t sRO = startRO + kROhalf;
        size_t eRO = endRO - kROhalf;
        size_t sE1 = std::abs(kE1[0]) + startE1;
        size_t eE1 = endE1 - kE1[kNE1 - 1];

        size_t lenRO = eRO - sRO + 1;

        size_t rowA = (eE1 - sE1 + 1)*lenRO;
        size_t colA = kRO * kNE1*srcCHA;
        size_t colB = dstCHA * oNE1;

        hoMatrix<T> A;
        hoMatrix<T> B;

        A_mem.create(rowA, colA);
        A.createMatrix(rowA, colA, A_mem.begin());
        T* pA = A.begin();

        B_mem.create(rowA, colB);
        B.createMatrix(A.rows(), colB, B_mem.begin());
        T* pB = B.begin();

        long long e1;
        for (e1 = (long long)sE1; e1 <= (long long)eE1; e1++)
        {
            for (long long ro = sRO; ro <= (long long)eRO; ro++)
            {
                long long rInd = (e1 - sE1)*lenRO + ro - kROhalf;

                size_t src, dst, ke1, oe1;
                long long kro;

                /// fill matrix A
                size_t col = 0;
                size_t offset = 0;
                for (src = 0; src<srcCHA; src++)
                {
                    for (ke1 = 0; ke1<kNE1; ke1++)
                    {
                        offset = src * RO*E1 + (e1 + kE1[ke1])*RO;
                        for (kro = -kROhalf; kro <= kROhalf; kro++)
                        {
                            /// A(rInd, col++) = acsSrc(ro+kro, e1+kE1[ke1], src);
                            pA[rInd + col * rowA] = pSrc[ro + kro + offset];
                            col++;
                        }
                    }
                }

                /// fill matrix B
                col = 0;
                for (oe1 = 0; oe1<oNE1; oe1++)
                {
                    for (dst = 0; dst<dstCHA; dst++)
                    {
                        B(rInd, col++) = acsDst(ro, e1 + oE1[oe1], dst);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_prepare_calib(...) ... ");
    }
}

template EXPORTMRICORE void grappa2d_prepare_calib(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<float> >& A_mem, hoNDArray< std::complex<float> >& B_mem);
template EXPORTMRICORE void grappa2d_prepare_calib(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<double> >& A_mem, hoNDArray< std::complex<double> >& B_mem);

// ------------------------------------------------------------------------

template <typename T> EXPORTMRICORE void grappa2d_perform_calib(const hoNDArray<T>& A, const hoNDArray<T>& B, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, double thres, hoNDArray<T>& ker)
{
    try
    {
        size_t M = A.get_size(0);
        size_t K = A.get_size(1);

        size_t KB = B.get_size(1);

        GADGET_CHECK_THROW(M == B.get_size(0));

        long long kROhalf = kRO / 2;
        if (2 * kROhalf == kRO)
        {
            GWARN_STREAM("grappa<T>::calib(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2 * kROhalf + 1;

        size_t kNE1 = kE1.size();
        size_t oNE1 = oE1.size();

        size_t srcCHA = K / (kRO*kNE1);
        size_t dstCHA = KB / oNE1;

        ker.create(kRO, kNE1, srcCHA, dstCHA, oNE1);

        hoNDArray<T> x;
        x.create(K, KB);

        SolveLinearSystem_Tikhonov( const_cast< hoNDArray<T>& >(A), const_cast< hoNDArray<T>& >(B), x, thres);
        memcpy(ker.begin(), x.begin(), ker.get_number_of_bytes());

        for (size_t kk = 0; kk>ker.get_number_of_elements(); kk++)
        {
            if (std::isnan(ker(kk).real()) || std::isnan(ker(kk).imag()))
            {
                GERROR_STREAM("nan detected in grappa2d_perform_calib ker ... ");
                throw;
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_perform_calib(...) ... ");
    }
}

template EXPORTMRICORE void grappa2d_perform_calib(const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& B, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, double thres, hoNDArray< std::complex<float> >& ker);
template EXPORTMRICORE void grappa2d_perform_calib(const hoNDArray< std::complex<double> >& A, const hoNDArray< std::complex<double> >& B, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, double thres, hoNDArray< std::complex<double> >& ker);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_calib(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, double thres, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray<T>& ker)
{
    try
    {
        GADGET_CHECK_THROW(acsSrc.get_size(0)==acsDst.get_size(0));
        GADGET_CHECK_THROW(acsSrc.get_size(1)==acsDst.get_size(1));
        GADGET_CHECK_THROW(acsSrc.get_size(2)>=acsDst.get_size(2));

        hoNDArray<T> A, B;
        Gadgetron::grappa2d_prepare_calib(acsSrc, acsDst, kRO, kE1, oE1, startRO, endRO, startE1, endE1, A, B);
        Gadgetron::grappa2d_perform_calib(A, B, kRO, kE1, oE1, thres, ker);
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa2d_calib(...) ... ");
    }
}

template EXPORTMRICORE void grappa2d_calib(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<float> >& ker);
template EXPORTMRICORE void grappa2d_calib(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<double> >& ker);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_prepare_recon(const hoNDArray<T>& kspace, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, bool periodic_boundary_condition, hoNDArray<T>& A, hoNDArray<unsigned short>& AInd)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t srcCHA = kspace.get_size(2);

        size_t kNE1 = kE1.size();
        size_t oNE1 = oE1.size();

        long long kROhalf = kRO / 2;
        if (2 * kROhalf == kRO)
        {
            GWARN_STREAM("grappa2d_prepare_recon(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2 * kROhalf + 1;

        size_t R = (size_t)(kE1[1] - kE1[0]);

        long long kro;
        size_t ro, e1, src, ke1;

        // first first sampled line along E1
        size_t startE1 = E1, endE1 = 0;
        for (e1 = 0; e1 < E1; e1++)
        {
            value_type v = std::abs(kspace(RO / 2, e1));

            if(v>0)
            {
                if (e1<startE1) startE1 = e1;
                if (e1>endE1) endE1 = e1;
            }
        }

        size_t startRO = RO, endRO = 0;
        for (ro = 0; ro < RO; ro++)
        {
            value_type v1 = std::abs(kspace(ro, startE1));
            value_type v2 = std::abs(kspace(ro, endE1));

            if (v1>0 || v2>0)
            {
                if (ro<startRO) startRO = ro;
                if (ro>endRO) endRO = ro;
            }
        }

        // start with first sampled line, assemble the data matrix
        size_t num_E1 = 0;
        for (e1 = startE1; e1 <= endE1; e1 += R) num_E1++;

        size_t rowA = num_E1 * (endRO-startRO+1); // for every position in the valid kspace loc
        size_t colA = kRO * kNE1 * srcCHA; // for every position in the valid kspace loc

        A.create(rowA, colA);
        AInd.create(rowA, 2);

        num_E1 = 0;
        for (e1 = startE1; e1 <= endE1; e1 += R)
        {
            for (ro=startRO; ro<=endRO; ro++)
            {
                size_t row = num_E1 * (endRO - startRO + 1) + (ro-startRO);

                size_t col = 0;
                // loop over srcCHA
                for (src = 0; src<srcCHA; src++)
                {
                    if (periodic_boundary_condition)
                    {
                        for (ke1 = 0; ke1 < kNE1; ke1++)
                        {
                            long long src_e1 = e1 + kE1[ke1];
                            if (src_e1 < 0) src_e1 += E1;
                            if (src_e1 >= E1) src_e1 -= E1;

                            for (kro = -kROhalf; kro <= kROhalf; kro++)
                            {
                                long long src_ro = ro + kro;
                                if (src_ro < 0) src_ro += RO;
                                if (src_ro >= RO) src_ro -= RO;

                                A(row, col) = kspace(src_ro, src_e1, src);
                                col++;
                            }
                        }
                    }
                    else
                    {
                        for (ke1 = 0; ke1 < kNE1; ke1++)
                        {
                            long long src_e1 = e1 + kE1[ke1];

                            for (kro = -kROhalf; kro <= kROhalf; kro++)
                            {
                                long long src_ro = ro + kro;

                                if ( (src_e1 < 0) || (src_e1 >= E1) || (src_ro < 0) || (src_ro >= RO) )
                                {
                                    A(row, col) = 0;
                                }
                                else
                                {
                                    A(row, col) = kspace(src_ro, src_e1, src);
                                }

                                col++;
                            }
                        }
                    }
                }

                AInd(row, 0) = ro;
                AInd(row, 1) = e1;
            }

            num_E1++;
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_prepare_recon(...) ... ");
    }
}

template EXPORTMRICORE void grappa2d_prepare_recon(const hoNDArray< std::complex<float> >& kspace, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, bool periodic_boundary_condition, hoNDArray< std::complex<float> >& A, hoNDArray<unsigned short>& AInd);
template EXPORTMRICORE void grappa2d_prepare_recon(const hoNDArray< std::complex<double> >& kspace, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, bool periodic_boundary_condition, hoNDArray< std::complex<double> >& A, hoNDArray<unsigned short>& AInd);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_perform_recon(const hoNDArray<T>& A, const hoNDArray<T>& ker, const hoNDArray<unsigned short>& AInd, const std::vector<int>& oE1, size_t RO, size_t E1, hoNDArray<T>& res)
{
    try
    {
        size_t M = A.get_size(0);
        size_t K = A.get_size(1);

        size_t KB = ker.get_number_of_elements() / K;

        hoNDArray<T> x;
        x.create(K, KB, const_cast<T*>(ker.begin()));

        hoNDArray<T> recon;
        Gadgetron::gemm(recon, A, false, x, false);

        Gadgetron::grappa2d_fill_reconed_kspace(AInd, recon, oE1, RO, E1, res);
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_perform_recon(...) ... ");
    }
}

template EXPORTMRICORE void grappa2d_perform_recon(const hoNDArray< std::complex<float> >& A, const hoNDArray< std::complex<float> >& ker, const hoNDArray<unsigned short>& AInd, const std::vector<int>& oE1, size_t RO, size_t E1, hoNDArray< std::complex<float> >& res);
template EXPORTMRICORE void grappa2d_perform_recon(const hoNDArray< std::complex<double> >& A, const hoNDArray< std::complex<double> >& ker, const hoNDArray<unsigned short>& AInd, const std::vector<int>& oE1, size_t RO, size_t E1, hoNDArray< std::complex<double> >& res);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_fill_reconed_kspace(const hoNDArray<unsigned short>& AInd, const hoNDArray<T>& recon, const std::vector<int>& oE1, size_t RO, size_t E1, hoNDArray<T>& res)
{
    try
    {
        size_t M = recon.get_size(0);
        size_t KB = recon.get_size(1);

        // for every row in data matrix, assign the recon data back to kspace
        size_t oNE1 = oE1.size();
        size_t dstCHA = KB / oNE1;

        if (res.get_size(0) != RO || res.get_size(1) != E1 || res.get_size(2) != dstCHA)
        {
            res.create(RO, E1, dstCHA);
            Gadgetron::clear(res);
        }

        size_t r, c, ro, e1, dst, oe1, de1;
        for (r = 0; r < M; r++)
        {
            ro = AInd(r, 0);
            e1 = AInd(r, 1);

            for (oe1 = 0; oe1<oNE1; oe1++)
            {
                de1 = e1 + oE1[oe1];
                if (de1 > E1) de1 -= E1;

                for (dst = 0; dst<dstCHA; dst++)
                {
                    res(ro, de1, dst) = recon(r, dst + oe1 * dstCHA);
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_fill_reconed_kspace(...) ... ");
    }
}

template EXPORTMRICORE void grappa2d_fill_reconed_kspace(const hoNDArray<unsigned short>& AInd, const hoNDArray< std::complex<float> >& recon, const std::vector<int>& oE1, size_t RO, size_t E1, hoNDArray< std::complex<float> >& res);
template EXPORTMRICORE void grappa2d_fill_reconed_kspace(const hoNDArray<unsigned short>& AInd, const hoNDArray< std::complex<double> >& recon, const std::vector<int>& oE1, size_t RO, size_t E1, hoNDArray< std::complex<double> >& res);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_recon(const hoNDArray<T>& kspace, const hoNDArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, bool periodic_boundary_condition, hoNDArray<T>& res)
{
    try
    {
        hoNDArray<T> A;
        hoNDArray<unsigned short> AInd;
        Gadgetron::grappa2d_prepare_recon(kspace, kRO, kE1, oE1, periodic_boundary_condition, A, AInd);
        Gadgetron::grappa2d_perform_recon(A, ker, AInd, oE1, kspace.get_size(0), kspace.get_size(1), res);
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_recon(...) ... ");
    }
}

template EXPORTMRICORE void grappa2d_recon(const hoNDArray< std::complex<float> >& kspace, const hoNDArray< std::complex<float> >& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, bool periodic_boundary_condition, hoNDArray< std::complex<float> >& res);
template EXPORTMRICORE void grappa2d_recon(const hoNDArray< std::complex<double> >& kspace, const hoNDArray< std::complex<double> >& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, bool periodic_boundary_condition, hoNDArray< std::complex<double> >& res);

// ------------------------------------------------------------------------

template <typename T>
void grappa2d_convert_to_convolution_kernel(const hoNDArray<T>& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, hoNDArray<T>& convKer)
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
        convKer.create(convKRO, convKE1, srcCHA, dstCHA);
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

template EXPORTMRICORE void grappa2d_convert_to_convolution_kernel(const hoNDArray< std::complex<float> >& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa2d_convert_to_convolution_kernel(const hoNDArray< std::complex<double> >& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T>
void grappa2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray<T>& convKer)
{
        std::vector<int> kE1, oE1;

        bool fitItself = false;
        if (&acsSrc != &acsDst) fitItself = true;

        size_t convkRO, convkE1;

        grappa2d_kerPattern(kE1, oE1, convkRO, convkE1, accelFactor, kRO, kNE1, fitItself);

        hoNDArray<T> ker;
        grappa2d_calib(acsSrc, acsDst, thres, kRO, kE1, oE1, startRO, endRO, startE1, endE1, ker);

        grappa2d_convert_to_convolution_kernel(ker, kRO, kE1, oE1, convKer);
}

template EXPORTMRICORE void grappa2d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa2d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, size_t startRO, size_t endRO, size_t startE1, size_t endE1, hoNDArray< std::complex<double> >& convKer);


// ------------------------------------------------------------------------

template <typename T>
void grappa2d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, hoNDArray<T>& convKer)
{
    size_t startRO = 0;
    size_t endRO = acsSrc.get_size(0) - 1;
    size_t startE1 = 0;
    size_t endE1 = acsSrc.get_size(1) - 1;

    grappa2d_calib_convolution_kernel(acsSrc, acsDst, accelFactor, thres, kRO, kNE1, startRO, endRO, startE1, endE1, convKer);
}

template EXPORTMRICORE void grappa2d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa2d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, size_t accelFactor, double thres, size_t kRO, size_t kNE1, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask, size_t accelFactor, double thres, size_t kRO, size_t kNE1, hoNDArray<T>& convKer)
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

        GADGET_CHECK_THROW(endRO>startRO);
        GADGET_CHECK_THROW(endE1>startE1 + accelFactor);

        GADGET_CATCH_THROW(grappa2d_calib_convolution_kernel(dataSrc, dataDst, accelFactor, thres, kRO, kNE1, startRO, endRO, startE1, endE1, convKer));
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_calib_convolution_kernel(dataMask) ... ");
    }
}

template EXPORTMRICORE void grappa2d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& dataSrc, const hoNDArray< std::complex<float> >& dataDst, hoNDArray<unsigned short>& dataMask, size_t accelFactor, double thres, size_t kRO, size_t kNE1, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa2d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& dataSrc, const hoNDArray< std::complex<double> >& dataDst, hoNDArray<unsigned short>& dataMask, size_t accelFactor, double thres, size_t kRO, size_t kNE1, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, hoNDArray<T>& kIm)
{
    try
    {
        hoNDArray<T> convKerScaled(convKer);
        Gadgetron::scal((typename realType<T>::Type)(std::sqrt((double)(RO*E1))), convKerScaled);
        Gadgetron::pad(RO, E1, convKerScaled, kIm);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kIm);
    }
    catch(...)
    {
        GADGET_THROW("Errors in grappa2d_image_domain_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void grappa2d_image_domain_kernel(const hoNDArray< std::complex<float> >& convKer, size_t RO, size_t E1, hoNDArray< std::complex<float> >& kIm);
template EXPORTMRICORE void grappa2d_image_domain_kernel(const hoNDArray< std::complex<double> >& convKer, size_t RO, size_t E1, hoNDArray< std::complex<double> >& kIm);

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

        std::vector<size_t> dimUnmixing(3);
        dimUnmixing[0] = RO; dimUnmixing[1] = E1; dimUnmixing[2] = srcCHA;
        if (!unmixCoeff.dimensions_equal(&dimUnmixing))
        {
            unmixCoeff.create(RO, E1, srcCHA);
        }
        Gadgetron::clear(&unmixCoeff);

        std::vector<size_t> dimGFactor(2);
        dimGFactor[0] = RO; dimGFactor[1] = E1;
        if (!gFactor.dimensions_equal(&dimGFactor))
        {
            gFactor.create(RO, E1);
        }
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
void grappa2d_image_domain_unwrapping(const hoNDArray<T>& kspace, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm)
{
    try
    {
        hoNDArray<T> aliasedIm(kspace);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspace, aliasedIm);

        Gadgetron::grappa2d_image_domain_unwrapping_aliased_image(aliasedIm, kerIm, complexIm);
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_image_domain_unwrapping(const hoNDArray<T>& kspace, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm) ... ");
    }
}

template EXPORTMRICORE void grappa2d_image_domain_unwrapping(const hoNDArray< std::complex<float> >& kspace, const hoNDArray< std::complex<float> >& kerIm, hoNDArray< std::complex<float> >& complexIm);
template EXPORTMRICORE void grappa2d_image_domain_unwrapping(const hoNDArray< std::complex<double> >& kspace, const hoNDArray< std::complex<double> >& kerIm, hoNDArray< std::complex<double> >& complexIm);

// ------------------------------------------------------------------------

template <typename T> 
void grappa2d_image_domain_unwrapping_aliased_image(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm)
{
    try
    {
        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t srcCHA = kerIm.get_size(2);
        size_t dstCHA = kerIm.get_size(3);

        GADGET_CHECK_THROW(aliasedIm.get_size(0) == RO);
        GADGET_CHECK_THROW(aliasedIm.get_size(1) == E1);
        GADGET_CHECK_THROW(aliasedIm.get_size(2) == srcCHA);

        std::vector<size_t> dim;
        aliasedIm.get_dimensions(dim);

        std::vector<size_t> dimIm(dim);
        dimIm[2] = dstCHA;

        if (!complexIm.dimensions_equal(&dimIm))
        {
            complexIm.create(&dimIm);
        }

        size_t num = aliasedIm.get_number_of_elements() / (RO*E1*srcCHA);

        long long n;

#pragma omp parallel default(none) private(n) shared(kerIm, num, aliasedIm, RO, E1, srcCHA, dstCHA, complexIm) if(num>=16)
        {
            hoNDArray<T> unwrappedBuffer;
            unwrappedBuffer.create(RO, E1, srcCHA);

            hoNDArray<T> unwrappedIm2D;
            unwrappedIm2D.create(RO, E1, 1);

#pragma omp for 
            for (n = 0; n < (long long)num; n++)
            {
                hoNDArray<T> bufIm(RO, E1, srcCHA, const_cast<T*>(aliasedIm.begin() + n*RO*E1*srcCHA));

                for (size_t dcha = 0; dcha < dstCHA; dcha++)
                {
                    hoNDArray<T> kerSrcCha(RO, E1, srcCHA, const_cast<T*>(kerIm.begin() + dcha*RO*E1*srcCHA));

                    Gadgetron::multiply(kerSrcCha, bufIm, unwrappedBuffer);
                    Gadgetron::sum_over_dimension(unwrappedBuffer, unwrappedIm2D, 2);

                    memcpy(complexIm.begin() + n*RO*E1*dstCHA + dcha*RO*E1, unwrappedIm2D.begin(), sizeof(T)*RO*E1);
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa2d_image_domain_unwrapping_aliased_image(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm) ... ");
    }
}

template EXPORTMRICORE void grappa2d_image_domain_unwrapping_aliased_image(const hoNDArray< std::complex<float> >& aliasedIm, const hoNDArray< std::complex<float> >& kerIm, hoNDArray< std::complex<float> >& complexIm);
template EXPORTMRICORE void grappa2d_image_domain_unwrapping_aliased_image(const hoNDArray< std::complex<double> >& aliasedIm, const hoNDArray< std::complex<double> >& kerIm, hoNDArray< std::complex<double> >& complexIm);

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
        GADGET_CATCH_THROW(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspace, buffer2DT));

        std::vector<size_t> dim;
        kspace.get_dimensions(dim);
        dim[2] = 1;

        if (!complexIm.dimensions_equal(&dim))
        {
            complexIm.create(&dim);
        }

        Gadgetron::multiply(buffer2DT, unmixCoeff, buffer2DT);
        Gadgetron::sum_over_dimension(buffer2DT, complexIm, 2);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_unmix_coeff_kspace(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
    }
}

template EXPORTMRICORE void apply_unmix_coeff_kspace(const hoNDArray< std::complex<float> >& kspace, const hoNDArray< std::complex<float> >& unmixCoeff, hoNDArray< std::complex<float> >& complexIm);
template EXPORTMRICORE void apply_unmix_coeff_kspace(const hoNDArray< std::complex<double> >& kspace, const hoNDArray< std::complex<double> >& unmixCoeff, hoNDArray< std::complex<double> >& complexIm);

// ------------------------------------------------------------------------

template <typename T>
void apply_unmix_coeff_aliased_image(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_THROW(aliasedIm.get_size(0) == unmixCoeff.get_size(0));
        GADGET_CHECK_THROW(aliasedIm.get_size(1) == unmixCoeff.get_size(1));
        GADGET_CHECK_THROW(aliasedIm.get_size(2) == unmixCoeff.get_size(2));

        std::vector<size_t> dim;
        aliasedIm.get_dimensions(dim);
        dim[2] = 1;

        if (!complexIm.dimensions_equal(&dim))
        {
            complexIm.create(&dim);
        }

        hoNDArray<T> buffer2DT(aliasedIm);

        Gadgetron::multiply(aliasedIm, unmixCoeff, buffer2DT);
        Gadgetron::sum_over_dimension(buffer2DT, complexIm, 2);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_unmix_coeff_aliased_image(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
    }
}

template EXPORTMRICORE void apply_unmix_coeff_aliased_image(const hoNDArray< std::complex<float> >& aliasedIm, const hoNDArray< std::complex<float> >& unmixCoeff, hoNDArray< std::complex<float> >& complexIm);
template EXPORTMRICORE void apply_unmix_coeff_aliased_image(const hoNDArray< std::complex<double> >& aliasedIm, const hoNDArray< std::complex<double> >& unmixCoeff, hoNDArray< std::complex<double> >& complexIm);

// ------------------------------------------------------------------------

void grappa3d_kerPattern(std::vector<int>& kE1, std::vector<int>& oE1,
                    std::vector<int>& kE2, std::vector<int>& oE2,
                    size_t& convKRO, size_t& convKE1, size_t& convKE2,
                    size_t accelFactorE1, size_t accelFactorE2,
                    size_t kRO, size_t kNE1, size_t kNE2, bool fitItself)
{
    grappa2d_kerPattern(kE1, oE1, convKRO, convKE1, accelFactorE1, kRO, kNE1, fitItself);
    grappa2d_kerPattern(kE2, oE2, convKRO, convKE2, accelFactorE2, kRO, kNE2, fitItself);
}

template <typename T> 
void grappa3d_calib(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                double thres, double overDetermineRatio, size_t kRO,
                const std::vector<int>& kE1, const std::vector<int>& oE1,
                const std::vector<int>& kE2, const std::vector<int>& oE2,
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

        const T* pSrc = acsSrc.begin();
        const T* pDst = acsDst.begin();

        long long kROhalf = (long long)kRO / 2;
        if (2 * kROhalf == kRO)
        {
            GWARN_STREAM("grappa3d_calib(...) - 2*kROhalf == kRO " << kRO);
        }

        kRO = 2 * kROhalf + 1;

        size_t kNE1 = kE1.size();
        size_t oNE1 = oE1.size();

        size_t kNE2 = kE2.size();
        size_t oNE2 = oE2.size();

        // allocate kernel
        ker.create(kRO, kNE1, kNE2, srcCHA, dstCHA, oNE1, oNE2);

        // loop over the calibration region and assemble the equation
        // Ax = b

        size_t sRO = kROhalf;
        size_t eRO = RO - kROhalf - 1;

        size_t sE1 = std::abs(kE1[0]);
        size_t eE1 = E1 - 1 - kE1[kNE1 - 1];

        size_t sE2 = std::abs(kE2[0]);
        size_t eE2 = E2 - 1 - kE2[kNE2 - 1];

        size_t lenRO = eRO - kROhalf + 1;
        size_t lenE1 = eE1 - sE1 + 1;
        size_t lenE2 = eE2 - sE2 + 1;

        size_t colA = kRO*kNE1*kNE2*srcCHA;
        size_t colB = dstCHA*oNE1*oNE2;

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
                        GDEBUG_STREAM("grappa3d_calib(...) - overDetermineRatio = " << overDetermineRatio << " ; RO data range used : [" << sRO << " " << eRO << "] ...");
                    }
                    catch (...)
                    {
                        GWARN_STREAM("grappa3d_calib(...) - overDetermineRatio is ignored ... ");
                        throw;
                    }
                }
                catch (...)
                {
                    GWARN_STREAM("grappa3d_calib(...) - overDetermineRatio is ignored ... ");
                    throw;
                }
            }
        }

        size_t rowA = lenRO*lenE1*lenE2;

        hoMatrix<T> A, B, x(colA, colB);

        hoNDArray<T> A_mem(rowA, colA);
        A.createMatrix(rowA, colA, A_mem.begin());
        T* pA = A.begin();

        hoNDArray<T> B_mem(rowA, colB);
        B.createMatrix(rowA, colB, B_mem.begin());
        T* pB = B.begin();

        long long e2;

#pragma omp parallel for default(none) private(e2) shared(sE2, eE2, sE1, eE1, kROhalf, sRO, eRO, lenRO, lenE1, srcCHA, kNE2, kNE1, A, rowA, pA, acsSrc, kE1, kE2, oNE2, oNE1, dstCHA, B, pB, acsDst, oE1, oE2)
        for (e2 = (long long)sE2; e2 <= (long long)eE2; e2++)
        {
            long long e1;
            for (e1 = (long long)sE1; e1 <= (long long)eE1; e1++)
            {
                for (long long ro = (long long)sRO; ro <= (long long)eRO; ro++)
                {
                    size_t rInd = (e2 - sE2)*lenRO*lenE1 + (e1 - sE1)*lenRO + ro - sRO;

                    size_t src, dst, ke1, ke2, oe1, oe2;
                    long long kro;

                    // fill matrix A
                    size_t col = 0;
                    for (src = 0; src<srcCHA; src++)
                    {
                        for (ke2 = 0; ke2<kNE2; ke2++)
                        {
                            for (ke1 = 0; ke1<kNE1; ke1++)
                            {
                                for (kro = -kROhalf; kro <= kROhalf; kro++)
                                {
                                    pA[rInd + col*rowA] = acsSrc(ro + kro, e1 + kE1[ke1], e2 + kE2[ke2], src);
                                    col++;
                                }
                            }
                        }
                    }

                    // fill matrix B
                    col = 0;
                    for (oe2 = 0; oe2<oNE2; oe2++)
                    {
                        for (oe1 = 0; oe1<oNE1; oe1++)
                        {
                            for (dst = 0; dst<dstCHA; dst++)
                            {
                                pB[rInd + col*rowA] = acsDst(ro, e1 + oE1[oe1], e2 + oE2[oe2], dst);
                                col++;
                            }
                        }
                    }
                }
            }
        }

        SolveLinearSystem_Tikhonov(A, B, x, thres);

        memcpy(ker.begin(), x.begin(), ker.get_number_of_bytes());

        for(size_t kk=0; kk>ker.get_number_of_elements(); kk++)
        {
            if(std::isnan(ker(kk).real()) || std::isnan(ker(kk).imag()))
            {
                GERROR_STREAM("nan detected in grappa3d_calib ker ... ");
                throw;
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa3d_calib(...) ... ");
    }

    return;
}

template EXPORTMRICORE void grappa3d_calib(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, double thres, double overDetermineRatio, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, const std::vector<int>& kE2, const std::vector<int>& oE2, hoNDArray< std::complex<float> >& ker);
template EXPORTMRICORE void grappa3d_calib(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, double thres, double overDetermineRatio, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, const std::vector<int>& kE2, const std::vector<int>& oE2, hoNDArray< std::complex<double> >& ker);

// ------------------------------------------------------------------------

template <typename T> 
void grappa3d_convert_to_convolution_kernel(const hoNDArray<T>& ker, 
                                        size_t kRO, 
                                        const std::vector<int>& kE1, const std::vector<int>& oE1, 
                                        const std::vector<int>& kE2, const std::vector<int>& oE2, 
                                        hoNDArray<T>& convKer)
{
    try
    {
        long long srcCHA = (long long)(ker.get_size(3));
        long long dstCHA = (long long)(ker.get_size(4));

        long long kNE1 = (long long)(kE1.size());
        long long oNE1 = (long long)(oE1.size());

        long long kNE2 = (long long)(kE2.size());
        long long oNE2 = (long long)(oE2.size());

        long long kROhalf = kRO / 2;
        if (2 * kROhalf == kRO)
        {
            GWARN_STREAM("grappa3d_convert_to_convolution_kernel(...) - 2*kROhalf == kRO " << kRO);
        }
        kRO = 2 * kROhalf + 1;

        /// fill the convolution kernels
        long long convKRO = 2 * kRO + 3;

        long long maxKE1 = std::abs(kE1[0]);
        if (std::abs(kE1[kNE1 - 1]) > maxKE1)
        {
            maxKE1 = std::abs(kE1[kNE1 - 1]);
        }
        long long convKE1 = 2 * maxKE1 + 1;

        long long maxKE2 = std::abs(kE2[0]);
        if (std::abs(kE2[kNE2 - 1]) > maxKE2)
        {
            maxKE2 = std::abs(kE2[kNE2 - 1]);
        }
        long long convKE2 = 2 * maxKE2 + 1;

        /// allocate the convolution kernel
        convKer.create(convKRO, convKE1, convKE2, srcCHA, dstCHA);
        Gadgetron::clear(&convKer);

        /// index
        long long oe1, oe2, kro, ke1, ke2, src, dst;

        /// fill the convolution kernel and sum up multiple kernels
        for (oe2 = 0; oe2<oNE2; oe2++)
        {
            for (oe1 = 0; oe1<oNE1; oe1++)
            {
                for (ke2 = 0; ke2<kNE2; ke2++)
                {
                    for (ke1 = 0; ke1<kNE1; ke1++)
                    {
                        for (kro = -kROhalf; kro <= kROhalf; kro++)
                        {
                            for (dst = 0; dst<dstCHA; dst++)
                            {
                                for (src = 0; src<srcCHA; src++)
                                {
                                    convKer(-kro + kRO + 1, oE1[oe1] - kE1[ke1] + maxKE1, oE2[oe2] - kE2[ke2] + maxKE2, src, dst) = ker(kro + kROhalf, ke1, ke2, src, dst, oe1, oe2);
                                }
                            }
                        }
                    }
                }
            }
        }

        if ((oE1[0] != 0) && (oE2[0] != 0) && (srcCHA == dstCHA))
        {
            for (dst = 0; dst<dstCHA; dst++)
            {
                for (src = 0; src<srcCHA; src++)
                {
                    if (src == dst)
                    {
                        convKer(kRO + 1, maxKE1, maxKE2, src, dst) = 1.0;
                    }
                    else
                    {
                        convKer(kRO + 1, maxKE1, maxKE2, src, dst) = 0.0;
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in grappa3d_convert_to_convolution_kernel(...) ... ");
        throw;
    }

    return;
}

template EXPORTMRICORE void grappa3d_convert_to_convolution_kernel(const hoNDArray< std::complex<float> >& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, const std::vector<int>& kE2, const std::vector<int>& oE2, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa3d_convert_to_convolution_kernel(const hoNDArray< std::complex<double> >& ker, size_t kRO, const std::vector<int>& kE1, const std::vector<int>& oE1, const std::vector<int>& kE2, const std::vector<int>& oE2, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void grappa3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                                size_t accelFactorE1, size_t accelFactorE2,
                                double thres, double overDetermineRatio,
                                size_t kRO, size_t kNE1, size_t kNE2,
                                size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
                                hoNDArray<T>& convKer)
{
    try
    {
        size_t RO = acsSrc.get_size(0);
        size_t E1 = acsSrc.get_size(1);
        size_t E2 = acsSrc.get_size(2);
        size_t srcCHA = acsSrc.get_size(3);
        size_t dstCHA = acsDst.get_size(3);

        if (startRO == 0 && endRO == RO - 1 && startE1 == 0 && endE1 == E1 - 1 && startE2 == 0 && endE2 == E2 - 1)
        {
            grappa3d_calib_convolution_kernel(acsSrc, acsDst, accelFactorE1, accelFactorE2, thres, overDetermineRatio, kRO, kNE1, kNE2, convKer);
        }
        else
        {
            hoNDArray<T> acsSrcFullSampled, acsDstFullSampled;

            vector_td<size_t, 4> cropOffset, cropSize;
            cropOffset[0] = startRO;
            cropOffset[1] = startE1;
            cropOffset[2] = startE2;
            cropOffset[3] = 0;

            cropSize[0] = endRO-startRO+1;
            cropSize[1] = endE1-startE1+1;
            cropSize[2] = endE2-startE2+1;
            cropSize[3] = srcCHA;

            Gadgetron::crop(cropOffset, cropSize, acsSrc, acsSrcFullSampled);

            cropSize[3] = dstCHA;
            Gadgetron::crop(cropOffset, cropSize, acsDst, acsDstFullSampled);

            grappa3d_calib_convolution_kernel(acsSrcFullSampled, acsDstFullSampled, accelFactorE1, accelFactorE2, thres, overDetermineRatio, kRO, kNE1, kNE2, convKer);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa3d_calib_convolution_kernel(start, end) ... ");
    }

    return;
}

template EXPORTMRICORE void grappa3d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, size_t accelFactorE1, size_t accelFactorE2, double thres, double overDetermineRatio, size_t kRO, size_t kNE1, size_t kNE2, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa3d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, size_t accelFactorE1, size_t accelFactorE2, double thres, double overDetermineRatio, size_t kRO, size_t kNE1, size_t kNE2, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void grappa3d_calib_convolution_kernel(const hoNDArray<T>& acsSrc, const hoNDArray<T>& acsDst,
                                size_t accelFactorE1, size_t accelFactorE2,
                                double thres, double overDetermineRatio,
                                size_t kRO, size_t kNE1, size_t kNE2,
                                hoNDArray<T>& convKer)
{
    try
    {
        std::vector<int> kE1, oE1, kE2, oE2;

        bool fitItself = false;
        if (&acsSrc != &acsDst) fitItself = true;

        size_t convkRO, convkE1, convkE2;

        grappa3d_kerPattern(kE1, oE1, kE2, oE2, convkRO, convkE1, convkE2, accelFactorE1, accelFactorE2, kRO, kNE1, kNE2, fitItself);

        hoNDArray<T> ker;
        grappa3d_calib(acsSrc, acsDst, thres, overDetermineRatio, kRO, kE1, oE1, kE2, oE2, ker);

        grappa3d_convert_to_convolution_kernel(ker, kRO, kE1, oE1, kE2, oE2, convKer);

    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa3d_calib_convolution_kernel(...) ... ");
    }

    return;
}

template EXPORTMRICORE void grappa3d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& acsSrc, const hoNDArray< std::complex<float> >& acsDst, size_t accelFactorE1, size_t accelFactorE2, double thres, double overDetermineRatio, size_t kRO, size_t kNE1, size_t kNE2, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa3d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& acsSrc, const hoNDArray< std::complex<double> >& acsDst, size_t accelFactorE1, size_t accelFactorE2, double thres, double overDetermineRatio, size_t kRO, size_t kNE1, size_t kNE2, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void grappa3d_calib_convolution_kernel(const hoNDArray<T>& dataSrc, const hoNDArray<T>& dataDst, hoNDArray<unsigned short>& dataMask,
                                size_t accelFactorE1, size_t accelFactorE2,
                                double thres, double overDetermineRatio,
                                size_t kRO, size_t kNE1, size_t kNE2,
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
        size_t dstCHA = dataDst.get_size(4);

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

        GADGET_CHECK_THROW(endRO>startRO);
        GADGET_CHECK_THROW(endE1>startE1 + accelFactorE1);
        GADGET_CHECK_THROW(endE2>startE2 + accelFactorE2);

        grappa3d_calib_convolution_kernel(dataSrc, dataDst, accelFactorE1, accelFactorE2, thres, overDetermineRatio, kRO, kNE1, kNE2, startRO, endRO, startE1, endE1, startE2, endE2, convKer);
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa3d_calib_convolution_kernel(dataMask) ... ");
    }

    return;
}

template EXPORTMRICORE void grappa3d_calib_convolution_kernel(const hoNDArray< std::complex<float> >& dataSrc, const hoNDArray< std::complex<float> >& dataDst, hoNDArray<unsigned short>& dataMask, size_t accelFactorE1, size_t accelFactorE2, double thres, double overDetermineRatio, size_t kRO, size_t kNE1, size_t kNE2, hoNDArray< std::complex<float> >& convKer);
template EXPORTMRICORE void grappa3d_calib_convolution_kernel(const hoNDArray< std::complex<double> >& dataSrc, const hoNDArray< std::complex<double> >& dataDst, hoNDArray<unsigned short>& dataMask, size_t accelFactorE1, size_t accelFactorE2, double thres, double overDetermineRatio, size_t kRO, size_t kNE1, size_t kNE2, hoNDArray< std::complex<double> >& convKer);

// ------------------------------------------------------------------------

template <typename T> 
void grappa3d_image_domain_kernel(const hoNDArray<T>& convKer, size_t RO, size_t E1, size_t E2, hoNDArray<T>& kIm, bool preset_kIm_with_zeros)
{
    try
    {
        size_t srcCHA = convKer.get_size(3);
        size_t dstCHA = convKer.get_size(4);

        if (kIm.get_size(0) != RO || kIm.get_size(1) != E1 || kIm.get_size(2) != E2 || kIm.get_size(3) != srcCHA || kIm.get_size(4) != dstCHA)
        {
            kIm.create(RO, E1, E2, srcCHA, dstCHA);
        }

        hoNDArray<T> convKerScaled;
        convKerScaled = convKer;

        Gadgetron::scal((typename realType<T>::Type)(std::sqrt((double)(RO*E1*E2))), convKerScaled);
        Gadgetron::pad(RO, E1, E2, convKerScaled, kIm, preset_kIm_with_zeros);

        long long n;

    #pragma omp parallel default(none) private(n) shared(RO, E1, E2, srcCHA, dstCHA, kIm)
        {
            hoNDArray<T> kImTmp(RO, E1, E2);
            hoNDArray<T> kImRes(RO, E1, E2);

    #pragma omp for 
            for (n = 0; n < (long long)(srcCHA*dstCHA); n++)
            {
                long long d = n / srcCHA;
                long long s = n - d*srcCHA;

                T* pkImCha = kIm.begin() + d*RO*E1*E2*srcCHA + s*RO*E1*E2;

                hoNDArray<T> kImCha(RO, E1, E2, pkImCha);
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kImCha, kImRes, kImTmp);
                memcpy(pkImCha, kImRes.begin(), kImRes.get_number_of_bytes());
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa3d_image_domain_kernel(...) ... ");
    }
}

template EXPORTMRICORE void grappa3d_image_domain_kernel(const hoNDArray< std::complex<float> >& convKer, size_t RO, size_t E1, size_t E2, hoNDArray< std::complex<float> >& kIm, bool preset_kIm_with_zeros);
template EXPORTMRICORE void grappa3d_image_domain_kernel(const hoNDArray< std::complex<double> >& convKer, size_t RO, size_t E1, size_t E2, hoNDArray< std::complex<double> >& kIm, bool preset_kIm_with_zeros);

// ------------------------------------------------------------------------

template <typename T> 
void grappa3d_unmixing_coeff(const hoNDArray<T>& convKer, const hoNDArray<T>& coilMap,
                        size_t acceFactorE1, size_t acceFactorE2, hoNDArray<T>& unmixCoeff, 
                        hoNDArray< typename realType<T>::Type >& gFactor)
{
    try
    {
        size_t kRO = convKer.get_size(0);
        size_t kE1 = convKer.get_size(1);
        size_t kE2 = convKer.get_size(2);

        size_t RO = coilMap.get_size(0);
        size_t E1 = coilMap.get_size(1);
        size_t E2 = coilMap.get_size(2);
        size_t dstCHA = coilMap.get_size(3);

        size_t srcCHA = convKer.get_size(3);

        GADGET_CHECK_THROW(convKer.get_size(4) == dstCHA);

        if (unmixCoeff.get_size(0)!=RO 
            || unmixCoeff.get_size(1) != E1 
            || unmixCoeff.get_size(2) != E2 
            || unmixCoeff.get_size(3) != srcCHA )
        {
            unmixCoeff.create(RO, E1, E2, srcCHA);
        }

        Gadgetron::clear(&unmixCoeff);

        if (gFactor.get_size(0) != RO
            || gFactor.get_size(1) != E1
            || gFactor.get_size(2) != E2)
        {
            gFactor.create(RO, E1, E2);
        }

        hoNDArray<T> convKerScaled;
        convKerScaled = convKer;

        Gadgetron::scal((typename realType<T>::Type)(std::sqrt((double)(RO*E1*E2))), convKerScaled);

        long long scha;

        #pragma omp parallel private(scha) shared(RO, E1, E2, srcCHA, dstCHA, kRO, kE1, kE2, convKerScaled, coilMap, unmixCoeff)
        {

            hoNDArray<T> convKerCha;
            hoNDArray<T> convKerChaPadded;
            convKerChaPadded.create(RO, E1, E2);

            hoNDArray<T> kImCha, kImTmp;
            kImCha.create(RO, E1, E2);
            kImTmp.create(RO, E1, E2);

            #pragma omp for 
            for (scha = 0; scha < (long long)srcCHA; scha++)
            {
                hoNDArray<T> unmixCha;
                unmixCha.create(RO, E1, E2, unmixCoeff.begin() + scha*RO*E1*E2);

                for (size_t dcha = 0; dcha < dstCHA; dcha++)
                {
                    convKerCha.create(kRO, kE1, kE2, convKerScaled.begin() + scha*kRO*kE1*kE2 + dcha*kRO*kE1*kE2*srcCHA);
                    Gadgetron::pad(RO, E1, E2, convKerCha, convKerChaPadded, true);
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(convKerChaPadded, kImCha, kImTmp);

                    hoNDArray<T> coilMapCha;
                    coilMapCha.create(RO, E1, E2, const_cast<T*>(coilMap.begin()) + dcha*RO*E1*E2);

                    Gadgetron::multiplyConj(kImCha, coilMapCha, kImTmp);
                    Gadgetron::add(unmixCha, kImTmp, unmixCha);
                }
            }
        }

        hoNDArray<T> conjUnmixCoeff(unmixCoeff);
        Gadgetron::multiplyConj(unmixCoeff, conjUnmixCoeff, conjUnmixCoeff);

        Gadgetron::clear(&gFactor);

        hoNDArray<T> gFactorBuf(RO, E1, E2, 1);
        Gadgetron::sum_over_dimension(conjUnmixCoeff, gFactorBuf, 3);

        complex_to_real(gFactorBuf, gFactor);

        Gadgetron::sqrt(gFactor, gFactor);
        Gadgetron::scal((typename realType<T>::Type)(1.0 / acceFactorE1 / acceFactorE2), gFactor);
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa3d_unmixing_coeff(...) ... ");
    }
}

template EXPORTMRICORE void grappa3d_unmixing_coeff(const hoNDArray< std::complex<float> >& convKer, const hoNDArray< std::complex<float> >& coilMap, size_t acceFactorE1, size_t acceFactorE2, hoNDArray< std::complex<float> >& unmixCoeff, hoNDArray< float >& gFactor);
template EXPORTMRICORE void grappa3d_unmixing_coeff(const hoNDArray< std::complex<double> >& convKer, const hoNDArray< std::complex<double> >& coilMap, size_t acceFactorE1, size_t acceFactorE2, hoNDArray< std::complex<double> >& unmixCoeff, hoNDArray< double >& gFactor);

// ------------------------------------------------------------------------

/// apply grappa convolution kernel to perform per-channel unwrapping
/// convKer: 3D kspace grappa convolution kernel
/// kspace: undersampled kspace [RO E1 E2 srcCHA]
/// unmixCoeff: [RO E1 E2 srcCHA] unmixing coefficient
/// gFactor: [RO E1 E2], gfactor
template <typename T> 
void grappa3d_image_domain_unwrapping(const hoNDArray<T>& convKer, const hoNDArray<T>& kspace,
                                size_t acceFactorE1, size_t acceFactorE2,
                                hoNDArray<T>& complexIm)
{
    try
    {
        size_t kRO = convKer.get_size(0);
        size_t kE1 = convKer.get_size(1);
        size_t kE2 = convKer.get_size(2);
        size_t srcCHA = convKer.get_size(3);
        size_t dstCHA = convKer.get_size(4);

        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);

        GADGET_CHECK_THROW(kspace.get_size(3) == srcCHA);

        size_t N = kspace.get_size(4);

        std::vector<size_t> dim;
        kspace.get_dimensions(dim);

        std::vector<size_t> dimRes(dim);
        dimRes[3] = dstCHA;

        if (!complexIm.dimensions_equal(&dimRes))
        {
            complexIm.create(dim);
        }

        Gadgetron::clear(&complexIm);

        // compute aliased image

        hoNDArray<T> aliasedIm;
        aliasedIm.create(dim);

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kspace, aliasedIm);

        // scaled conv kernel
        hoNDArray<T> convKerScaled;
        convKerScaled = convKer;

        Gadgetron::scal((typename realType<T>::Type)(std::sqrt((double)(RO*E1*E2))), convKerScaled);

        long long dcha;

#pragma omp parallel private(dcha) shared(RO, E1, E2, srcCHA, dstCHA, kRO, kE1, kE2, convKerScaled, aliasedIm, complexIm)
        {
            hoNDArray<T> convKerCha;
            hoNDArray<T> convKerChaPadded;
            convKerChaPadded.create(RO, E1, E2);

            hoNDArray<T> kImCha, kImTmp;
            kImCha.create(RO, E1, E2);
            kImTmp.create(RO, E1, E2);

            hoNDArray<T> aliasedImCha;
            hoNDArray<T> complexImCha;

#pragma omp for 
            for (dcha = 0; dcha < (long long)dstCHA; dcha++)
            {
                for (size_t scha = 0; scha < srcCHA; scha++)
                {
                    convKerCha.create(kRO, kE1, kE2, convKerScaled.begin() + scha*kRO*kE1*kE2 + dcha*kRO*kE1*kE2*srcCHA);
                    Gadgetron::pad(RO, E1, E2, convKerCha, convKerChaPadded, true);
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(convKerChaPadded, kImCha, kImTmp);

                    for (size_t n = 0; n < N; n++)
                    {
                        complexImCha.create(RO, E1, E2, complexIm.begin() + dcha*RO*E1*E2 + n*RO*E1*E2*dstCHA);

                        aliasedImCha.create(RO, E1, E2, aliasedIm.begin() + scha*RO*E1*E2 + n*RO*E1*E2*srcCHA);

                        Gadgetron::multiply(kImCha, aliasedImCha, kImTmp);
                        Gadgetron::add(complexImCha, kImTmp, complexImCha);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa3d_image_domain_unwrapping(...) ... ");
    }
}

template EXPORTMRICORE void grappa3d_image_domain_unwrapping(const hoNDArray< std::complex<float> >& convKer, const hoNDArray< std::complex<float> >& kspace, size_t acceFactorE1, size_t acceFactorE2, hoNDArray< std::complex<float> >& complexIm);
template EXPORTMRICORE void grappa3d_image_domain_unwrapping(const hoNDArray< std::complex<double> >& convKer, const hoNDArray< std::complex<double> >& kspace, size_t acceFactorE1, size_t acceFactorE2, hoNDArray< std::complex<double> >& complexIm);

template <typename T>
void grappa3d_image_domain_unwrapping_aliasedImage(const hoNDArray<T>& convKer, const hoNDArray<T>& aliasedIm,
                                            size_t acceFactorE1, size_t acceFactorE2,
                                            hoNDArray<T>& complexIm)
{
    try
    {
        size_t kRO = convKer.get_size(0);
        size_t kE1 = convKer.get_size(1);
        size_t kE2 = convKer.get_size(2);
        size_t srcCHA = convKer.get_size(3);
        size_t dstCHA = convKer.get_size(4);

        size_t RO = aliasedIm.get_size(0);
        size_t E1 = aliasedIm.get_size(1);
        size_t E2 = aliasedIm.get_size(2);

        GADGET_CHECK_THROW(aliasedIm.get_size(3) == srcCHA);

        size_t N = aliasedIm.get_size(4);

        std::vector<size_t> dim;
        aliasedIm.get_dimensions(dim);

        std::vector<size_t> dimRes(dim);
        dimRes[3] = dstCHA;

        if (!complexIm.dimensions_equal(&dimRes))
        {
            complexIm.create(dim);
        }

        Gadgetron::clear(&complexIm);

        // scaled conv kernel
        hoNDArray<T> convKerScaled;
        convKerScaled = convKer;

        Gadgetron::scal((typename realType<T>::Type)(std::sqrt((double)(RO*E1*E2))), convKerScaled);

        long long dcha;

#pragma omp parallel private(dcha) shared(RO, E1, E2, srcCHA, dstCHA, kRO, kE1, kE2, convKerScaled, aliasedIm, complexIm)
        {
            hoNDArray<T> convKerCha;
            hoNDArray<T> convKerChaPadded;
            convKerChaPadded.create(RO, E1, E2);

            hoNDArray<T> kImCha, kImTmp;
            kImCha.create(RO, E1, E2);
            kImTmp.create(RO, E1, E2);

            hoNDArray<T> aliasedImCha;
            hoNDArray<T> complexImCha;

#pragma omp for 
            for (dcha = 0; dcha < (long long)dstCHA; dcha++)
            {
                for (size_t scha = 0; scha < srcCHA; scha++)
                {
                    convKerCha.create(kRO, kE1, kE2, convKerScaled.begin() + scha*kRO*kE1*kE2 + dcha*kRO*kE1*kE2*srcCHA);
                    Gadgetron::pad(RO, E1, E2, convKerCha, convKerChaPadded, true);
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(convKerChaPadded, kImCha, kImTmp);

                    for (size_t n = 0; n < N; n++)
                    {
                        complexImCha.create(RO, E1, E2, complexIm.begin() + dcha*RO*E1*E2 + n*RO*E1*E2*dstCHA);

                        aliasedImCha.create(RO, E1, E2, const_cast<T*>(aliasedIm.begin()) + scha*RO*E1*E2 + n*RO*E1*E2*srcCHA);

                        Gadgetron::multiply(kImCha, aliasedImCha, kImTmp);
                        Gadgetron::add(complexImCha, kImTmp, complexImCha);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in grappa3d_image_domain_unwrapping_aliasedImage(...) ... ");
    }
}

template EXPORTMRICORE void grappa3d_image_domain_unwrapping_aliasedImage(const hoNDArray< std::complex<float> >& convKer, const hoNDArray< std::complex<float> >& kspace, size_t acceFactorE1, size_t acceFactorE2, hoNDArray< std::complex<float> >& complexIm);
template EXPORTMRICORE void grappa3d_image_domain_unwrapping_aliasedImage(const hoNDArray< std::complex<double> >& convKer, const hoNDArray< std::complex<double> >& kspace, size_t acceFactorE1, size_t acceFactorE2, hoNDArray< std::complex<double> >& complexIm);

// ------------------------------------------------------------------------

template <typename T> 
void apply_unmix_coeff_kspace_3D(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);
        size_t srcCHA = kspace.get_size(3);

        size_t N = kspace.get_size(4);

        std::vector<size_t> dim;
        kspace.get_dimensions(dim);

        GADGET_CHECK_THROW(unmixCoeff.get_size(0) == RO);
        GADGET_CHECK_THROW(unmixCoeff.get_size(1) == E1);
        GADGET_CHECK_THROW(unmixCoeff.get_size(2) == E2);
        GADGET_CHECK_THROW(unmixCoeff.get_size(3) == srcCHA);

        if (complexIm.get_size(0) != RO
            || complexIm.get_size(1) != E1
            || complexIm.get_size(2) != E2
            || complexIm.get_number_of_elements() != RO*E1*E2*N)
        {
            complexIm.create(RO, E1, E2, N);
        }

        hoNDArray<T> aliasedIm, buffer;
        aliasedIm.create(dim);
        buffer.create(dim);
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kspace, aliasedIm, buffer);

        Gadgetron::multiply(aliasedIm, unmixCoeff, buffer);

        hoNDArray<T> bufferIm;
        bufferIm.create(RO, E1, E2, 1, N, complexIm.begin());

        Gadgetron::sum_over_dimension(buffer, bufferIm, 3);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_unmix_coeff_kspace_3D(...) ... ");
    }
}

template EXPORTMRICORE void apply_unmix_coeff_kspace_3D(const hoNDArray< std::complex<float> >& kspace, const hoNDArray< std::complex<float> >& unmixCoeff, hoNDArray< std::complex<float> >& complexIm);
template EXPORTMRICORE void apply_unmix_coeff_kspace_3D(const hoNDArray< std::complex<double> >& kspace, const hoNDArray< std::complex<double> >& unmixCoeff, hoNDArray< std::complex<double> >& complexIm);

// ------------------------------------------------------------------------

template <typename T> 
void apply_unmix_coeff_aliased_image_3D(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        size_t RO = aliasedIm.get_size(0);
        size_t E1 = aliasedIm.get_size(1);
        size_t E2 = aliasedIm.get_size(2);
        size_t srcCHA = aliasedIm.get_size(3);

        size_t N = aliasedIm.get_size(4);

        std::vector<size_t> dim;
        aliasedIm.get_dimensions(dim);

        GADGET_CHECK_THROW(unmixCoeff.get_size(0) == RO);
        GADGET_CHECK_THROW(unmixCoeff.get_size(1) == E1);
        GADGET_CHECK_THROW(unmixCoeff.get_size(2) == E2);
        GADGET_CHECK_THROW(unmixCoeff.get_size(3) == srcCHA);

        if (complexIm.get_size(0) != RO
            || complexIm.get_size(1) != E1
            || complexIm.get_size(2) != E2
            || complexIm.get_number_of_elements() != RO*E1*E2*N)
        {
            complexIm.create(RO, E1, E2, N);
        }

        hoNDArray<T> buffer;
        buffer.create(dim);

        Gadgetron::multiply(aliasedIm, unmixCoeff, buffer);

        hoNDArray<T> bufferIm;
        bufferIm.create(RO, E1, E2, 1, N, complexIm.begin());

        Gadgetron::sum_over_dimension(buffer, bufferIm, 3);
    }
    catch (...)
    {
        GADGET_THROW("Errors in apply_unmix_coeff_aliased_image_3D(...) ... ");
    }
}

template EXPORTMRICORE void apply_unmix_coeff_aliased_image_3D(const hoNDArray< std::complex<float> >& aliasedIm, const hoNDArray< std::complex<float> >& unmixCoeff, hoNDArray< std::complex<float> >& complexIm);
template EXPORTMRICORE void apply_unmix_coeff_aliased_image_3D(const hoNDArray< std::complex<double> >& aliasedIm, const hoNDArray< std::complex<double> >& unmixCoeff, hoNDArray< std::complex<double> >& complexIm);

// ------------------------------------------------------------------------

}
