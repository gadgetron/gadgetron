
/** \file   mri_core_partial_fourier.cpp
    \brief  Implementation partial fourier handling functionalities for 2D and 3D MRI
    \author Hui Xue
*/

#include "mri_core_partial_fourier.h"
#include "mri_core_kspace_filter.h"
#include "hoNDFFT.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_reductions.h"
#include "ho2DArray.h"
#include "ho3DArray.h"
#include "hoMatrix.h"
#include "hoNDArray_utils.h"

namespace Gadgetron
{
    // ------------------------------------------------------------------------

    std::string get_ismrmrd_partial_fourier_handling_algo_name(ismrmrdPARTIALFOURIERHANDLINGALGO v)
    {
        std::string name;

        switch (v)
        {
        case ISMRMRD_PF_None:
            name = "None";
            break;

        case ISMRMRD_PF_Zerofilling_filter:
            name = "ZeroFillingFilter";
            break;

        case ISMRMRD_PF_Pocs:
            name = "POCS";
            break;

        case ISMRMRD_PF_FengHuang:
            name = "FengHuang";
            break;

        default:
            GERROR_STREAM("Unrecognized ISMRMRD partial fourier handling algorithm type : " << v);
        }

        return name;
    }

    ismrmrdPARTIALFOURIERHANDLINGALGO get_ismrmrd_partial_fourier_handling_algo(const std::string& name)
    {
        ismrmrdPARTIALFOURIERHANDLINGALGO v;

        if (name == "None")
        {
            v = ISMRMRD_PF_None;
        }
        else if (name == "ZeroFillingFilter")
        {
            v = ISMRMRD_PF_Zerofilling_filter;
        }
        else if (name == "POCS")
        {
            v = ISMRMRD_PF_Pocs;
        }
        else if (name == "FengHuang")
        {
            v = ISMRMRD_PF_FengHuang;
        }
        else
        {
            GERROR_STREAM("Unrecognized ISMRMRD partial fourier handling algorithm name : " << name);
        }

        return v;
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_reset_kspace_2d(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1)
    {
        try
        {
            size_t NDim = src.get_number_of_dimensions();
            GADGET_CHECK_THROW(NDim >= 2);

            size_t RO = dst.get_size(0);
            size_t E1 = dst.get_size(1);

            size_t RO_src = src.get_size(0);
            size_t E1_src = src.get_size(1);

            GADGET_CHECK_THROW(RO == RO_src);
            GADGET_CHECK_THROW(E1 == E1_src);
            GADGET_CHECK_THROW(src.get_number_of_elements() == dst.get_number_of_elements());

            if ((startRO >= RO) || (endRO >= RO) || (startRO>endRO))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_reset_kspace_2d(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
                return;
            }

            if ((startE1 >= E1) || (endE1 >= E1) || (startE1>endE1))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_reset_kspace_2d(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
                return;
            }

            size_t N = dst.get_number_of_elements() / (RO*E1);
            const T* pSrc = src.begin();
            T* pDst = dst.begin();

            long long n;

#pragma omp parallel for default(none) private(n) shared(N, pSrc, pDst, RO, E1, startRO, endRO, startE1, endE1)
            for (n = 0; n<(long long)N; n++)
            {
                for (size_t e1 = startE1; e1 <= endE1; e1++)
                {
                    size_t offset = n*RO*E1 + e1*RO + startRO;
                    memcpy(pDst + offset, pSrc + offset, sizeof(T)*(endRO - startRO + 1));
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in partial_fourier_reset_kspace_2d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_reset_kspace_3d(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2)
    {
        try
        {
            size_t NDim = src.get_number_of_dimensions();
            GADGET_CHECK_THROW(NDim >= 2);

            size_t RO = dst.get_size(0);
            size_t E1 = dst.get_size(1);
            size_t E2 = dst.get_size(2);

            size_t RO_src = src.get_size(0);
            size_t E1_src = src.get_size(1);
            size_t E2_src = src.get_size(2);

            GADGET_CHECK_THROW(RO == RO_src);
            GADGET_CHECK_THROW(E1 == E1_src);
            GADGET_CHECK_THROW(E2 == E2_src);
            GADGET_CHECK_THROW(src.get_number_of_elements() == dst.get_number_of_elements());

            if ((startRO >= RO) || (endRO >= RO) || (startRO>endRO))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_reset_kspace_3d(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
                return;
            }

            if ((startE1 >= E1) || (endE1 >= E1) || (startE1>endE1))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_reset_kspace_3d(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
                return;
            }

            if ((startE2 >= E2) || (endE2 >= E2) || (startE2>endE2))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_reset_kspace_3d(...) : (startE2>=E2) || (endE2>=E2) || (startE2>endE2) ... ");
                return;
            }

            size_t N = dst.get_number_of_elements() / (RO*E1*E2);
            const T* pSrc = src.begin();
            T* pDst = dst.begin();

            long long n;

#pragma omp parallel for default(none) private(n) shared(N, pSrc, pDst, RO, E1, E2, startRO, endRO, startE1, endE1, startE2, endE2)
            for (n = 0; n<(long long)N; n++)
            {
                for (size_t e2 = startE2; e2 <= endE2; e2++)
                {
                    for (size_t e1 = startE1; e1 <= endE1; e1++)
                    {
                        size_t offset = n*RO*E1*E2 + e2*E1*RO + e1*RO + startRO;
                        memcpy(pDst + offset, pSrc + offset, sizeof(T)*(endRO - startRO + 1));
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in partial_fourier_reset_kspace_3d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_transition_band_2d(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t transBandRO, size_t transBandE1)
    {
        try
        {
            size_t NDim = src.get_number_of_dimensions();
            GADGET_CHECK_THROW(NDim >= 2);

            size_t RO = dst.get_size(0);
            size_t E1 = dst.get_size(1);

            size_t RO_src = src.get_size(0);
            size_t E1_src = src.get_size(1);

            GADGET_CHECK_THROW(RO == RO_src);
            GADGET_CHECK_THROW(E1 == E1_src);
            GADGET_CHECK_THROW(src.get_number_of_elements() == dst.get_number_of_elements());

            if ((startRO >= RO) || (endRO >= RO) || (startRO>endRO))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_transition_band_2d(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
                return;
            }

            if ((startE1 >= E1) || (endE1 >= E1) || (startE1>endE1))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_transition_band_2d(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
                return;
            }

            while (transBandRO>1 && startRO + transBandRO > RO / 2)
            {
                transBandRO--;
            }

            while (transBandRO>1 && endRO - transBandRO < RO / 2)
            {
                transBandRO--;
            }

            while (transBandE1>1 && startE1 + transBandE1 > E1 / 2)
            {
                transBandE1--;
            }

            while (transBandE1>1 && endE1 - transBandE1 < E1 / 2)
            {
                transBandE1--;
            }

            ISMRMRDKSPACEFILTER filterType = ISMRMRD_FILTER_TAPERED_HANNING;
            bool densityComp = false;

            hoNDArray<T> filter_src_RO, filter_src_E1;

            if (startRO == 0 && endRO == RO - 1)
            {
                Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_NONE, transBandRO, densityComp);
            }
            else
            {
                Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_TAPERED_HANNING, transBandRO, densityComp);
            }

            if (startE1 == 0 && endE1 == E1 - 1)
            {
                Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_NONE, transBandE1, densityComp);
            }
            else
            {
                Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_TAPERED_HANNING, transBandE1, densityComp);
            }

            // in this way, the SNR unit scale property is perserved
            T midValue = filter_src_RO(RO / 2);
            T scalFactor = T(1.0) / midValue;
            Gadgetron::scal(scalFactor, filter_src_RO);

            midValue = filter_src_E1(E1 / 2);
            scalFactor = T(1.0) / midValue;
            Gadgetron::scal(scalFactor, filter_src_E1);

            hoNDArray<T> filter_dst_RO(RO), filter_dst_E1(E1);

            size_t ii;
            for (ii = 0; ii<RO; ii++)
            {
                filter_dst_RO(ii) = T(1.0) - filter_src_RO(ii);
            }

            for (ii = 0; ii<E1; ii++)
            {
                filter_dst_E1(ii) = T(1.0) - filter_src_E1(ii);
            }

            hoNDArray<T> srcFiltered(src), dstFiltered(dst);
            if (startRO == 0 && endRO == RO - 1)
            {
                Gadgetron::apply_kspace_filter_E1(src, filter_src_E1, srcFiltered);
                Gadgetron::apply_kspace_filter_E1(dst, filter_dst_E1, dstFiltered);
            }
            else if (startE1 == 0 && endE1 == E1 - 1)
            {
                Gadgetron::apply_kspace_filter_RO(src, filter_src_RO, srcFiltered);
                Gadgetron::apply_kspace_filter_RO(dst, filter_dst_RO, dstFiltered);
            }
            else
            {
                Gadgetron::apply_kspace_filter_ROE1(src, filter_src_RO, filter_src_E1, srcFiltered);

                hoNDArray<T> fxy;
                Gadgetron::compute_2d_filter(filter_src_RO, filter_src_E1, fxy);

                size_t Nxy = RO*E1;
                for (ii = 0; ii<Nxy; ii++)
                {
                    fxy(ii) = T(1.0) - fxy(ii);
                }

                Gadgetron::apply_kspace_filter_ROE1(dst, fxy, dstFiltered);
            }

            Gadgetron::add(srcFiltered, dstFiltered, dst);
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in partial_fourier_transition_band_2d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_transition_band_3d(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, size_t transBandRO, size_t transBandE1, size_t transBandE2)
    {
        try
        {
            size_t NDim = src.get_number_of_dimensions();
            GADGET_CHECK_THROW(NDim >= 3);

            size_t RO = dst.get_size(0);
            size_t E1 = dst.get_size(1);
            size_t E2 = dst.get_size(2);

            size_t RO_src = src.get_size(0);
            size_t E1_src = src.get_size(1);
            size_t E2_src = src.get_size(2);

            GADGET_CHECK_THROW(RO == RO_src);
            GADGET_CHECK_THROW(E1 == E1_src);
            GADGET_CHECK_THROW(E2 == E2_src);
            GADGET_CHECK_THROW(src.get_number_of_elements() == dst.get_number_of_elements());

            if ((startRO >= RO) || (endRO >= RO) || (startRO>endRO))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_transition_band_3d(...) : (startRO>=RO) || (endRO>=RO) || (startRO>endRO) ... ");
                return;
            }

            if ((startE1 >= E1) || (endE1 >= E1) || (startE1>endE1))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_transition_band_3d(...) : (startE1>=E1) || (endE1>=E1) || (startE1>endE1) ... ");
                return;
            }

            if ((startE2 >= E2) || (endE2 >= E2) || (startE2>endE2))
            {
                dst = src;
                GWARN_STREAM("partial_fourier_transition_band_3d(...) : (startE2>=E2) || (endE2>=E2) || (startE2>endE2) ... ");
                return;
            }

            while (transBandRO>1 && startRO + transBandRO > RO / 2)
            {
                transBandRO--;
            }

            while (transBandRO>1 && endRO - transBandRO < RO / 2)
            {
                transBandRO--;
            }

            while (transBandE1>1 && startE1 + transBandE1 > E1 / 2)
            {
                transBandE1--;
            }

            while (transBandE1>1 && endE1 - transBandE1 < E1 / 2)
            {
                transBandE1--;
            }

            while (transBandE2>1 && startE2 + transBandE2 > E2 / 2)
            {
                transBandE2--;
            }

            while (transBandE2>1 && endE2 - transBandE2 < E2 / 2)
            {
                transBandE2--;
            }

            ISMRMRDKSPACEFILTER filterType = ISMRMRD_FILTER_TAPERED_HANNING;
            bool densityComp = false;

            hoNDArray<T> filter_src_RO, filter_src_E1, filter_src_E2;

            if (startRO == 0 && endRO == RO - 1)
            {
                Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_NONE, transBandRO, densityComp);
            }
            else
            {
                Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_TAPERED_HANNING, transBandRO, densityComp);
            }

            if (startE1 == 0 && endE1 == E1 - 1)
            {
                Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_NONE, transBandE1, densityComp);
            }
            else
            {
                Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_TAPERED_HANNING, transBandE1, densityComp);
            }

            if (startE2 == 0 && endE2 == E2 - 1)
            {
                Gadgetron::generate_asymmetric_filter(E2, startE2, endE2, filter_src_E2, ISMRMRD_FILTER_NONE, transBandE2, densityComp);
            }
            else
            {
                Gadgetron::generate_asymmetric_filter(E2, startE2, endE2, filter_src_E2, ISMRMRD_FILTER_TAPERED_HANNING, transBandE2, densityComp);
            }

            // in this way, the SNR unit scale property is perserved
            T midValue = filter_src_RO(RO / 2);
            T scalFactor = T(1.0) / midValue;
            Gadgetron::scal(scalFactor, filter_src_RO);

            midValue = filter_src_E1(E1 / 2);
            scalFactor = T(1.0) / midValue;
            Gadgetron::scal(scalFactor, filter_src_E1);

            midValue = filter_src_E2(E2 / 2);
            scalFactor = T(1.0) / midValue;
            Gadgetron::scal(scalFactor, filter_src_E2);

            hoNDArray<T> filter_dst_RO(RO), filter_dst_E1(E1), filter_dst_E2(E2);

            size_t ii;
            for (ii = 0; ii<RO; ii++)
            {
                filter_dst_RO(ii) = T(1.0) - filter_src_RO(ii);
            }

            for (ii = 0; ii<E1; ii++)
            {
                filter_dst_E1(ii) = T(1.0) - filter_src_E1(ii);
            }

            for (ii = 0; ii<E2; ii++)
            {
                filter_dst_E2(ii) = T(1.0) - filter_src_E2(ii);
            }

            hoNDArray<T> srcFiltered(src), dstFiltered(dst);
            if (endRO <= RO - 1 && startE1 == 0 && endE1 == E1 - 1 && startE2 == 0 && endE1 == E2 - 1)
            {
                Gadgetron::apply_kspace_filter_RO(src, filter_src_E1, srcFiltered);
                Gadgetron::apply_kspace_filter_RO(dst, filter_dst_E1, dstFiltered);
            }
            else if (startRO == 0 && endRO == RO - 1 && endE1 <= E1 - 1 && startE2 == 0 && endE1 == E2 - 1)
            {
                Gadgetron::apply_kspace_filter_E1(src, filter_src_RO, srcFiltered);
                Gadgetron::apply_kspace_filter_E1(dst, filter_dst_RO, dstFiltered);
            }
            else if (startRO == 0 && endRO == RO - 1 && startE1 == 0 && endE1 == E1 - 1 && endE1 <= E2 - 1)
            {
                Gadgetron::apply_kspace_filter_E2(src, filter_src_RO, srcFiltered);
                Gadgetron::apply_kspace_filter_E2(dst, filter_dst_RO, dstFiltered);
            }
            else
            {
                Gadgetron::apply_kspace_filter_ROE1E2(src, filter_src_RO, filter_src_E1, filter_src_E2, srcFiltered);

                hoNDArray<T> fxyz;
                Gadgetron::compute_3d_filter(filter_src_RO, filter_src_E1, filter_src_E2, fxyz);

                size_t Nxyz = RO*E1*E2;
                for (ii = 0; ii<Nxyz; ii++)
                {
                    fxyz(ii) = T(1.0) - fxyz(ii);
                }

                Gadgetron::apply_kspace_filter_ROE1E2(dst, fxyz, dstFiltered);
            }

            Gadgetron::add(srcFiltered, dstFiltered, dst);
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in partial_fourier_transition_band_3d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_conjugate_symmetry_2d(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj)
    {
        try
        {
            if (!kspaceConj.dimensions_equal(&kspace))
            {
                kspaceConj.create(kspace.get_dimensions());
            }

            long long RO = kspace.get_size(0);
            long long E1 = kspace.get_size(1);
            long long num = kspace.get_number_of_elements() / (RO*E1);

            long long centerRO = RO / 2;
            long long centerE1 = E1 / 2;

            long long ii;

#pragma omp parallel for default(none) private(ii) shared(RO, E1, num, centerRO, centerE1, kspace, kspaceConj)
            for (ii = 0; ii<num; ii++)
            {
                ho2DArray<T> src(RO, E1, const_cast<T*>(kspace.begin() + ii*RO*E1));
                ho2DArray<T> dst(RO, E1, const_cast<T*>(kspaceConj.begin() + ii*RO*E1));

                long long ro, e1;
                long long cro, ce1;

                for (e1 = 0; e1<E1; e1++)
                {
                    ce1 = 2 * centerE1 - e1;
                    if (ce1 > E1 - 1)
                    {
                        ce1 -= E1;
                    }

                    for (ro = 0; ro<RO; ro++)
                    {
                        cro = 2 * centerRO - ro;
                        if (cro > RO - 1)
                        {
                            cro -= RO;
                        }

                        dst(ro, e1) = std::conj(src(cro, ce1));
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in partial_fourier_conjugate_symmetry_2d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_conjugate_symmetry_3d(const hoNDArray<T>& kspace, hoNDArray<T>& kspaceConj)
    {
        try
        {
            if (!kspaceConj.dimensions_equal(&kspace))
            {
                kspaceConj.create(kspace.get_dimensions());
            }

            long long RO = kspace.get_size(0);
            long long E1 = kspace.get_size(1);
            long long E2 = kspace.get_size(2);
            long long num = kspace.get_number_of_elements() / (RO*E1*E2);

            long long centerRO = RO / 2;
            long long centerE1 = E1 / 2;
            long long centerE2 = E2 / 2;

            long long ii;

#pragma omp parallel for default(none) private(ii) shared(RO, E1, E2, num, centerRO, centerE1, centerE2, kspace, kspaceConj)
            for (ii = 0; ii<num; ii++)
            {
                ho3DArray<T> src(RO, E1, E2, const_cast<T*>(kspace.begin() + ii*RO*E1*E2));
                ho3DArray<T> dst(RO, E1, E2, const_cast<T*>(kspaceConj.begin() + ii*RO*E1*E2));

                long long ro, e1, e2;
                long long cro, ce1, ce2;

                for (e2 = 0; e2<E2; e2++)
                {
                    ce2 = 2 * centerE2 - e2;
                    if (ce2 > E2 - 1)
                    {
                        ce2 -= E2;
                    }

                    for (e1 = 0; e1<E1; e1++)
                    {
                        ce1 = 2 * centerE1 - e1;
                        if (ce1 > E1 - 1)
                        {
                            ce1 -= E1;
                        }

                        for (ro = 0; ro<RO; ro++)
                        {
                            cro = 2 * centerRO - ro;
                            if (cro > RO - 1)
                            {
                                cro -= RO;
                            }

                            dst(ro, e1, e2) = std::conj(src(cro, ce1, ce2));
                        }
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in partial_fourier_conjugate_symmetry_3d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    /// kspace: input kspae [RO E1 E2 ...]
    /// if is3D==false, 2D POCS is performed, otherwise 3D POCS is performed
    /// startRO, endRO, startE1, endE1, startE2, endE2: acquired kspace range
    /// transit_band_RO/E1/E2: transition band width in pixel for RO/E1/E2
    /// iter: number of maximal iterations for POCS
    /// thres: iteration threshold
    template <typename T> 
    void partial_fourier_POCS(const hoNDArray<T>& kspace,
                            size_t startRO, size_t endRO, 
                            size_t startE1, size_t endE1, 
                            size_t startE2, size_t endE2,
                            size_t transit_band_RO, 
                            size_t transit_band_E1, 
                            size_t transit_band_E2,
                            size_t iter, double thres, 
                            bool is3D, 
                            hoNDArray<T>& res,
                            bool verbose)
    {
        try
        {
            size_t RO = kspace.get_size(0);
            size_t E1 = kspace.get_size(1);
            size_t E2 = kspace.get_size(2);

            res = kspace;

            // check whether partial fourier is used
            if (startRO >= RO || endRO >= RO || startRO >= endRO)
            {
                GWARN_STREAM("partial_fourier_POCS(...) : (startRO >= RO || endRO >= RO || startRO >= endRO) ... ");
                startRO = 0;
                endRO = RO - 1;
            }

            if (startE1 >= E1 || endE1 >= E1 || startE1 >= endE1)
            {
                GWARN_STREAM("partial_fourier_POCS(...) : (startE1 >= E1 || endE1 >= E1 || startE1 >= endE1) ... ");
                startE1 = 0;
                endE1 = E1 - 1;
            }

            if (is3D)
            {
                if (startE2 >= E2 || endE2 >= E2 || startE2 >= endE2)
                {
                    GWARN_STREAM("partial_fourier_POCS(...) : (startE2 >= E2 || endE2 >= E2 || startE2 >= endE2) ... ");
                    startE2 = 0;
                    endE2 = E2 - 1;
                }
            }

            if ((endRO - startRO + 1 == RO) && (endE1 - startE1 + 1 == E1) && (endE2 - startE2 + 1 == E2))
            {
                GWARN_STREAM("partial_fourier_POCS(...) : (endRO - startRO + 1 == RO) && (endE1 - startE1 + 1 == E1) && (endE2 - startE2 + 1 == E2) ... ");
                return;
            }

            // create kspace filter for homodyne phase estimation
            hoNDArray<T> filterRO(RO);
            Gadgetron::generate_symmetric_filter_ref(RO, startRO, endRO, filterRO);

            hoNDArray<T> filterE1(E1);
            Gadgetron::generate_symmetric_filter_ref(E1, startE1, endE1, filterE1);

            hoNDArray<T> filterE2(E2);
            if (is3D)
            {
                Gadgetron::generate_symmetric_filter_ref(E2, startE2, endE2, filterE2);
            }

            hoNDArray<T> kspaceIter(kspace);
            // magnitude of complex images
            hoNDArray<typename realType<T>::Type> mag(kspace.get_dimensions());
            hoNDArray<T> magComplex(kspace.get_dimensions());

            // kspace filter
            hoNDArray<T> buffer_partial_fourier(kspaceIter), buffer(kspaceIter);

            if (is3D)
            {
                Gadgetron::apply_kspace_filter_ROE1E2(kspaceIter, filterRO, filterE1, filterE2, buffer_partial_fourier);

                // go to image domain
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(buffer_partial_fourier);
            }
            else
            {
                Gadgetron::apply_kspace_filter_ROE1(kspaceIter, filterRO, filterE1, buffer_partial_fourier);

                // go to image domain
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(buffer_partial_fourier);
            }

            // get the complex image phase for the filtered kspace
            Gadgetron::abs(buffer_partial_fourier, mag);
            Gadgetron::addEpsilon(mag);
            magComplex.copyFrom(mag);
            Gadgetron::divide(buffer_partial_fourier, magComplex, buffer);

            // complex images, initialized as not filtered complex image
            hoNDArray<T> complexIm(kspaceIter);

            if (is3D)
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kspaceIter, complexIm);
            }
            else
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspaceIter, complexIm);
            }

            hoNDArray<T> complexImPOCS(complexIm);

            // the kspace during iteration is buffered here
            hoNDArray<T> buffer_partial_fourier_Iter(kspaceIter);

            size_t ii;
            for (ii = 0; ii<iter; ii++)
            {
                Gadgetron::abs(complexImPOCS, mag);
                magComplex.copyFrom(mag);
                Gadgetron::multiply(magComplex, buffer, complexImPOCS);

                // go back to kspace
                if (is3D)
                {
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(complexImPOCS, kspaceIter);
                }
                else
                {
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(complexImPOCS, kspaceIter);
                }

                // buffer kspace during iteration
                buffer_partial_fourier_Iter = kspaceIter;

                // restore the acquired region
                if (is3D)
                {
                    partial_fourier_reset_kspace_3d(kspace, kspaceIter, startRO, endRO, startE1, endE1, startE2, endE2);

                    // update complex image
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kspaceIter, complexImPOCS);
                }
                else
                {
                    partial_fourier_reset_kspace_2d(kspace, kspaceIter, startRO, endRO, startE1, endE1);

                    // update complex image
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspaceIter, complexImPOCS);
                }

                // compute threshold to stop the iteration
                Gadgetron::subtract(complexImPOCS, complexIm, buffer_partial_fourier);
                typename realType<T>::Type diff, prev;
                Gadgetron::norm2(complexIm, prev);
                Gadgetron::norm2(buffer_partial_fourier, diff);

                typename realType<T>::Type t = diff / prev;

                if (verbose)
                {
                    GDEBUG_STREAM("POCS iter : " << ii << " - thres : " << t << " ... ");
                }

                if (t < thres)
                {
                    break;
                }

                complexIm = complexImPOCS;
            }

            if (transit_band_RO == 0 && transit_band_E1 == 0 && transit_band_E2==0)
            {
                res = kspaceIter;
            }
            else
            {
                if (is3D)
                {
                    Gadgetron::partial_fourier_transition_band_3d(kspace, buffer_partial_fourier_Iter, startRO, endRO, startE1, endE1, startE2, endE2, transit_band_RO, transit_band_E1, transit_band_E2);
                }
                else
                {
                    Gadgetron::partial_fourier_transition_band_2d(kspace, buffer_partial_fourier_Iter, startRO, endRO, startE1, endE1, transit_band_RO, transit_band_E1);
                }

                res = buffer_partial_fourier_Iter;
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in partial_fourier_POCS(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_FengHuang_calib_2d(const hoNDArray<T>& src, const hoNDArray<T>& dst, size_t kRO, size_t kE1, double thres, hoNDArray<T>& kernel)
    {
        try
        {
            GADGET_CHECK_THROW(src.dimensions_equal(&dst));

            long long RO = (long long)src.get_size(0);
            long long E1 = (long long)src.get_size(1);

            long long kx = (long long)kRO;
            long long ky = (long long)kE1;

            if (kx % 2 == 0) kx++;
            if (ky % 2 == 0) ky++;

            long long halfKx = (long long)kx / 2;
            long long halfKy = (long long)ky / 2;

            /// the cross-channel kernel is not estimated
            std::vector<size_t> dim;
            src.get_dimensions(dim);

            dim[0] = kx;
            dim[1] = ky;

            kernel.create(dim);

            long long ii = 0;
            long long num = kernel.get_number_of_elements() / (kx*ky);

            size_t startRO = halfKx;
            size_t endRO = RO - halfKx - 1;

            size_t startE1 = halfKy;
            size_t endE1 = E1 - halfKy - 1;

            long long rowA, colA, rowB, colB;
            rowA = (endE1 - startE1 + 1)*(endRO - startRO + 1);
            colA = kx*ky;

            rowB = rowA;
            colB = 1;

#pragma omp parallel default(none) private(ii) shared(num, RO, E1, kx, ky, src, dst, kernel, rowA, colA, rowB, colB, startRO, endRO, startE1, endE1, halfKx, halfKy, thres)
            {
                hoMatrix<T> A(rowA, colA);
                T* pA = A.begin();

                hoMatrix<T> B(rowB, colB);
                T* pB = B.begin();

                hoMatrix<T> K(colA, colB);

#pragma omp for
                for (ii = 0; ii<num; ii++)
                {
                    T* pSrc2D = const_cast<T*>(src.begin()) + ii*RO*E1;
                    T* pDst2D = const_cast<T*>(dst.begin()) + ii*RO*E1;

                    size_t ro, e1, row(0);
                    long long x, y;

                    for (e1 = startE1; e1 <= endE1; e1++)
                    {
                        for (ro = startRO; ro <= endRO; ro++)
                        {

                            size_t colInd(0);
                            for (y = -halfKy; y <= halfKy; y++)
                            {
                                for (x = -halfKx; x <= halfKx; x++)
                                {
                                    pA[row + colInd*rowA] = pSrc2D[ro + x + (e1 + y)*RO];
                                    colInd++;
                                }
                            }

                            pB[row] = pDst2D[ro + e1*RO];

                            row++;
                        }
                    }

                    Gadgetron::SolveLinearSystem_Tikhonov(A, B, K, thres);

                    memcpy(kernel.begin() + ii*kx*ky, K.begin(), sizeof(T)*kx*ky);
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in partial_fourier_FengHuang_calib_2d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_FengHuang_apply_kernel_2d(const hoNDArray<T>& kspaceConj, size_t startRO, size_t endRO, size_t startE1, size_t endE1, const hoNDArray<T>& kernel, hoNDArray<T>& kspace)
    {
        try
        {
            GADGET_CHECK_THROW(kspaceConj.dimensions_equal(&kspace));

            long long RO = (long long)kspace.get_size(0);
            long long E1 = (long long)kspace.get_size(1);

            long long kx = (long long)kernel.get_size(0);
            long long ky = (long long)kernel.get_size(1);

            long long halfKx = kx / 2;
            long long halfKy = ky / 2;

            long long num = kspace.get_number_of_elements() / (RO*E1);

            long long rowD = RO*E1 - ((endE1 - startE1 + 1) * (endRO - startRO + 1));
            long long colD = kx*ky;

            ho2DArray<size_t> coeffX(rowD, colD);
            ho2DArray<size_t> coeffY(rowD, colD);

            long long ro, e1, row(0);
            long long x, y, dx, dy;

            for (e1 = 0; e1<E1; e1++)
            {
                for (ro = 0; ro<RO; ro++)
                {
                    if ((ro >= (long long)startRO) && (ro <= (long long)endRO) && (e1 >= (long long)startE1) && (e1 <= (long long)endE1))
                    {
                        continue;
                    }

                    size_t colInd(0);
                    for (y = -halfKy; y <= halfKy; y++)
                    {
                        dy = e1 + y;
                        if (dy < 0) dy += E1;
                        if (dy > E1 - 1) dy -= E1;

                        for (x = -halfKx; x <= halfKx; x++)
                        {
                            dx = ro + x;
                            if (dx < 0) dx += RO;
                            if (dx > RO - 1) dx -= RO;

                            coeffX(row, colInd) = dx;
                            coeffY(row, colInd) = dy;
                            colInd++;
                        }
                    }

                    row++;
                }
            }

            long long ii;
#pragma omp parallel default(none) private(ii) shared(num, RO, E1, kspaceConj, kspace, kernel, rowD, colD, coeffX, coeffY)
            {
                hoMatrix<T> D(rowD, colD);
                hoMatrix<T> K(colD, 1);
                hoMatrix<T> R(rowD, 1);

                Gadgetron::clear(D);
                Gadgetron::clear(K);
                Gadgetron::clear(R);

#pragma omp for
                for (ii = 0; ii<num; ii++)
                {
                    ho2DArray<T> src2D(RO, E1, const_cast<T*>(kspaceConj.begin()) + ii*RO*E1);
                    ho2DArray<T> dst2D(RO, E1, kspace.begin() + ii*RO*E1);

                    long long row, col;
                    for (col = 0; col<colD; col++)
                    {
                        for (row = 0; row<rowD; row++)
                        {
                            D(row, col) = src2D(coeffX(row, col), coeffY(row, col));
                        }
                    }

                    memcpy(K.begin(), kernel.begin() + ii*colD, sizeof(T)*colD);

                    // R = D*K
                    Gadgetron::gemm(R, D, false, K, false);

                    for (row = 0; row<rowD; row++)
                    {
                        dst2D(coeffX(row, colD / 2), coeffY(row, colD / 2)) = R(row, 0);
                    }
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in partial_fourier_FengHuang_apply_kernel_2d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T> 
    void partial_fourier_FengHuang_2d(const hoNDArray<T>& kspace,
                                    size_t startRO, size_t endRO,
                                    size_t startE1, size_t endE1,
                                    size_t transit_band_RO,
                                    size_t transit_band_E1,
                                    size_t kRO, size_t kE1,
                                    double thres,
                                    hoNDArray<T>& res)
    {
        try
        {
            size_t RO = kspace.get_size(0);
            size_t E1 = kspace.get_size(1);

            res = kspace;

            if (startRO >= RO || endRO >= RO || startRO >= endRO)
            {
                GWARN_STREAM("partial_fourier_FengHuang_2d(...) : (startRO >= RO || endRO >= RO || startRO >= endRO) ... ");
                startRO = 0;
                endRO = RO - 1;
            }

            if (startE1 >= E1 || endE1 >= E1 || startE1 >= endE1)
            {
                GWARN_STREAM("partial_fourier_FengHuang_2d(...) : (startE1 >= E1 || endE1 >= E1 || startE1 >= endE1) ... ");
                startE1 = 0;
                endE1 = E1 - 1;
            }

            if ((endRO - startRO + 1 == RO) && (endE1 - startE1 + 1 == E1) )
            {
                GWARN_STREAM("partial_fourier_FengHuang_2d(...) : (endRO - startRO + 1 == RO) && (endE1 - startE1 + 1 == E1) ... ");
                return;
            }

            // compute the conjugate symmetric kspace
            hoNDArray<T> buffer2DT_;
            Gadgetron::partial_fourier_conjugate_symmetry_2d(kspace, buffer2DT_);

            // find the symmetric region in the kspace
            size_t startSymRO, endSymRO;
            Gadgetron::find_symmetric_sampled_region(startRO, endRO, RO / 2, startSymRO, endSymRO);

            size_t startSymE1, endSymE1;
            Gadgetron::find_symmetric_sampled_region(startE1, endE1, E1 / 2, startSymE1, endSymE1);

            // the reference kspace for kernel estimation
            hoNDArray<T> src, dst;
            Gadgetron::vector_td<size_t, 2> start, size;

            start[0] = startSymRO;
            start[1] = startSymE1;

            size[0] = endSymRO - startSymRO + 1;
            size[1] = endSymE1 - startSymE1 + 1;

            Gadgetron::crop(start, size, &buffer2DT_, &src);
            Gadgetron::crop(start, size, const_cast< hoNDArray<T>* >(&kspace), &dst);

            // estimate the kernels
            hoNDArray<T> kernel;
            Gadgetron::partial_fourier_FengHuang_calib_2d(src, dst, kRO, kE1, thres, kernel);

            // perform the recon
            if (transit_band_RO == 0 && transit_band_E1 == 0)
            {
                partial_fourier_FengHuang_apply_kernel_2d(buffer2DT_, startRO, endRO, startE1, endE1, kernel, res);
            }
            else
            {
                size_t sRO(startRO), eRO(endRO), sE1(startE1), eE1(endE1);

                if (startRO > 0)
                {
                    startRO += transit_band_RO;
                    if (startRO > RO) startRO = 0;
                }

                if (endRO < RO - 1)
                {
                    endRO -= transit_band_RO;
                }

                if (startRO > endRO)
                {
                    startRO = 0;
                    endRO = RO - 1;
                }

                if (startE1 > 0)
                {
                    startE1 += transit_band_E1;
                    if (startE1 > E1) startE1 = 0;
                }

                if (endE1 < E1 - 1)
                {
                    endE1 -= transit_band_E1;
                }

                if (startE1 > endE1)
                {
                    startE1 = 0;
                    endE1 = E1 - 1;
                }

                Gadgetron::partial_fourier_FengHuang_apply_kernel_2d(buffer2DT_, startRO, endRO, startE1, endE1, kernel, res);

                Gadgetron::partial_fourier_transition_band_2d(kspace, res, sRO, eRO, sE1, eE1, transit_band_RO, transit_band_E1);
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in partial_fourier_FengHuang_2d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_FengHuang_calib_3d(const hoNDArray<T>& src, const hoNDArray<T>& dst, size_t kRO, size_t kE1, size_t kE2, double thres, hoNDArray<T>& kernel)
    {
        try
        {
            GADGET_CHECK_THROW(src.dimensions_equal(&dst));

            long long RO = (long long)src.get_size(0);
            long long E1 = (long long)src.get_size(1);
            long long E2 = (long long)src.get_size(2);

            long long kx = (long long)kRO;
            long long ky = (long long)kE1;
            long long kz = (long long)kE2;

            if (kx % 2 == 0) kx++;
            if (ky % 2 == 0) ky++;
            if (kz % 2 == 0) kz++;

            long long halfKx = (long long)kx / 2;
            long long halfKy = (long long)ky / 2;
            long long halfKz = (long long)kz / 2;

            // the cross-channel kernel is not estimated
            std::vector<size_t> dim;
            src.get_dimensions(dim);

            dim[0] = kx;
            dim[1] = ky;
            dim[2] = kz;

            kernel.create(dim);

            long long ii = 0;
            long long num = kernel.get_number_of_elements()/(kx*ky*kz);

            long long startRO = halfKx;
            long long endRO = RO - halfKx - 1;

            long long startE1 = halfKy;
            long long endE1 = E1 - halfKy - 1;

            long long startE2 = halfKz;
            long long endE2 = E2 - halfKz - 1;

            long long rowA, colA, rowB, colB;
            rowA = (endE2 - startE2 + 1)*(endE1 - startE1 + 1)*(endRO - startRO + 1);
            colA = kx*ky*kz;

            rowB = rowA;
            colB = 1;

#pragma omp parallel default(none) private(ii) shared(num, RO, E1, E2, kx, ky, kz, src, dst, kernel, rowA, colA, rowB, colB, startRO, endRO, startE1, endE1, startE2, endE2, halfKx, halfKy, halfKz, thres) if ( num > 1 ) num_threads( (int)(num<16 ? num : 16) )
            {
                hoNDArray<T> A_mem(rowA, colA);
                hoNDArray<T> B_mem(rowB, colB);
                hoNDArray<T> K_mem(colA, colB);

                hoMatrix<T> A(rowA, colA, A_mem.begin());
                hoMatrix<T> B(rowB, colB, B_mem.begin());
                hoMatrix<T> K(colA, colB, K_mem.begin());

#pragma omp for
                for (ii = 0; ii<num; ii++)
                {
                    ho3DArray<T> src3D(RO, E1, E2, const_cast<T*>(src.begin()) + ii*RO*E1*E2);
                    ho3DArray<T> dst3D(RO, E1, E2, const_cast<T*>(dst.begin()) + ii*RO*E1*E2);

                    long long ro, e1, e2, row(0);
                    long long x, y, z;

                    for (e2 = startE2; e2 <= endE2; e2++)
                    {
                        for (e1 = startE1; e1 <= endE1; e1++)
                        {
                            for (ro = startRO; ro <= endRO; ro++)
                            {

                                size_t colInd(0);
                                for (z = -halfKz; z <= halfKz; z++)
                                {
                                    for (y = -halfKy; y <= halfKy; y++)
                                    {
                                        for (x = -halfKx; x <= halfKx; x++)
                                        {
                                            A(row, colInd++) = src3D(ro + x, e1 + y, e2 + z);
                                        }
                                    }
                                }

                                B(row, 0) = dst3D(ro, e1, e2);

                                row++;
                            }
                        }
                    }

                    Gadgetron::SolveLinearSystem_Tikhonov(A, B, K, thres);

                    memcpy(kernel.begin() + ii*kx*ky*kz, K.begin(), sizeof(T)*kx*ky*kz);
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in partial_fourier_FengHuang_calib_3d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_FengHuang_apply_kernel_3d(const hoNDArray<T>& kspaceConj, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, const hoNDArray<T>& kernel, hoNDArray<T>& kspace)
    {
        try
        {
            GADGET_CHECK_THROW(kspaceConj.dimensions_equal(&kspace));

            long long RO = (long long)kspace.get_size(0);
            long long E1 = (long long)kspace.get_size(1);
            long long E2 = (long long)kspace.get_size(2);

            long long kx = (long long)kernel.get_size(0);
            long long ky = (long long)kernel.get_size(1);
            long long kz = (long long)kernel.get_size(2);

            long long halfKx = kx / 2;
            long long halfKy = ky / 2;
            long long halfKz = kz / 2;

            long long num = kspace.get_number_of_elements()/(RO*E1*E2);

            long long rowD = RO*E1*E2 - ((endE2 - startE2 + 1) * (endE1 - startE1 + 1) * (endRO - startRO + 1));
            long long colD = kx*ky*kz;

            ho2DArray<long long> coeffX(colD, rowD);
            long long* pCx = coeffX.begin();

            ho2DArray<long long> coeffY(colD, rowD);
            long long* pCy = coeffY.begin();

            ho2DArray<long long> coeffZ(colD, rowD);
            long long* pCz = coeffZ.begin();

            long long ro, e1, e2;
            long long row(0);
            long long x, y, z;

            ho2DArray<long long> rowInd(3, rowD);
            long long* pRowInd = rowInd.begin();

            hoNDArray<long long> offsetX(colD);
            long long* pOffsetX = offsetX.begin();

            hoNDArray<long long> offsetY(colD);
            long long* pOffsetY = offsetY.begin();

            hoNDArray<long long> offsetZ(colD);
            long long* pOffsetZ = offsetZ.begin();

            long long colInd(0);
            for (z = -halfKz; z <= halfKz; z++)
            {
                for (y = -halfKy; y <= halfKy; y++)
                {
                    for (x = -halfKx; x <= halfKx; x++)
                    {
                        offsetX(colInd) = x;
                        offsetY(colInd) = y;
                        offsetZ(colInd) = z;
                        colInd++;
                    }
                }
            }

            long long* pRowIndCurr;
            for (e2 = 0; e2<E2; e2++)
            {
                for (e1 = 0; e1<E1; e1++)
                {
                    for (ro = 0; ro<RO; ro++)
                    {
                        if ((ro >= (long long)startRO) && (ro <= (long long)endRO) && (e1 >= (long long)startE1) && (e1 <= (long long)endE1) && (e2 >= (long long)startE2) && (e2 <= (long long)endE2))
                        {
                            continue;
                        }

                        pRowIndCurr = pRowInd + row * 3;

                        pRowIndCurr[0] = ro;
                        pRowIndCurr[1] = e1;
                        pRowIndCurr[2] = e2;

                        row++;
                    }
                }
            }

            long long r;
#pragma omp parallel for default(none) private(r) shared(rowD, colD, pCx, pCy, pCz, pRowInd, pRowIndCurr, pOffsetX, pOffsetY, pOffsetZ)
            for (r = 0; r<rowD; r++)
            {
                long long offsetC = r*colD;
                pRowIndCurr = pRowInd + r * 3;

                for (int colInd = 0; colInd<colD; colInd++)
                {
                    pCx[offsetC + colInd] = pRowIndCurr[0] + pOffsetX[colInd];
                    pCy[offsetC + colInd] = pRowIndCurr[1] + pOffsetY[colInd];
                    pCz[offsetC + colInd] = pRowIndCurr[2] + pOffsetZ[colInd];
                }
            }

#pragma omp parallel for default(none) private(r) shared(rowD, colD, pCx, pCy, pCz, RO, E1, E2)
            for (r = 0; r<rowD; r++)
            {
                for (int c = 0; c<colD; c++)
                {
                    long long offset = c + r*colD;

                    if (pCx[offset] < 0)
                    {
                        pCx[offset] += RO;
                    }
                    else if (pCx[offset] > RO - 1)
                    {
                        pCx[offset] -= RO;
                    }

                    if (pCy[offset] < 0)
                    {
                        pCy[offset] += E1;
                    }
                    else if (pCy[offset] > E1 - 1)
                    {
                        pCy[offset] -= E1;
                    }

                    if (pCz[offset] < 0)
                    {
                        pCz[offset] += E2;
                    }
                    else if (pCz[offset] > E2 - 1)
                    {
                        pCz[offset] -= E2;
                    }
                }
            }

            long long ii;
            int numOfThreads = (int)((num>4) ? 4 : num);
#pragma omp parallel default(none) private(ii) shared(num, RO, E1, E2, kspaceConj, kspace, kernel, rowD, colD, coeffX, coeffY, coeffZ, pCx, pCy, pCz) if ( num > 1 ) num_threads( numOfThreads )
            {
                hoNDArray<T> D_mem(rowD, colD);

                hoMatrix<T> D(rowD, colD, D_mem.begin());
                T* pD = D.begin();

                hoMatrix<T> K(colD, 1);
                hoMatrix<T> R(rowD, 1);

                Gadgetron::clear(D);
                Gadgetron::clear(K);
                Gadgetron::clear(R);

#pragma omp for
                for (ii = 0; ii<num; ii++)
                {
                    ho3DArray<T> src3D(RO, E1, E2, const_cast<T*>(kspaceConj.begin()) + ii*RO*E1*E2);
                    ho3DArray<T> dst3D(RO, E1, E2, kspace.begin() + ii*RO*E1*E2);

                    long long row;

                    for (row = 0; row<rowD; row++)
                    {
                        for (long long col = 0; col<colD; col++)
                        {
                            long long offset = col + row*colD;
                            pD[offset] = src3D(pCx[offset], pCy[offset], pCz[offset]);
                        }
                    }

                    memcpy(K.begin(), kernel.begin() + ii*colD, sizeof(T)*colD);

                    // R = D*K
                    Gadgetron::gemm(R, D, false, K, false);

                    size_t colCenter = colD / 2;

                    for (row = 0; row<rowD; row++)
                    {
                        dst3D(coeffX(colCenter, row), coeffY(colCenter, row), coeffZ(colCenter, row)) = R(row, 0);
                    }
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in partial_fourier_FengHuang_apply_kernel_3d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    void partial_fourier_FengHuang_3d(const hoNDArray<T>& kspace,
        size_t startRO, size_t endRO,
        size_t startE1, size_t endE1,
        size_t startE2, size_t endE2,
        size_t transit_band_RO,
        size_t transit_band_E1,
        size_t transit_band_E2,
        size_t kRO, size_t kE1, size_t kE2,
        double thres,
        hoNDArray<T>& res)
    {
        try
        {
            size_t RO = kspace.get_size(0);
            size_t E1 = kspace.get_size(1);
            size_t E2 = kspace.get_size(2);

            res = kspace;

            if (startRO >= RO || endRO >= RO || startRO >= endRO)
            {
                GWARN_STREAM("partial_fourier_FengHuang_3d(...) : (startRO >= RO || endRO >= RO || startRO >= endRO) ... ");
                startRO = 0;
                endRO = RO - 1;
            }

            if (startE1 >= E1 || endE1 >= E1 || startE1 >= endE1)
            {
                GWARN_STREAM("partial_fourier_FengHuang_3d(...) : (startE1 >= E1 || endE1 >= E1 || startE1 >= endE1) ... ");
                startE1 = 0;
                endE1 = E1 - 1;
            }

            if (startE2 >= E2 || endE2 >= E2 || startE2 >= endE2)
            {
                GWARN_STREAM("partial_fourier_FengHuang_3d(...) : (startE2 >= E2 || endE2 >= E2 || startE2 >= endE2) ... ");
                startE2 = 0;
                endE2 = E2 - 1;
            }

            if ((endRO - startRO + 1 == RO) && (endE1 - startE1 + 1 == E1) && (endE2 - startE2 + 1 == E2))
            {
                GWARN_STREAM("partial_fourier_FengHuang_3d(...) : (endRO - startRO + 1 == RO) && (endE1 - startE1 + 1 == E1) && (endE2 - startE2 + 1 == E2) ... ");
                return;
            }

            // compute the conjugate symmetric kspace
            hoNDArray<T> buffer3DT(kspace.get_dimensions());
            Gadgetron::partial_fourier_conjugate_symmetry_3d(kspace, buffer3DT);

            // find the symmetric region in the kspace
            size_t startSymRO, endSymRO;
            Gadgetron::find_symmetric_sampled_region(startRO, endRO, RO / 2, startSymRO, endSymRO);

            size_t startSymE1, endSymE1;
            Gadgetron::find_symmetric_sampled_region(startE1, endE1, E1 / 2, startSymE1, endSymE1);

            size_t startSymE2, endSymE2;
            Gadgetron::find_symmetric_sampled_region(startE2, endE2, E2 / 2, startSymE2, endSymE2);

            // the reference kspace for kernel estimation
            hoNDArray<T> src, dst;
            Gadgetron::vector_td<size_t, 3> start, size;

            start[0] = startSymRO;
            start[1] = startSymE1;
            start[2] = startSymE2;

            size[0] = endSymRO - startSymRO + 1;
            size[1] = endSymE1 - startSymE1 + 1;
            size[2] = endSymE2 - startSymE2 + 1;;

            Gadgetron::crop(start, size, &buffer3DT, &src);
            Gadgetron::crop(start, size, const_cast< hoNDArray<T>* >(&kspace), &dst);

            // estimate the kernels
            hoNDArray<T> kernel;
            Gadgetron::partial_fourier_FengHuang_calib_3d(src, dst, kRO, kE1, kE2, thres, kernel);

            // perform the recon
            if (transit_band_RO == 0 && transit_band_E1 == 0 && transit_band_E2 == 0)
            {
                Gadgetron::partial_fourier_FengHuang_apply_kernel_3d(buffer3DT, startRO, endRO, startE1, endE1, startE2, endE2, kernel, res);
            }
            else
            {
                long long sRO(startRO), eRO(endRO), sE1(startE1), eE1(endE1), sE2(startE2), eE2(endE2);

                if (startRO > 0)
                {
                    startRO += transit_band_RO;
                    if (startRO > RO) startRO = 0;
                }

                if (endRO < RO - 1)
                {
                    endRO -= transit_band_RO;
                }

                if (startRO > endRO)
                {
                    startRO = 0;
                    endRO = RO - 1;
                }

                if (startE1 > 0)
                {
                    startE1 += transit_band_E1;
                    if (startE1 > E1) startE1 = 0;
                }

                if (endE1 < E1 - 1)
                {
                    endE1 -= transit_band_E1;
                }

                if (startE1 > endE1)
                {
                    startE1 = 0;
                    endE1 = E1 - 1;
                }

                if (startE2 > 0)
                {
                    startE2 += transit_band_E2;
                    if (startE2 > E2) startE2 = 0;
                }

                if (endE2 < E2 - 1)
                {
                    endE2 -= transit_band_E2;
                }

                if (startE2 > endE2)
                {
                    startE2 = 0;
                    endE2 = E2 - 1;
                }

                Gadgetron::partial_fourier_FengHuang_apply_kernel_3d(buffer3DT, startRO, endRO, startE1, endE1, startE2, endE2, kernel, res);


                Gadgetron::partial_fourier_transition_band_3d(kspace, res, sRO, eRO, sE1, eE1, sE2, eE2, transit_band_RO, transit_band_E1, transit_band_E2);
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in partial_fourier_FengHuang_3d(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    template <typename T>
    partialFourierCartesianHandler<T>::partialFourierCartesianHandler()
    {
        start_RO_ = 0;
        end_RO_ = 0;

        start_E1_ = 0;
        end_E1_ = 0;

        start_E2_ = 0;
        end_E2_ = 0;

        transit_band_RO_ = 24;
        transit_band_E1_ = 24;
        transit_band_E2_ = 24;

        verbose_ = false;
    }

    template <typename T>
    partialFourierCartesianHandler<T>::~partialFourierCartesianHandler()
    {
    }

    template <typename T>
    void partialFourierCartesianHandler<T>::dump(std::ostream& os) const
    {
        using namespace std;
        os << "-------------------------------------------------------------------------------" << endl;
        os << "partialFourierCartesianHandler will fill the unsampled region in case of partial fourier or asymmetric echo for the signle/multi-channel complex images ... " << endl;
        os << "partial_fourier_2d and partial_fourier_3d functions should be implemented ... " << endl;
        os << "-------------------------------------------------------------------------------" << endl;
    }

    // ------------------------------------------------------------------------

    template <typename T> 
    partialFourierCartesianFilterHandler<T>::partialFourierCartesianFilterHandler() : BaseClass()
    {
        filter_pf_density_comp_ = false;
        filter_pf_width_RO_ = 0.15;
        filter_pf_width_E1_ = 0.15;
        filter_pf_width_E2_ = 0.15;
    }

    template <typename T>
    partialFourierCartesianFilterHandler<T>::~partialFourierCartesianFilterHandler()
    {
    }

    template <typename T>
    void partialFourierCartesianFilterHandler<T>::partial_fourier(const hoNDArray<T>& kspace, hoNDArray<T>& res)
    {
        try
        {
            size_t RO = kspace.get_size(0);
            size_t E1 = kspace.get_size(1);
            size_t E2 = kspace.get_size(2);

            size_t lenRO = end_RO_ - start_RO_ + 1;
            if (filter_RO_pf_.get_size(0) != RO && lenRO<RO)
            {
                Gadgetron::generate_asymmetric_filter(RO, start_RO_, end_RO_, filter_RO_pf_, ISMRMRD_FILTER_TAPERED_HANNING, (size_t)(RO*filter_pf_width_RO_), filter_pf_density_comp_);
            }

            size_t lenE1 = end_E1_ - start_E1_ + 1;
            if (filter_E1_pf_.get_size(0) != E1 && lenE1<E1)
            {
                Gadgetron::generate_asymmetric_filter(E1, start_E1_, end_E1_, filter_E1_pf_, ISMRMRD_FILTER_TAPERED_HANNING, (size_t)(E1*filter_pf_width_E1_), filter_pf_density_comp_);
            }

            res = kspace;

            size_t lenE2 = end_E2_ - start_E2_ + 1;

            if (E2 > 1)
            {
                if (filter_E2_pf_.get_size(0) != E2 && lenE2 < E2)
                {
                    Gadgetron::generate_asymmetric_filter(E2, start_E2_, end_E2_, filter_E2_pf_, ISMRMRD_FILTER_TAPERED_HANNING, (size_t)(E2*filter_pf_width_E2_), filter_pf_density_comp_);
                }

                if ((filter_RO_pf_.get_number_of_elements() == RO) && (filter_E1_pf_.get_number_of_elements() == E1) && (filter_E2_pf_.get_number_of_elements() == E2))
                {
                    Gadgetron::apply_kspace_filter_ROE1E2(kspace, filter_RO_pf_, filter_E1_pf_, filter_E2_pf_, res);
                }
                else
                {
                    if ((filter_RO_pf_.get_number_of_elements() == RO)
                        || (filter_E1_pf_.get_number_of_elements() == E1)
                        || (filter_E2_pf_.get_number_of_elements() == E2))
                    {
                        hoNDArray<T> kspace_copy(kspace);

                        hoNDArray<T>* pSrc = const_cast<hoNDArray<T>*>(&kspace_copy);
                        hoNDArray<T>* pDst = &res;

                        bool filterPerformed = false;

                        if (filter_RO_pf_.get_number_of_elements() == RO)
                        {
                            Gadgetron::apply_kspace_filter_RO(*pSrc, filter_RO_pf_, *pDst);
                            std::swap(pSrc, pDst);
                            filterPerformed = true;
                        }

                        if (filter_E1_pf_.get_number_of_elements() == E1)
                        {
                            Gadgetron::apply_kspace_filter_E1(*pSrc, filter_E1_pf_, *pDst);
                            std::swap(pSrc, pDst);
                            filterPerformed = true;
                        }

                        if (filter_E2_pf_.get_number_of_elements() == E2)
                        {
                            Gadgetron::apply_kspace_filter_E2(*pSrc, filter_E2_pf_, *pDst);
                            std::swap(pSrc, pDst);
                            filterPerformed = true;
                        }

                        if (filterPerformed && pDst != &res)
                        {
                            res = *pDst;
                        }
                    }
                }
            }
            else
            {
                if ((filter_RO_pf_.get_number_of_elements() == RO) && (filter_E1_pf_.get_number_of_elements() == E1))
                {
                    Gadgetron::apply_kspace_filter_ROE1(kspace, filter_RO_pf_, filter_E1_pf_, res);
                }
                else
                {
                    if (filter_RO_pf_.get_number_of_elements() == RO && filter_E1_pf_.get_number_of_elements() != E1)
                    {
                        Gadgetron::apply_kspace_filter_RO(kspace, filter_RO_pf_, res);
                    }
                    else if (filter_RO_pf_.get_number_of_elements() != RO && filter_E1_pf_.get_number_of_elements() == E1)
                    {
                        Gadgetron::apply_kspace_filter_E1(kspace, filter_E1_pf_, res);
                    }
                }
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in partialFourierCartesianFilterHandler<T>::partial_fourier(...) ... ");
        }
    }

    template <typename T>
    void partialFourierCartesianFilterHandler<T>::dump(std::ostream& os) const
    {
        using namespace std;
        os << "-------------------------------------------------------------------------------" << endl;
        os << "partialFourierCartesianFilterHandler ... " << endl;
        os << "filter_pf_density_comp_ is " << filter_pf_density_comp_ << endl;
        os << "filter_pf_width_RO_ is " << filter_pf_width_RO_ << endl;
        os << "filter_pf_width_E1_ is " << filter_pf_width_E1_ << endl;
        os << "filter_pf_width_E2_ is " << filter_pf_width_E2_ << endl;
        os << "-------------------------------------------------------------------------------" << endl;
    }

    template class EXPORTMRICORE partialFourierCartesianFilterHandler < std::complex<float> >;
    template class EXPORTMRICORE partialFourierCartesianFilterHandler < std::complex<double> >;

    // ------------------------------------------------------------------------

    template <typename T>
    partialFourierCartesianPOCSHandler<T>::partialFourierCartesianPOCSHandler() : BaseClass()
    {
        iter_ = 5;
        thres_ = 0.001;
    }

    template <typename T>
    partialFourierCartesianPOCSHandler<T>::~partialFourierCartesianPOCSHandler()
    {
    }

    template <typename T>
    void partialFourierCartesianPOCSHandler<T>::partial_fourier(const hoNDArray<T>& kspace, hoNDArray<T>& res)
    {
        try
        {
            size_t E2 = kspace.get_size(2);

            if (E2 > 1)
            {
                Gadgetron::partial_fourier_POCS(kspace, start_RO_, end_RO_, start_E1_, end_E1_, start_E2_, end_E2_, transit_band_RO_, transit_band_E1_, transit_band_E2_, iter_, thres_, true, res, verbose_);
            }
            else
            {
                Gadgetron::partial_fourier_POCS(kspace, start_RO_, end_RO_, start_E1_, end_E1_, 0, 0, transit_band_RO_, transit_band_E1_, 0, iter_, thres_, false, res, verbose_);
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in partialFourierCartesianPOCSHandler<T>::partial_fourier(...) ... ");
        }
    }

    template <typename T>
    void partialFourierCartesianPOCSHandler<T>::dump(std::ostream& os) const
    {
        using namespace std;
        os << "-------------------------------------------------------------------------------" << endl;
        os << "partialFourierCartesianPOCSHandler ... " << endl;
        os << "iter_ is " << iter_ << endl;
        os << "thres_ is " << thres_ << endl;
        os << "-------------------------------------------------------------------------------" << endl;
    }

    template class EXPORTMRICORE partialFourierCartesianPOCSHandler < std::complex<float> >;
    template class EXPORTMRICORE partialFourierCartesianPOCSHandler < std::complex<double> >;

    // ------------------------------------------------------------------------

    template <typename T>
    partialFourierCartesianFengHuangHandler<T>::partialFourierCartesianFengHuangHandler() : BaseClass()
    {
        kRO_ = 5;
        kE1_ = 5;
        kE2_ = 5;
        thres_ = 0.01;
    }

    template <typename T>
    partialFourierCartesianFengHuangHandler<T>::~partialFourierCartesianFengHuangHandler()
    {
    }

    template <typename T>
    void partialFourierCartesianFengHuangHandler<T>::partial_fourier(const hoNDArray<T>& kspace, hoNDArray<T>& res)
    {
        try
        {
            size_t E2 = kspace.get_size(2);

            if (E2 > 1)
            {
                Gadgetron::partial_fourier_FengHuang_3d(kspace, start_RO_, end_RO_, start_E1_, end_E1_, start_E2_, end_E2_, transit_band_RO_, transit_band_E1_, transit_band_E2_, kRO_, kE1_, kE2_, thres_, res);
            }
            else
            {
                Gadgetron::partial_fourier_FengHuang_2d(kspace, start_RO_, end_RO_, start_E1_, end_E1_, transit_band_RO_, transit_band_E1_, kRO_, kE1_, thres_, res);
            }
        }
        catch (...)
        {
            GADGET_THROW("Errors happened in partialFourierCartesianFengHuangHandler<T>::partial_fourier(...) ... ");
        }
    }

    template <typename T>
    void partialFourierCartesianFengHuangHandler<T>::dump(std::ostream& os) const
    {
        using namespace std;
        os << "-------------------------------------------------------------------------------" << endl;
        os << "partialFourierCartesianFengHuangHandler ... " << endl;
        os << "kRO_ is " << kRO_ << endl;
        os << "kE1_ is " << kE1_ << endl;
        os << "kE2_ is " << kE2_ << endl;
        os << "thres_ is " << thres_ << endl;
        os << "-------------------------------------------------------------------------------" << endl;
    }

    template class EXPORTMRICORE partialFourierCartesianFengHuangHandler < std::complex<float> >;
    template class EXPORTMRICORE partialFourierCartesianFengHuangHandler < std::complex<double> >;

    // ------------------------------------------------------------------------
}
