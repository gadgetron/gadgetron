
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

    /// reset the kspace for acquired region
    template <typename T>
    void partial_fourier_reset_kspace(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2)
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
            GADGET_THROW("Errors happened in partial_fourier_reset_kspace(...) ... ");
        }
    }

    // ------------------------------------------------------------------------

    /// create and apply the transition band
    template <typename T>
    void partial_fourier_transition_band(const hoNDArray<T>& src, hoNDArray<T>& dst, size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2, size_t transBandRO, size_t transBandE1, size_t transBandE2)
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
            hoNDArray<T> filter_dst_RO(RO), filter_dst_E1(E1), filter_dst_E2(E2);
            size_t ii;

            if (startRO == 0 && endRO == RO - 1)
            {
                Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_NONE, transBandRO, densityComp);
            }
            else
            {
                Gadgetron::generate_asymmetric_filter(RO, startRO, endRO, filter_src_RO, ISMRMRD_FILTER_TAPERED_HANNING, transBandRO, densityComp);
            }

            T midValue = filter_src_RO(RO / 2);
            T scalFactor = T(1.0) / midValue;
            Gadgetron::scal(scalFactor, filter_src_RO);

            for (ii = 0; ii<RO; ii++)
            {
                filter_dst_RO(ii) = T(1.0) - filter_src_RO(ii);
            }

            if (startE1 == 0 && endE1 == E1 - 1)
            {
                Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_NONE, transBandE1, densityComp);
            }
            else
            {
                Gadgetron::generate_asymmetric_filter(E1, startE1, endE1, filter_src_E1, ISMRMRD_FILTER_TAPERED_HANNING, transBandE1, densityComp);
            }

            midValue = filter_src_E1(E1 / 2);
            scalFactor = T(1.0) / midValue;
            Gadgetron::scal(scalFactor, filter_src_E1);

            for (ii = 0; ii<E1; ii++)
            {
                filter_dst_E1(ii) = T(1.0) - filter_src_E1(ii);
            }

            if (E2 > 1)
            {
                if (startE2 == 0 && endE2 == E2 - 1)
                {
                    Gadgetron::generate_asymmetric_filter(E2, startE2, endE2, filter_src_E2, ISMRMRD_FILTER_NONE, transBandE2, densityComp);
                }
                else
                {
                    Gadgetron::generate_asymmetric_filter(E2, startE2, endE2, filter_src_E2, ISMRMRD_FILTER_TAPERED_HANNING, transBandE2, densityComp);
                }

                midValue = filter_src_E2(E2 / 2);
                scalFactor = T(1.0) / midValue;
                Gadgetron::scal(scalFactor, filter_src_E2);

                for (ii = 0; ii<E2; ii++)
                {
                    filter_dst_E2(ii) = T(1.0) - filter_src_E2(ii);
                }
            }
            else
            {
                filter_src_E2.create(1);
                filter_src_E2(0) = 1;

                filter_dst_E2.create(1);
                filter_dst_E2(0) = 1;
            }

            hoNDArray<T> srcFiltered(src), dstFiltered(dst);
            if (endRO <= RO - 1 && startE1 == 0 && endE1 == E1 - 1 && startE2 == 0 && endE1 == E2 - 1)
            {
                Gadgetron::apply_kspace_filter_RO(src, filter_src_RO, srcFiltered);
                Gadgetron::apply_kspace_filter_RO(dst, filter_dst_RO, dstFiltered);
            }
            else if (startRO == 0 && endRO == RO - 1 && endE1 <= E1 - 1 && startE2 == 0 && endE1 == E2 - 1)
            {
                Gadgetron::apply_kspace_filter_E1(src, filter_src_E1, srcFiltered);
                Gadgetron::apply_kspace_filter_E1(dst, filter_dst_E1, dstFiltered);
            }
            else if (startRO == 0 && endRO == RO - 1 && startE1 == 0 && endE1 == E1 - 1 && endE1 <= E2 - 1)
            {
                Gadgetron::apply_kspace_filter_E2(src, filter_src_E2, srcFiltered);
                Gadgetron::apply_kspace_filter_E2(dst, filter_dst_E2, dstFiltered);
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
            GADGET_THROW("Errors happened in partial_fourier_transition_band(...) ... ");
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
                            hoNDArray<T>& res)
    {
        try
        {
            size_t RO = kspace.get_size(0);
            size_t E1 = kspace.get_size(1);
            size_t E2 = kspace.get_size(2);

            bool is3D = (E2 > 1);

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
            else
            {
                startE2 = 0;
                endE2 = 0;
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
            hoNDArray<typename realType<T>::Type> mag(kspace.dimensions());
            hoNDArray<T> magComplex(kspace.dimensions());

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
                partial_fourier_reset_kspace(kspace, kspaceIter, startRO, endRO, startE1, endE1, startE2, endE2);

                if (is3D)
                {
                    // update complex image
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kspaceIter, complexImPOCS);
                }
                else
                {
                    // update complex image
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspaceIter, complexImPOCS);
                }

                // compute threshold to stop the iteration
                Gadgetron::subtract(complexImPOCS, complexIm, buffer_partial_fourier);
                auto prev = Gadgetron::nrm2(complexIm);
                auto diff = Gadgetron::nrm2(buffer_partial_fourier);

                auto t = diff / prev;

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
                Gadgetron::partial_fourier_transition_band(kspace, buffer_partial_fourier_Iter, startRO, endRO, startE1, endE1, startE2, endE2, transit_band_RO, transit_band_E1, transit_band_E2);
                res = buffer_partial_fourier_Iter;
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in partial_fourier_POCS(...) ... ");
        }
    }

    template EXPORTMRICORE void partial_fourier_POCS(const hoNDArray< std::complex<float> >& kspace,
        size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
        size_t transit_band_RO, size_t transit_band_E1, size_t transit_band_E2,
        size_t iter, double thres, hoNDArray< std::complex<float> >& res);

    template EXPORTMRICORE void partial_fourier_POCS(const hoNDArray< std::complex<double> >& kspace,
        size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
        size_t transit_band_RO, size_t transit_band_E1, size_t transit_band_E2,
        size_t iter, double thres, hoNDArray< std::complex<double> >& res);

    // ------------------------------------------------------------------------

    template <typename T>
    void compute_partial_fourier_filter(size_t len, size_t start, size_t end, double filter_pf_width, bool filter_pf_density_comp, hoNDArray<T>& filter_pf)
    {
        try
        {
            if (start == 0 || end == len - 1)
            {
                Gadgetron::generate_asymmetric_filter(len, start, end, filter_pf, ISMRMRD_FILTER_TAPERED_HANNING, (size_t)(len*filter_pf_width), filter_pf_density_comp);
            }
            else
            {
                size_t fil_len = end - start + 1;
                size_t len_end = 2 * (end - len / 2 + 1);
                size_t len_start = 2 * (len / 2 - start + 1);

                if (len_end > len_start)
                {
                    hoNDArray<T> fil(len_end);
                    Gadgetron::generate_asymmetric_filter(len_end, len_end - fil_len, len_end - 1, fil, ISMRMRD_FILTER_TAPERED_HANNING, (size_t)(len_end*filter_pf_width), filter_pf_density_comp);
                    Gadgetron::pad(len, fil, filter_pf);
                }
                else
                {
                    hoNDArray<T> fil(len_start);
                    Gadgetron::generate_asymmetric_filter(len_start, 0, fil_len - 1, fil, ISMRMRD_FILTER_TAPERED_HANNING, (size_t)(len_start*filter_pf_width), filter_pf_density_comp);
                    Gadgetron::pad(len, fil, filter_pf);
                }
            }
        }
        catch (...)
        {

        }
    }

    template <typename T>
    void partial_fourier_filter(const hoNDArray<T>& kspace,
                                size_t startRO, size_t endRO,
                                size_t startE1, size_t endE1,
                                size_t startE2, size_t endE2,
                                double filter_pf_width_RO,
                                double filter_pf_width_E1,
                                double filter_pf_width_E2,
                                bool filter_pf_density_comp,
                                hoNDArray<T>& filter_pf_RO,
                                hoNDArray<T>& filter_pf_E1,
                                hoNDArray<T>& filter_pf_E2,
                                hoNDArray<T>& res)
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);

        size_t lenRO = endRO - startRO + 1;
        if (filter_pf_RO.get_size(0) != RO && lenRO<RO)
        {
            Gadgetron::compute_partial_fourier_filter(RO, startRO, endRO, filter_pf_width_RO, filter_pf_density_comp, filter_pf_RO);
        }

        size_t lenE1 = endE1 - startE1 + 1;
        if (filter_pf_E1.get_size(0) != E1 && lenE1<E1)
        {
            Gadgetron::compute_partial_fourier_filter(E1, startE1, endE1, filter_pf_width_E1, filter_pf_density_comp, filter_pf_E1);
        }

        res = kspace;

        size_t lenE2 = endE2 - startE2 + 1;

        if (E2 > 1)
        {
            if (filter_pf_E2.get_size(0) != E2 && lenE2 < E2)
            {
                Gadgetron::compute_partial_fourier_filter(E2, startE2, endE2, filter_pf_width_E2, filter_pf_density_comp, filter_pf_E2);
            }

            if ((filter_pf_RO.get_number_of_elements() == RO) && (filter_pf_E1.get_number_of_elements() == E1) && (filter_pf_E2.get_number_of_elements() == E2))
            {
                Gadgetron::apply_kspace_filter_ROE1E2(kspace, filter_pf_RO, filter_pf_E1, filter_pf_E2, res);
            }
            else
            {
                if ((filter_pf_RO.get_number_of_elements() == RO)
                    || (filter_pf_E1.get_number_of_elements() == E1)
                    || (filter_pf_E2.get_number_of_elements() == E2))
                {
                    hoNDArray<T> kspace_copy(kspace);

                    hoNDArray<T>* pSrc = &kspace_copy;
                    hoNDArray<T>* pDst = &res;

                    bool filterPerformed = false;

                    if (filter_pf_RO.get_number_of_elements() == RO)
                    {
                        Gadgetron::apply_kspace_filter_RO(*pSrc, filter_pf_RO, *pDst);
                        std::swap(pSrc, pDst);
                        filterPerformed = true;
                    }

                    if (filter_pf_E1.get_number_of_elements() == E1)
                    {
                        Gadgetron::apply_kspace_filter_E1(*pSrc, filter_pf_E1, *pDst);
                        std::swap(pSrc, pDst);
                        filterPerformed = true;
                    }

                    if (filter_pf_E2.get_number_of_elements() == E2)
                    {
                        Gadgetron::apply_kspace_filter_E2(*pSrc, filter_pf_E2, *pDst);
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
            if ((filter_pf_RO.get_number_of_elements() == RO) && (filter_pf_E1.get_number_of_elements() == E1))
            {
                Gadgetron::apply_kspace_filter_ROE1(kspace, filter_pf_RO, filter_pf_E1, res);
            }
            else
            {
                if (filter_pf_RO.get_number_of_elements() == RO && filter_pf_E1.get_number_of_elements() != E1)
                {
                    Gadgetron::apply_kspace_filter_RO(kspace, filter_pf_RO, res);
                }
                else if (filter_pf_RO.get_number_of_elements() != RO && filter_pf_E1.get_number_of_elements() == E1)
                {
                    Gadgetron::apply_kspace_filter_E1(kspace, filter_pf_E1, res);
                }
            }
        }

    }

    template EXPORTMRICORE void partial_fourier_filter(const hoNDArray< std::complex<float> >& kspace,
        size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
        double filter_pf_width_RO, double filter_pf_width_E1, double filter_pf_width_E2, bool filter_pf_density_comp,
        hoNDArray< std::complex<float> >& filter_pf_RO, hoNDArray< std::complex<float> >& filter_pf_E1, hoNDArray< std::complex<float> >& filter_pf_E2,
        hoNDArray< std::complex<float> >& res);

    template EXPORTMRICORE void partial_fourier_filter(const hoNDArray< std::complex<double> >& kspace,
        size_t startRO, size_t endRO, size_t startE1, size_t endE1, size_t startE2, size_t endE2,
        double filter_pf_width_RO, double filter_pf_width_E1, double filter_pf_width_E2, bool filter_pf_density_comp,
        hoNDArray< std::complex<double> >& filter_pf_RO, hoNDArray< std::complex<double> >& filter_pf_E1, hoNDArray< std::complex<double> >& filter_pf_E2,
        hoNDArray< std::complex<double> >& res);

    // ------------------------------------------------------------------------
}
