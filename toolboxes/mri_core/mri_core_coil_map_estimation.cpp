
/** \file   mri_core_coil_map_estimation.cpp
    \brief  Implementation MRI coil sensitivity map estimation.
    \author Hui Xue
*/

#include "mri_core_coil_map_estimation.h"

namespace Gadgetron
{

template<typename T> 
void coil_map_2d_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t power)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long CHA = data.get_size(2);

        long long N = data.get_number_of_elements() / (RO*E1*CHA);
        GADGET_CHECK_THROW(N == 1);

        const T* pData = data.begin();

        if (!data.dimensions_equal(&coilMap))
        {
            coilMap = data;
        }
        T* pSen = coilMap.begin();

        if (ks % 2 != 1)
        {
            ks++;
        }

        size_t kss = ks*ks;
        long long halfKs = (long long)ks / 2;

        int e1;

        #pragma omp parallel private(e1) shared(ks, RO, E1, CHA, pSen, pData, halfKs, power, kss)
        {
            hoNDArray<T> D(ks*ks, CHA);
            T* pD = D.begin();

            hoNDArray<T> DC(ks*ks, CHA);
            T* pDC = DC.begin();

            hoNDArray<T> DH_D(CHA, CHA);
            Gadgetron::clear(DH_D);

            hoNDArray<T> U1(ks*ks, 1);
            T* pU1 = U1.begin();

            hoNDArray<T> V1(CHA, 1);
            T* pV1 = V1.begin();

            hoNDArray<T> V(CHA, 1);

            Gadgetron::clear(D);
            Gadgetron::clear(DC);
            Gadgetron::clear(DH_D);
            Gadgetron::clear(U1);
            Gadgetron::clear(V1);
            Gadgetron::clear(V);

            T phaseU1;

            value_type v1Norm(1), u1Norm(1);

            long long cha, ro, kro, ke1, de1, dro;
            size_t po;

            #pragma omp for
            for (e1 = 0; e1<(int)E1; e1++)
            {
                for (ro = 0; ro<(long long)RO; ro++)
                {
                    // fill the data matrix D
                    if (e1 >= halfKs && e1<E1 - halfKs && ro >= halfKs && ro<RO - halfKs)
                    {
                        for (cha = 0; cha<CHA; cha++)
                        {
                            const T* pDataCurr = pData + cha*RO*E1;
                            int ind = 0;
                            for (ke1 = -halfKs; ke1 <= halfKs; ke1++)
                            {
                                de1 = e1 + ke1;
                                for (kro = -halfKs; kro <= halfKs; kro++)
                                {
                                    pD[ind + cha*kss] = pDataCurr[de1*RO + ro + kro];
                                    ind++;
                                }
                            }
                        }
                    }
                    else
                    {
                        for (cha = 0; cha<CHA; cha++)
                        {
                            const T* pDataCurr = pData + cha*RO*E1;
                            int ind = 0;
                            for (ke1 = -halfKs; ke1 <= halfKs; ke1++)
                            {
                                de1 = e1 + ke1;
                                if (de1 < 0) de1 += E1;
                                if (de1 >= E1) de1 -= E1;

                                for (kro = -halfKs; kro <= halfKs; kro++)
                                {
                                    dro = ro + kro;
                                    if (dro < 0) dro += RO;
                                    if (dro >= RO) dro -= RO;

                                    pD[ind + cha*kss] = pDataCurr[de1*RO + dro];
                                    ind++;
                                }
                            }
                        }
                    }

                    T* pTmp;
                    for (cha = 0; cha<CHA; cha++)
                    {
                        pTmp = pD + cha*kss;
                        pV1[cha] = pTmp[0];
                        for (po = 1; po<kss; po++)
                        {
                            pV1[cha] += pTmp[po];
                        }
                    }

                    value_type sum(0);
                    for (cha = 0; cha<CHA; cha++)
                    {
                        const T& c = pV1[cha];
                        const value_type re = c.real();
                        const value_type im = c.imag();
                        sum += ((re*re) + (im * im));
                    }
                    v1Norm = std::sqrt(sum);

                    value_type v1NormInv = (value_type)1.0 / v1Norm;
                    for (cha = 0; cha<CHA; cha++)
                    {
                        pV1[cha] *= v1NormInv;
                    }

                    memcpy(pDC, pD, sizeof(T)*ks*ks*CHA);
                    gemm(DH_D, DC, true, D, false);

                    for (po = 0; po<power; po++)
                    {
                        gemm(V, DH_D, false, V1, false);
                        memcpy(V1.begin(), V.begin(), V.get_number_of_bytes());

                        sum = 0;
                        for (cha = 0; cha<CHA; cha++)
                        {
                            const T& c = pV1[cha];
                            const value_type re = c.real();
                            const value_type im = c.imag();
                            sum += ((re*re) + (im * im));
                        }
                        v1Norm = std::sqrt(sum);

                        value_type v1NormInv = (value_type)1.0 / v1Norm;
                        for (cha = 0; cha<CHA; cha++)
                        {
                            pV1[cha] *= v1NormInv;
                        }
                    }

                    gemm(U1, D, false, V1, false);

                    phaseU1 = pU1[0];
                    for (po = 1; po<kss; po++)
                    {
                        phaseU1 += pU1[po];
                    }
                    phaseU1 /= std::abs(phaseU1);

                    const value_type c = phaseU1.real();
                    const value_type d = phaseU1.imag();

                    for (cha = 0; cha<CHA; cha++)
                    {
                        const T& v = pV1[cha];
                        const value_type a = v.real();
                        const value_type b = v.imag();

                        reinterpret_cast< value_type(&)[2] >(pV1[cha])[0] = a*c + b*d;
                        reinterpret_cast< value_type(&)[2] >(pV1[cha])[1] = a*d - b*c;
                    }

                    for (cha = 0; cha<CHA; cha++)
                    {
                        pSen[cha*RO*E1 + e1*RO + ro] = V1(cha, 0);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in coil_map_2d_Inati(...) ... ");
        throw;
    }
}

template EXPORTMRICORE void coil_map_2d_Inati(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t power);
template EXPORTMRICORE void coil_map_2d_Inati(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t power);

// ------------------------------------------------------------------------

template<typename T> 
void coil_map_3d_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t power)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long E2 = data.get_size(2);
        long long CHA = data.get_size(3);

        long long N = data.get_number_of_elements() / (RO*E1*E2*CHA);
        GADGET_CHECK_THROW(N == 1);

        const T* pData = data.begin();

        if (!data.dimensions_equal(&coilMap))
        {
            coilMap = data;
        }
        T* pSen = coilMap.begin();

        if (ks % 2 != 1)
        {
            ks++;
        }

        size_t kss = ks*ks*ks;
        long long halfKs = (long long)ks / 2;

        long long e2;

        #pragma omp parallel default(none) private(e2) shared(ks, RO, E1, E2, CHA, pSen, pData, halfKs, power, kss)
        {
            hoMatrix<T> D(kss, CHA);
            hoMatrix<T> DC(kss, CHA);
            hoMatrix<T> DH_D(CHA, CHA);

            hoMatrix<T> U1(kss, 1);
            hoMatrix<T> V1(CHA, 1);
            hoMatrix<T> V(CHA, 1);

            Gadgetron::clear(D);
            Gadgetron::clear(DC);
            Gadgetron::clear(DH_D);
            Gadgetron::clear(U1);
            Gadgetron::clear(V1);
            Gadgetron::clear(V);

            T phaseU1;

            value_type v1Norm(1);

            long long cha, ro, e1, kro, dro, ke1, de1, ke2, de2;
            size_t po;

            #pragma omp for
            for (e2 = 0; e2<(long long)E2; e2++)
            {
                for (e1 = 0; e1<(long long)E1; e1++)
                {
                    for (ro = 0; ro<(long long)RO; ro++)
                    {
                        // fill the data matrix D
                        if (e2 >= halfKs && e2<E2 - halfKs && e1 >= halfKs && e1<E1 - halfKs && ro >= halfKs && ro<RO - halfKs)
                        {
                            for (cha = 0; cha<CHA; cha++)
                            {
                                const T* pDataCurr = pData + cha*RO*E1*E2;
                                long long ind = 0;
                                for (ke2 = -halfKs; ke2 <= halfKs; ke2++)
                                {
                                    de2 = e2 + ke2;
                                    for (ke1 = -halfKs; ke1 <= halfKs; ke1++)
                                    {
                                        de1 = e1 + ke1;
                                        for (kro = -halfKs; kro <= halfKs; kro++)
                                        {
                                            D(ind++, cha) = pDataCurr[de2*RO*E1 + de1*RO + ro + kro];
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            for (cha = 0; cha<CHA; cha++)
                            {
                                const T* pDataCurr = pData + cha*RO*E1*E2;
                                long long ind = 0;
                                for (ke2 = -halfKs; ke2 <= halfKs; ke2++)
                                {
                                    de2 = e2 + ke2;
                                    if (de2 < 0) de2 += E2;
                                    if (de2 >= E2) de2 -= E2;

                                    for (ke1 = -halfKs; ke1 <= halfKs; ke1++)
                                    {
                                        de1 = e1 + ke1;
                                        if (de1 < 0) de1 += E1;
                                        if (de1 >= E1) de1 -= E1;

                                        for (kro = -halfKs; kro <= halfKs; kro++)
                                        {
                                            dro = ro + kro;
                                            if (dro < 0) dro += RO;
                                            if (dro >= RO) dro -= RO;

                                            D(ind++, cha) = pDataCurr[de2*RO*E1 + de1*RO + dro];
                                        }
                                    }
                                }
                            }
                        }

                        // compute V1
                        D.sumOverCol(V1);
                        norm2(V1, v1Norm);
                        scal((value_type)1.0 / v1Norm, V1);

                        memcpy(DC.begin(), D.begin(), sizeof(T)*kss*CHA);
                        gemm(DH_D, DC, true, D, false);
                        // gemm(DH_D, D, true, D, false);

                        for (po = 0; po<power; po++)
                        {
                            gemm(V, DH_D, false, V1, false);
                            V1 = V;
                            norm2(V1, v1Norm);
                            scal((value_type)1.0 / v1Norm, V1);
                        }

                        // compute U1
                        gemm(U1, D, false, V1, false);

                        phaseU1 = U1(0, 0);
                        for (po = 1; po<kss; po++)
                        {
                            phaseU1 += U1(po, 0);
                        }
                        phaseU1 /= std::abs(phaseU1);

                        // put the mean object phase to coil map
                        conjugate(V1, V1);
                        scal(phaseU1, V1);

                        for (cha = 0; cha<CHA; cha++)
                        {
                            pSen[cha*RO*E1*E2 + e2*RO*E1 + e1*RO + ro] = V1(cha, 0);
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in coil_map_3d_Inati(...) ... ");
        throw;
    }
}

template EXPORTMRICORE void coil_map_3d_Inati(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t power);
template EXPORTMRICORE void coil_map_3d_Inati(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t power);

// ------------------------------------------------------------------------

template<typename T> 
void coil_map_2d_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t iterNum, typename realType<T>::Type thres)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long CHA = data.get_size(2);

        long long N = data.get_number_of_elements() / (RO*E1*CHA);
        GADGET_CHECK_THROW(N == 1);

        const T* pData = data.begin();

        if (!data.dimensions_equal(&coilMap))
        {
            coilMap = data;
        }

        // create convolution kernel
        hoNDArray<T> ker(ks, ks);
        Gadgetron::fill(ker, T((value_type)1.0 / (ks*ks)));

        hoNDArray<T> prevR(RO, E1, 1), R(RO, E1, 1), imT(RO, E1, 1), magT(RO, E1, 1), diffR(RO, E1, 1);
        hoNDArray<T> coilMapConv(RO, E1, CHA);
        hoNDArray<T> D(RO, E1, CHA);
        hoNDArray<T> D_sum(1, E1, CHA);
        hoNDArray<T> D_sum_1st_2nd(1, 1, CHA);
        typename realType<T>::Type v, vR, vDiffR;
        T vCha;
        size_t iter;
        long long cha;

        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, D_sum, 0));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D_sum, D_sum_1st_2nd, 1));
        Gadgetron::norm2(D_sum_1st_2nd, v);
        Gadgetron::scal((value_type)1.0 / v, D_sum_1st_2nd);

        Gadgetron::clear(R);
        for (cha = 0; cha<CHA; cha++)
        {
            hoNDArray<T> dataCHA(RO, E1, const_cast<T*>(data.begin()) + cha*RO*E1);
            vCha = D_sum_1st_2nd(cha);
            Gadgetron::axpy(std::conj(vCha), dataCHA, R, R);
        }

        for (iter = 0; iter<iterNum; iter++)
        {
            prevR = R;

            Gadgetron::conjugate(R, R);

            GADGET_CATCH_THROW(Gadgetron::multiply(data, R, coilMap));

            Gadgetron::conv2(coilMap, ker, coilMapConv);

            Gadgetron::multiplyConj(coilMapConv, coilMapConv, D);

            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, R, 2));

            Gadgetron::sqrt(R, R);

            Gadgetron::addEpsilon(R);
            Gadgetron::inv(R, R);

            GADGET_CATCH_THROW(Gadgetron::multiply(coilMapConv, R, coilMap));

            Gadgetron::multiplyConj(data, coilMap, D);
            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, R, 2));

            GADGET_CATCH_THROW(Gadgetron::multiply(coilMap, R, D));

            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, D_sum, 0));
            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D_sum, D_sum_1st_2nd, 1));

            Gadgetron::norm2(D_sum_1st_2nd, v);
            Gadgetron::scal((value_type)1.0 / v, D_sum_1st_2nd);

            Gadgetron::clear(imT);
            for (cha = 0; cha<CHA; cha++)
            {
                hoNDArray<T> coilMapCHA(RO, E1, coilMap.begin() + cha*RO*E1);
                vCha = D_sum_1st_2nd(cha);
                Gadgetron::axpy(std::conj(vCha), coilMapCHA, imT, imT);
            }

            Gadgetron::abs(imT, magT);
            Gadgetron::divide(imT, magT, imT);

            Gadgetron::multiply(R, imT, R);
            Gadgetron::conjugate(imT, imT);
            GADGET_CATCH_THROW(Gadgetron::multiply(coilMap, imT, coilMap));

            Gadgetron::subtract(prevR, R, diffR);
            Gadgetron::norm2(diffR, vDiffR);
            Gadgetron::norm2(R, vR);

            if (vDiffR / vR < thres) break;
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in coil_map_2d_Inati_Iter(...) ... ");
        throw;
    }
}

template EXPORTMRICORE void coil_map_2d_Inati_Iter(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t iterNum, float thres);
template EXPORTMRICORE void coil_map_2d_Inati_Iter(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t iterNum, double thres);

// ------------------------------------------------------------------------

template<typename T> 
void coil_map_3d_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t kz, size_t iterNum, typename realType<T>::Type thres)
{
    try
    {
        typedef typename realType<T>::Type value_type;

        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t CHA = data.get_size(3);

        size_t N = data.get_number_of_elements() / (RO*E1*E2*CHA);
        GADGET_CHECK_THROW(N == 1);

        const T* pData = data.begin();

        if (!data.dimensions_equal(&coilMap))
        {
            coilMap = data;
        }

        // create convolution kernel
        hoNDArray<T> ker(ks, ks, kz);
        Gadgetron::fill(&ker, T((value_type)1.0 / (ks*ks*kz)));

        hoNDArray<T> R(RO, E1, E2, 1), imT(RO, E1, E2, 1), magT(RO, E1, E2, 1);
        hoNDArray<T> coilMapConv(RO, E1, E2, CHA);
        hoNDArray<T> D(RO, E1, E2, CHA);
        hoNDArray<T> D_sum(1, CHA);
        typename realType<T>::Type v;
        T vCha;
        size_t iter, cha;

        hoNDArray<T> dataByCha(RO*E1*E2, CHA, const_cast<T*>(data.begin()));
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(data, D_sum, 0));
        Gadgetron::norm2(D_sum, v);
        Gadgetron::scal((value_type)1.0 / v, D_sum);

        Gadgetron::clear(R);
        for (cha = 0; cha<CHA; cha++)
        {
            hoNDArray<T> dataCHA(RO, E1, E2, const_cast<T*>(data.begin()) + cha*RO*E1*E2);
            vCha = D_sum(cha);
            Gadgetron::axpy(std::conj(vCha), dataCHA, R, R);
        }

        for (iter = 0; iter<iterNum; iter++)
        {
            Gadgetron::conjugate(R, R);

            Gadgetron::multiply(data, R, coilMap);

            Gadgetron::conv2(coilMap, ker, coilMapConv);

            Gadgetron::multiplyConj(coilMapConv, coilMapConv, D);

            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, R, 3));

            Gadgetron::sqrt(R, R);

            Gadgetron::addEpsilon(R);
            Gadgetron::inv(R, R);

            Gadgetron::multiply(coilMapConv, R, coilMap);

            Gadgetron::multiplyConj(data, coilMap, D);
            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(D, R, 3));

            Gadgetron::multiply(coilMap, R, D);

            hoNDArray<T> DByCha(RO*E1*E2, CHA, D.begin());
            GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(DByCha, D_sum, 0));

            Gadgetron::norm2(D_sum, v);
            Gadgetron::scal((value_type)1.0 / v, D_sum);

            Gadgetron::clear(imT);
            for (cha = 0; cha<CHA; cha++)
            {
                hoNDArray<T> coilMapCHA(RO, E1, E2, 1, coilMap.begin() + cha*RO*E1*E2);
                vCha = D_sum(cha);
                Gadgetron::axpy(std::conj(vCha), coilMapCHA, imT, imT);
            }

            Gadgetron::abs(imT, magT);
            Gadgetron::divide(imT, magT, imT);

            Gadgetron::multiply(R, imT, R);
            Gadgetron::conjugate(imT, imT);
            Gadgetron::multiply(coilMap, imT, coilMap);
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in coil_map_3d_Inati_Iter(...) ... ");
        throw;
    }
}

template EXPORTMRICORE void coil_map_3d_Inati_Iter(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t kz, size_t iterNum, float thres);
template EXPORTMRICORE void coil_map_3d_Inati_Iter(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t kz, size_t iterNum, double thres);

// ------------------------------------------------------------------------

}
