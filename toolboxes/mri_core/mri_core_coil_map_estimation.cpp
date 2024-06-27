
/** \file   mri_core_coil_map_estimation.cpp
    \brief  Implementation MRI coil sensitivity map estimation.
    \author Hui Xue
*/

#include "mri_core_coil_map_estimation.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "complext.h"
#include "GadgetronTimer.h"
#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

namespace Gadgetron
{

template<typename T> 
void coil_map_2d_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t power)
{
    try
    {
        typedef typename realType<T>::Type value_type;
        using std::abs;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long CHA = data.get_size(2);

        long long N = data.get_number_of_elements() / (RO*E1*CHA);
        GADGET_CHECK_THROW(N == 1);

        const T* pData = data.begin();

        if (!data.dimensions_equal(coilMap))
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

        long long e1;

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
                    phaseU1 /= abs(phaseU1);

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

template void coil_map_2d_Inati(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t power);
template void coil_map_2d_Inati(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t power);

template void coil_map_2d_Inati(const hoNDArray< complext<float> >& data, hoNDArray< complext<float> >& coilMap, size_t ks, size_t power);
template void coil_map_2d_Inati(const hoNDArray< complext<double> >& data, hoNDArray< complext<double> >& coilMap, size_t ks, size_t power);
// ------------------------------------------------------------------------

template<typename T> 
void coil_map_3d_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t kz, size_t power)
{
    try
    {
        typedef typename realType<T>::Type value_type;
        using std::abs;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long E2 = data.get_size(2);
        long long CHA = data.get_size(3);

        long long N = data.get_number_of_elements() / (RO*E1*E2*CHA);
        GADGET_CHECK_THROW(N == 1);

        const T* pData = data.begin();

        if (!data.dimensions_equal(coilMap))
        {
            coilMap = data;
        }
        T* pSen = coilMap.begin();

        if (ks % 2 != 1)
        {
            ks++;
        }

        if (kz % 2 != 1)
        {
            kz++;
        }

        size_t kss = ks*ks*kz;
        long long halfKs = (long long)ks / 2;
        long long halfKz = (long long)kz / 2;

        long long e2;

        #pragma omp parallel default(none) private(e2) shared(ks, kz, RO, E1, E2, CHA, pSen, pData, halfKs, halfKz, power, kss)
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
                        if (e2 >= halfKz && e2<E2 - halfKz && e1 >= halfKs && e1<E1 - halfKs && ro >= halfKs && ro<RO - halfKs)
                        {
                            for (cha = 0; cha<CHA; cha++)
                            {
                                const T* pDataCurr = pData + cha*RO*E1*E2;
                                long long ind = 0;
                                for (ke2 = -halfKz; ke2 <= halfKz; ke2++)
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
                                for (ke2 = -halfKz; ke2 <= halfKz; ke2++)
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
                        v1Norm = nrm2(V1);
                        scal((value_type)1.0 / v1Norm, V1);

                        memcpy(DC.begin(), D.begin(), sizeof(T)*kss*CHA);
                        gemm(DH_D, DC, true, D, false);
                        // gemm(DH_D, D, true, D, false);

                        for (po = 0; po<power; po++)
                        {
                            gemm(V, DH_D, false, V1, false);
                            V1 = V;
                            v1Norm = nrm2(V1);
                            scal((value_type)1.0 / v1Norm, V1);
                        }

                        // compute U1
                        gemm(U1, D, false, V1, false);

                        phaseU1 = U1(0, 0);
                        for (po = 1; po<kss; po++)
                        {
                            phaseU1 += U1(po, 0);
                        }
                        phaseU1 /= abs(phaseU1);

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

template void coil_map_3d_Inati(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t kz, size_t power);
template void coil_map_3d_Inati(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t kz, size_t power);

template void coil_map_3d_Inati(const hoNDArray< complext<float> >& data, hoNDArray< complext<float> >& coilMap, size_t ks, size_t kz, size_t power);
template void coil_map_3d_Inati(const hoNDArray< complext<double> >& data, hoNDArray< complext<double> >& coilMap, size_t ks, size_t kz, size_t power);
// ------------------------------------------------------------------------

template<typename T> 
void coil_map_2d_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t iterNum, typename realType<T>::Type thres)
{
    using std::conj;
    try
    {
        typedef typename realType<T>::Type value_type;

        long long RO = data.get_size(0);
        long long E1 = data.get_size(1);
        long long CHA = data.get_size(2);

        long long N = data.get_number_of_elements() / (RO*E1*CHA);
        GADGET_CHECK_THROW(N == 1);

        const T* pData = data.begin();

        if (!data.dimensions_equal(coilMap))
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

        Gadgetron::sum_over_dimension(data, D_sum, 0);
        Gadgetron::sum_over_dimension(D_sum, D_sum_1st_2nd, 1);
        v = Gadgetron::nrm2(D_sum_1st_2nd);
        Gadgetron::scal((value_type)1.0 / v, D_sum_1st_2nd);

        Gadgetron::clear(R);
        for (size_t cha = 0; cha<CHA; cha++)
        {
            hoNDArray<T> dataCHA(RO, E1, const_cast<T*>(data.begin()) + cha*RO*E1);
            T vCha = D_sum_1st_2nd(cha);
            Gadgetron::axpy(conj(vCha), dataCHA, R);
        }

        for (size_t iter = 0; iter<iterNum; iter++)
        {
            prevR = R;

            Gadgetron::conjugate(R, R);

            Gadgetron::multiply(data, R, coilMap);

            Gadgetron::conv2(coilMap, ker, coilMapConv);

            Gadgetron::multiplyConj(coilMapConv, coilMapConv, D);

            Gadgetron::sum_over_dimension(D, R, 2);

            Gadgetron::sqrt(R, R);

            Gadgetron::addEpsilon(R);
            Gadgetron::inv(R, R);

            Gadgetron::multiply(coilMapConv, R, coilMap);

            Gadgetron::multiplyConj(data, coilMap, D);
            Gadgetron::sum_over_dimension(D, R, 2);

            Gadgetron::multiply(coilMap, R, D);

            Gadgetron::sum_over_dimension(D, D_sum, 0);
            Gadgetron::sum_over_dimension(D_sum, D_sum_1st_2nd, 1);

            v = Gadgetron::nrm2(D_sum_1st_2nd);
            Gadgetron::scal((value_type)1.0 / v, D_sum_1st_2nd);

            Gadgetron::clear(imT);
            for (size_t cha = 0; cha<CHA; cha++)
            {
                hoNDArray<T> coilMapCHA(RO, E1, coilMap.begin() + cha*RO*E1);
                T vCha = D_sum_1st_2nd(cha);
                Gadgetron::axpy(conj(vCha), coilMapCHA, imT);
            }

            Gadgetron::abs(imT, magT);
            Gadgetron::divide(imT, magT, imT);

            Gadgetron::multiply(R, imT, R);
            Gadgetron::conjugate(imT, imT);
            GADGET_CATCH_THROW(Gadgetron::multiply(coilMap, imT, coilMap));

            Gadgetron::subtract(prevR, R, diffR);
            vDiffR = Gadgetron::nrm2(diffR);
            vR = Gadgetron::nrm2(R);

            if (vDiffR / vR < thres) break;
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in coil_map_2d_Inati_Iter(...) ... ");
        throw;
    }
}

template void coil_map_2d_Inati_Iter(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t iterNum, float thres);
template void coil_map_2d_Inati_Iter(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t iterNum, double thres);

template void coil_map_2d_Inati_Iter(const hoNDArray< complext<float> >& data, hoNDArray< complext<float> >& coilMap, size_t ks, size_t iterNum, float thres);
template void coil_map_2d_Inati_Iter(const hoNDArray< complext<double> >& data, hoNDArray< complext<double> >& coilMap, size_t ks, size_t iterNum, double thres);
// ------------------------------------------------------------------------

template<typename T> 
void coil_map_3d_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t kz, size_t iterNum, typename realType<T>::Type thres)
{
    using std::conj;
    typedef typename realType<T>::Type value_type;

    size_t RO = data.get_size(0);
    size_t E1 = data.get_size(1);
    size_t E2 = data.get_size(2);
    size_t CHA = data.get_size(3);

    size_t N = data.get_number_of_elements() / (RO*E1*E2*CHA);
    GADGET_CHECK_THROW(N == 1);

    const T* pData = data.begin();

    if (!data.dimensions_equal(coilMap))
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

    hoNDArray<T> dataByCha(RO*E1*E2, CHA, const_cast<T*>(data.begin()));
    Gadgetron::sum_over_dimension(dataByCha, D_sum, 0);
    auto v = Gadgetron::nrm2(D_sum);
    Gadgetron::scal((value_type)1.0 / v, D_sum);

    Gadgetron::clear(R);
    for (size_t cha = 0; cha<CHA; cha++)
    {
        hoNDArray<T> dataCHA(RO, E1, E2, const_cast<T*>(data.begin()) + cha*RO*E1*E2);
        T vCha = D_sum(cha);
        Gadgetron::axpy(vCha, dataCHA, R);
    }

    for (size_t iter = 0; iter<iterNum; iter++)
    {
        Gadgetron::conjugate(R, R);

        Gadgetron::multiply(data, R, coilMap);

        Gadgetron::conv3(coilMap, ker, coilMapConv);

        Gadgetron::multiplyConj(coilMapConv, coilMapConv, D);

        Gadgetron::sum_over_dimension(D, R, 3);

        Gadgetron::sqrt(R, R);

        Gadgetron::addEpsilon(R);
        Gadgetron::inv(R, R);

        Gadgetron::multiply(coilMapConv, R, coilMap);

        Gadgetron::multiplyConj(data, coilMap, D);
        Gadgetron::sum_over_dimension(D, R, 3);

        Gadgetron::multiply(coilMap, R, D);

        hoNDArray<T> DByCha(RO*E1*E2, CHA, D.begin());
        Gadgetron::sum_over_dimension(DByCha, D_sum, 0);

        auto v= Gadgetron::nrm2(D_sum);
        Gadgetron::scal((value_type)1.0 / v, D_sum);

        Gadgetron::clear(imT);
        for (size_t cha = 0; cha<CHA; cha++)
        {
            hoNDArray<T> coilMapCHA(RO, E1, E2, 1, coilMap.begin() + cha*RO*E1*E2);
            T vCha = D_sum(cha);
            Gadgetron::axpy(conj(vCha), coilMapCHA, imT);
        }

        Gadgetron::abs(imT, magT);
        Gadgetron::divide(imT, magT, imT);

        Gadgetron::multiply(R, imT, R);
        Gadgetron::conjugate(imT, imT);
        Gadgetron::multiply(coilMap, imT, coilMap);
    }


}

template void coil_map_3d_Inati_Iter(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t kz, size_t iterNum, float thres);
template void coil_map_3d_Inati_Iter(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t kz, size_t iterNum, double thres);

template void coil_map_3d_Inati_Iter(const hoNDArray< complext<float> >& data, hoNDArray< complext<float> >& coilMap, size_t ks, size_t kz, size_t iterNum, float thres);
template void coil_map_3d_Inati_Iter(const hoNDArray< complext<double> >& data, hoNDArray< complext<double> >& coilMap, size_t ks, size_t kz, size_t iterNum, double thres);
// ------------------------------------------------------------------------

template<typename T> void coil_map_Inati(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t kz, size_t power)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t CHA = data.get_size(3);

        if (!data.dimensions_equal(coilMap))
        {
            coilMap = data;
        }

        if (CHA <= 1)
        {
            GWARN_STREAM("coil_map_Inati, CHA <= 1");
            return;
        }

        size_t num = data.get_number_of_elements() / (RO*E1*E2*CHA);

        long long n;

        if (E2 > 1)
        {
            for (n = 0; n < (long long)num; n++)
            {

                hoNDArray<T> im(RO, E1, E2, CHA, const_cast<T*>(data.begin() + n*RO*E1*E2*CHA));
                hoNDArray<T> cmap(RO, E1, E2, CHA, coilMap.begin() + n*RO*E1*E2*CHA);

                Gadgetron::coil_map_3d_Inati(im, cmap, ks, kz, power);
            }
        }
        else
        {
#ifdef USE_OMP
            int num_procs = omp_get_num_procs();
#pragma omp parallel for default(none) private(n) shared(num, RO, E1, CHA, data, coilMap, ks, power) if(num>num_procs/2)
#endif // USE_OMP
            for (n = 0; n < (long long)num; n++)
            {
                hoNDArray<T> im(RO, E1, CHA, const_cast<T*>(data.begin()) + n*RO*E1*CHA);
                hoNDArray<T> cmap(RO, E1, CHA, coilMap.begin() + n*RO*E1*CHA);

                Gadgetron::coil_map_2d_Inati(im, cmap, ks, power);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in coil_map_Inati(...) ... ")
    }
}

template void coil_map_Inati(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t kz, size_t power);
template void coil_map_Inati(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t kz, size_t power);

template <typename T> hoNDArray<T> coil_map_Inati(const hoNDArray<T>& data, size_t ks, size_t kz, size_t power) {
    auto coilMap = hoNDArray<T>(data.dimensions());
    coil_map_Inati(data, coilMap, ks, kz, power);
    return coilMap;
}
template hoNDArray<std::complex<float>> coil_map_Inati(const hoNDArray< std::complex<float> >& data, size_t ks, size_t kz, size_t power);
template hoNDArray<std::complex<double>> coil_map_Inati(const hoNDArray< std::complex<double> >& data, size_t ks, size_t kz, size_t power);
// ------------------------------------------------------------------------

template<typename T> void coil_map_Inati_Iter(const hoNDArray<T>& data, hoNDArray<T>& coilMap, size_t ks, size_t kz, size_t iterNum, typename realType<T>::Type thres)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t E2 = data.get_size(2);
        size_t CHA = data.get_size(3);

        if (!data.dimensions_equal(coilMap))
        {
            coilMap = data;
        }

        if (CHA <= 1)
        {
            GWARN_STREAM("coil_map_Inati_Iter, CHA <= 1");
            return;
        }

        size_t num = data.get_number_of_elements() / (RO*E1*E2*CHA);

        long long n;

        if (E2 > 1)
        {
            for (n = 0; n < (long long)num; n++)
            {

                hoNDArray<T> im(RO, E1, E2, CHA, const_cast<T*>(data.begin() + n*RO*E1*E2*CHA));
                hoNDArray<T> cmap(RO, E1, E2, CHA, coilMap.begin() + n*RO*E1*E2*CHA);

                Gadgetron::coil_map_3d_Inati_Iter(im, cmap, ks, kz, iterNum, thres);
            }
        }
        else
        {
#ifdef USE_OMP
            int num_procs = omp_get_num_procs();
#pragma omp parallel for default(none) private(n) shared(num, RO, E1, CHA, data, coilMap, ks, iterNum, thres) if(num>num_procs/2)
#endif // USE_OMP
            for (n = 0; n < (long long)num; n++)
            {
                hoNDArray<T> im(RO, E1, CHA, const_cast<T*>(data.begin()) + n*RO*E1*CHA);
                hoNDArray<T> cmap(RO, E1, CHA, coilMap.begin() + n*RO*E1*CHA);

                Gadgetron::coil_map_2d_Inati_Iter(im, cmap, ks, iterNum, thres);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors happened in coil_map_Inati_Iter(...) ... ")
    }
}

template void coil_map_Inati_Iter(const hoNDArray< std::complex<float> >& data, hoNDArray< std::complex<float> >& coilMap, size_t ks, size_t kz, size_t iterNum, float thres);
template void coil_map_Inati_Iter(const hoNDArray< std::complex<double> >& data, hoNDArray< std::complex<double> >& coilMap, size_t ks, size_t kz, size_t iterNum, double thres);


// ------------------------------------------------------------------------

template <typename T>
void coil_combine(const hoNDArray<T>& data, const hoNDArray<T>& coilMap, size_t cha_dim, hoNDArray<T>& combined)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        size_t NDimCoil = coilMap.get_number_of_dimensions();

        GADGET_CHECK_THROW(data.get_number_of_elements() >= coilMap.get_number_of_elements());

        GADGET_CHECK_THROW(NDim>cha_dim);
        GADGET_CHECK_THROW(NDimCoil>cha_dim);

        size_t n, perChaSize = 1;
        for (n = 0; n < cha_dim+1; n++)
        {
            GADGET_CHECK_THROW(data.get_size(n) == coilMap.get_size(n));

            perChaSize *= data.get_size(n);
        }

        size_t perCombinedSize = perChaSize / data.get_size(cha_dim);

        std::vector<size_t> dim;
        data.get_dimensions(dim);

        std::vector<size_t> dimCoil;
        coilMap.get_dimensions(dimCoil);

        std::vector<size_t> dimCombined(dim);
        dimCombined.erase(dimCombined.begin() + cha_dim);
        combined.create(dimCombined);

        size_t N = data.get_size(cha_dim+1);
        size_t coilN = coilMap.get_size(cha_dim + 1);

        size_t num = data.get_number_of_elements() / (perChaSize*N);
        size_t numCoilMap = coilMap.get_number_of_elements() / (perChaSize*coilN);

        std::vector<size_t> dimChaN(cha_dim+2, 1);
        for (n = 0; n < cha_dim+2; n++)
        {
            if (n < NDim)
            {
                dimChaN[n] = dim[n];
            }
        }

        std::vector<size_t> dimCombinedN(dimChaN);
        dimCombinedN[cha_dim] = 1;

        std::vector<size_t> dimCoilMapChaN(dimChaN);
        dimCoilMapChaN[cha_dim+1] = coilN;

        std::vector<size_t> dimCha(cha_dim + 1, 1);
        for (n = 0; n < cha_dim + 1; n++)
        {
            if (n < NDim) dimCha[n] = dim[n];
        }

        std::vector<size_t> dimCombinedChaOne(dimCha);
        dimCombinedChaOne[cha_dim] = 1;

        size_t nn;
        hoNDArray<T> dataTmp(dimChaN);
        hoNDArray<T> combinedCurr(dimCombinedN);

        for (nn = 0; nn < num; nn++)
        {
            hoNDArray<T> dataCurr(dimChaN, const_cast<T*>(data.begin()) + nn*perChaSize*N);

            size_t nn_coil = nn;
            if (nn_coil >= numCoilMap) nn_coil = numCoilMap - 1;

            hoNDArray<T> coilMapCurr(dimCoilMapChaN, const_cast<T*>(coilMap.begin()) + nn_coil*perChaSize*coilN);

            if (coilN == N)
            {
                Gadgetron::multiplyConj(dataCurr, coilMapCurr, dataTmp);
                GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(dataTmp, combinedCurr, cha_dim));

                memcpy(combined.begin() + nn*perCombinedSize*N, combinedCurr.begin(), combinedCurr.get_number_of_bytes());
            }
            else
            {
                long long d;

#pragma omp parallel default(none) private(d)shared(nn, N, coilN, cha_dim, dimCha, dimCombinedChaOne, perChaSize, perCombinedSize, dataCurr, coilMapCurr, combined) if(N>6)
                {
                    hoNDArray<T> dataTmpN(dimCha);
                    hoNDArray<T> combinedCurrN(dimCombinedChaOne);

#pragma omp for 
                    for (d = 0; d < (long long)N; d++)
                    {
                        size_t d_coil = d;
                        if (d_coil >= coilN) d_coil = coilN - 1;

                        hoNDArray<T> dataCurrN(dimCha, dataCurr.begin() + d*perChaSize);
                        hoNDArray<T> coilMapCurrN(dimCha, coilMapCurr.begin() + d_coil*perChaSize);

                        Gadgetron::multiplyConj(dataCurrN, coilMapCurrN, dataTmpN);
                        Gadgetron::sum_over_dimension(dataTmpN, combinedCurrN, cha_dim);

                        memcpy(combined.begin() + nn*perCombinedSize*N + d*perCombinedSize, combinedCurrN.begin(), combinedCurrN.get_number_of_bytes());
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in coil_combine(...) ... ");
    }
}

template void coil_combine(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& coilMap, size_t cha_dim, hoNDArray< std::complex<float> >& combined);
template void coil_combine(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& coilMap, size_t cha_dim, hoNDArray< std::complex<double> >& combined);

template<class T> hoNDArray<T> coil_combine(const hoNDArray< T >& data, const hoNDArray< T >& coilMap, size_t cha_dim){
    auto combined = hoNDArray<T>{};
    coil_combine(data,coilMap,cha_dim,combined);
    return combined;
}

template hoNDArray<std::complex<float>> coil_combine(const hoNDArray< std::complex<float> >& data, const hoNDArray< std::complex<float> >& coilMap, size_t cha_dim);
template hoNDArray<std::complex<double>> coil_combine(const hoNDArray< std::complex<double> >& data, const hoNDArray< std::complex<double> >& coilMap, size_t cha_dim);

namespace {
    template<class REAL, unsigned int D>
    struct coil_algorithm_wrapper {};

    template<class REAL> struct coil_algorithm_wrapper<REAL,2> {
        static hoNDArray<complext<REAL>> estimate_b1_map(const hoNDArray<complext<REAL>>& data){
            hoNDArray<float_complext> output(data.dimensions());
            coil_map_2d_Inati(data,output);
            return output;
        }
    };
    template<class REAL> struct coil_algorithm_wrapper<REAL,3> {
        static hoNDArray<complext<REAL>> estimate_b1_map(const hoNDArray<complext<REAL>>& data){
            hoNDArray<float_complext> output(data.dimensions());
            coil_map_3d_Inati_Iter(data,output);
            return output;
        }
    };
}


template<class REAL, unsigned int D>
hoNDArray<complext<REAL>> estimate_b1_map(const hoNDArray<complext<REAL>>& data) {
    GadgetronTimer timer("Estimate_b1_map");
    return std::move(coil_algorithm_wrapper<REAL,D>::estimate_b1_map(data));
}


template hoNDArray<complext<float>> estimate_b1_map<float,2>(const hoNDArray<float_complext>& data);
template hoNDArray<complext<float>> estimate_b1_map<float,3>(const hoNDArray<float_complext>& data);




}
