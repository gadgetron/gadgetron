
#include "hoNDHarrWavelet.h"

namespace Gadgetron{

template<typename T> 
hoNDHarrWavelet<T>::hoNDHarrWavelet()
{
}

template<typename T>
hoNDHarrWavelet<T>::~hoNDHarrWavelet()
{
}

template<typename T>
void hoNDHarrWavelet<T>::dwt1D(const T* const in, T* out, size_t RO, size_t level)
{
    memcpy(out, in, sizeof(T)*RO);

    for (size_t n = 0; n < level; n++)
    {
        T* l = out;
        T* h = l + n * RO + RO;

        T v1 = l[0];
        size_t ro;
        for (ro = 0; ro < RO - 1; ro++)
        {
            const T& t = l[ro + 1];
            h[ro] = l[ro] - t;
            l[ro] += t;
        }

        // periodic boundary condition
        h[RO - 1] = l[RO - 1] - v1;
        l[RO - 1] += v1;

        for (ro = 0; ro < RO; ro++)
        {
            l[ro] *= (value_type)(0.5);
            h[ro] *= (value_type)(0.5);
        }
    }
}

template<typename T>
void hoNDHarrWavelet<T>::idwt1D(const T* const in, T* out, size_t RO, size_t level)
{
    memcpy(out, in, sizeof(T)*RO);

    long long n;
    for (n = (long long)level - 1; n >= 0; n--)
    {
        T* l = out;
        const T* const h = in + n * RO + RO;

        T v1 = l[RO - 1];

        size_t ro;
        for (ro = RO - 1; ro > 0; ro--)
        {
            l[ro] = (l[ro] + l[ro - 1]) + (h[ro] - h[ro - 1]);
        }

        l[0] = (l[0] + v1) + (h[0] - h[RO - 1]);

        for (ro = 0; ro < RO; ro++)
        {
            l[ro] *= (value_type)(0.5);
        }
    }
}

template<typename T>
void hoNDHarrWavelet<T>::dwt2D(const T* const in, T* out, size_t RO, size_t E1, size_t level)
{
    value_type scaleFactor = 0.5;

    memcpy(out, in, sizeof(T)*RO*E1);

    for (size_t n = 0; n<level; n++)
    {
        T* LH = out + (3 * n + 1)*RO*E1;

        long long ro;
        for (ro = 0; ro<(long long)RO; ro++)
        {
            T v1 = out[ro];

            long long ii = ro, e1;
            for (e1 = 0; e1<(long long)E1 - 1; e1++)
            {
                LH[ii] = out[ii] - out[ii + RO];
                out[ii] += out[ii + RO];
                ii += RO;
            }

            LH[ii] = out[ii] - v1;
            out[ii] += v1;
        }

        this->apply_harr_scal(RO*E1, scaleFactor, out);
        this->apply_harr_scal(RO*E1, scaleFactor, LH);

        T* HL = LH + RO*E1;
        T* HH = HL + RO*E1;

        long long e1;
        for (e1 = 0; e1<(long long)E1; e1++)
        {
            T v1 = out[e1*RO];
            T v2 = LH[e1*RO];

            size_t ii = e1*RO;
            for (long long ro = 0; ro<(long long)RO - 1; ro++)
            {
                HH[ii] = LH[ii] - LH[ii + 1];
                LH[ii] += LH[ii + 1];

                HL[ii] = out[ii] - out[ii + 1];
                out[ii] += out[ii + 1];

                ii++;
            }

            HH[ii] = LH[ii] - v2;
            LH[ii] += v2;

            HL[ii] = out[ii] - v1;
            out[ii] += v1;
        }

        this->apply_harr_scal(RO*E1, scaleFactor, out);
        this->apply_harr_scal(RO*E1, scaleFactor, LH);
        this->apply_harr_scal(RO*E1, scaleFactor, HL);
        this->apply_harr_scal(RO*E1, scaleFactor, HH);
    }
}

template<typename T>
void hoNDHarrWavelet<T>::idwt2D(const T* const in, T* out, size_t RO, size_t E1, size_t level)
{
    memcpy(out, in, sizeof(T)*RO*E1);

    hoNDArray<T> tmp(RO*E1);
    T* pTmp = tmp.begin();

    value_type scaleFactor = 0.5;

    long long n;
    for (n = (long long)level - 1; n >= 0; n--)
    {
        const T* const LH = in + (3 * n + 1)*RO*E1;
        const T* const HL = LH + RO*E1;
        const T* const HH = HL + RO*E1;

        long long e1;
        for (e1 = 0; e1<(long long)E1; e1++)
        {
            size_t ii = e1*RO + RO - 1;

            T vLL = out[ii];
            T vLH = LH[ii];
            T vHL = HL[ii];
            T vHH = HH[ii];

            for (long long ro = RO - 1; ro>0; ro--)
            {
                out[ii] += out[ii - 1] + HL[ii] - HL[ii - 1];
                pTmp[ii] = LH[ii] + LH[ii - 1] + HH[ii] - HH[ii - 1];

                ii--;
            }

            out[ii] += vLL + HL[ii] - vHL;
            pTmp[ii] = LH[ii] + vLH + HH[ii] - vHH;
        }

        this->apply_harr_scal(RO*E1, scaleFactor, out);
        this->apply_harr_scal(RO*E1, scaleFactor, pTmp);

        long long ro;
        for (ro = 0; ro<(long long)RO; ro++)
        {
            size_t ii = (E1 - 1)*RO + ro;
            T vLL = out[ii];
            T vLH = pTmp[ii];

            for (long long e1 = E1 - 1; e1>0; e1--)
            {
                out[ii] += pTmp[ii] + out[ii - RO] - pTmp[ii - RO];
                ii -= RO;
            }

            out[ro] += pTmp[ro] + vLL - vLH;
        }

        this->apply_harr_scal(RO*E1, scaleFactor, out);
    }
}

template<typename T>
void hoNDHarrWavelet<T>::dwt3D(const T* const in, T* out, size_t RO, size_t E1, size_t E2, size_t level)
{
    try
    {
        memcpy(out, in, sizeof(T)*RO*E1*E2);

        long long N2D = RO*E1;
        long long N3D = RO*E1*E2;

        // process order E2, E1, RO

        for (size_t n = 0; n<level; n++)
        {
            T* lll = out;
            T* llh = lll + n * 7 * N3D + N3D;
            T* lhl = llh + N3D;
            T* lhh = lhl + N3D;
            T* hll = lhh + N3D;
            T* hlh = hll + N3D;
            T* hhl = hlh + N3D;
            T* hhh = hhl + N3D;

            // ------------------------------------------
            // E2
            // ------------------------------------------
            long long e1;
#pragma omp parallel for default(none) private(e1) shared(RO, E1, E2, N2D, lll, hll)
            for (e1 = 0; e1<(long long)E1; e1++)
            {
                for (size_t ro = 0; ro<RO; ro++)
                {
                    size_t ind = ro + e1*RO;

                    T v1 = lll[ind];
                    for (size_t e2 = 0; e2<E2 - 1; e2++)
                    {
                        size_t ind3d = ind + e2*N2D;

                        const T& t = lll[ind3d + N2D];
                        hll[ind3d] = lll[ind3d] - t;
                        lll[ind3d] += t;
                    }

                    if (E2 > 1)
                    {
                        size_t ind3d = ind + (E2 - 1)*N2D;

                        hll[ind3d] = lll[ind3d] - v1;
                        lll[ind3d] += v1;
                    }
                }
            }

            this->apply_harr_scal(N3D, (value_type)(0.5), lll);
            this->apply_harr_scal(N3D, (value_type)(0.5), hll);

            // ------------------------------------------
            // E1
            // ------------------------------------------

            long long e2;

#pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, N2D, lll, lhl, hll, hhl)
            for (e2 = 0; e2<(long long)E2; e2++)
            {
                size_t ind3D = e2*N2D;
                for (size_t ro = 0; ro<RO; ro++)
                {
                    T v1 = lll[ro + ind3D];
                    T v2 = hll[ro + ind3D];

                    size_t ind;
                    for (size_t e1 = 0; e1<E1 - 1; e1++)
                    {
                        ind = ro + ind3D + e1*RO;

                        const T& t1 = lll[ind + RO];
                        lhl[ind] = lll[ind] - t1;
                        lll[ind] += t1;

                        const T& t2 = hll[ind + RO];
                        hhl[ind] = hll[ind] - t2;
                        hll[ind] += t2;
                    }

                    ind = ro + ind3D + (E1 - 1)*RO;

                    lhl[ind] = lll[ind] - v1;
                    lll[ind] += v1;

                    hhl[ind] = hll[ind] - v2;
                    hll[ind] += v2;
                }
            }

#pragma omp parallel sections
            {
#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), lll);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), lhl);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), hll);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), hhl);
            }

            // ------------------------------------------
            // RO
            // ------------------------------------------

#pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, N2D, lll, hll, lhl, hhl, llh, hlh, lhh, hhh)
            for (e2 = 0; e2<(long long)E2; e2++)
            {
                for (size_t e1 = 0; e1<E1; e1++)
                {
                    size_t ind3D = e1*RO + e2*N2D;

                    T v1 = lll[ind3D];
                    T v2 = lhl[ind3D];
                    T v3 = hll[ind3D];
                    T v4 = hhl[ind3D];

                    size_t ind;
                    for (size_t ro = 0; ro<RO - 1; ro++)
                    {
                        ind = ind3D + ro;

                        const T& t1 = lll[ind + 1];
                        llh[ind] = lll[ind] - t1;
                        lll[ind] += t1;

                        const T& t2 = lhl[ind + 1];
                        lhh[ind] = lhl[ind] - t2;
                        lhl[ind] += t2;

                        const T& t3 = hll[ind + 1];
                        hlh[ind] = hll[ind] - t3;
                        hll[ind] += t3;

                        const T& t4 = hhl[ind + 1];
                        hhh[ind] = hhl[ind] - t4;
                        hhl[ind] += t4;
                    }

                    ind = ind3D + RO - 1;

                    llh[ind] = lll[ind] - v1;
                    lll[ind] += v1;

                    lhh[ind] = lhl[ind] - v2;
                    lhl[ind] += v2;

                    hlh[ind] = hll[ind] - v3;
                    hll[ind] += v3;

                    hhh[ind] = hhl[ind] - v4;
                    hhl[ind] += v4;
                }
            }

#pragma omp parallel sections
            {
#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), lll);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), llh);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), lhl);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), lhh);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), hll);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), hlh);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), hhl);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), hhh);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDWavelet<T>::dwt3D(...) ... ");
    }
}

template<typename T>
void hoNDHarrWavelet<T>::idwt3D(const T* const in, T* out, size_t RO, size_t E1, size_t E2, size_t level)
{
    try
    {
        memcpy(out, in, sizeof(T)*RO*E1*E2);

        long long N2D = RO*E1;
        long long N3D = RO*E1*E2;

        hoNDArray<T> LL(N3D);
        T* pLL = LL.begin();

        hoNDArray<T> HL(N3D);
        T* pHL = HL.begin();

        hoNDArray<T> LH(N3D);
        T* pLH = LH.begin();

        hoNDArray<T> HH(N3D);
        T* pHH = HH.begin();

        long long n;
        for (n = (long long)level - 1; n >= 0; n--)
        {
            T* lll = out;
            const T* const llh = in + n * 7 * N3D + N3D;
            const T* const lhl = llh + N3D;
            const T* const lhh = lhl + N3D;
            const T* const hll = lhh + N3D;
            const T* const hlh = hll + N3D;
            const T* const hhl = hlh + N3D;
            const T* const hhh = hhl + N3D;

            // ------------------------------------------
            // RO
            // ------------------------------------------

            long long e2;

#pragma omp parallel for private(e2) shared(RO, E1, E2, N2D, pLL, pHL, pLH, pHH) 
            for (e2 = 0; e2<(long long)E2; e2++)
            {
                for (size_t e1 = 0; e1<E1; e1++)
                {
                    size_t ind3D = e1*RO + e2*N2D;

                    size_t ind;
                    for (size_t ro = RO - 1; ro>0; ro--)
                    {
                        ind = ind3D + ro;
                        pLL[ind] = (lll[ind] + lll[ind - 1]) + (llh[ind] - llh[ind - 1]);
                        pLH[ind] = (lhl[ind] + lhl[ind - 1]) + (lhh[ind] - lhh[ind - 1]);
                        pHL[ind] = (hll[ind] + hll[ind - 1]) + (hlh[ind] - hlh[ind - 1]);
                        pHH[ind] = (hhl[ind] + hhl[ind - 1]) + (hhh[ind] - hhh[ind - 1]);
                    }

                    ind = ind3D + RO - 1;

                    pLL[ind3D] = (lll[ind3D] + lll[ind]) + (llh[ind3D] - llh[ind]);
                    pLH[ind3D] = (lhl[ind3D] + lhl[ind]) + (lhh[ind3D] - lhh[ind]);
                    pHL[ind3D] = (hll[ind3D] + hll[ind]) + (hlh[ind3D] - hlh[ind]);
                    pHH[ind3D] = (hhl[ind3D] + hhl[ind]) + (hhh[ind3D] - hhh[ind]);
                }
            }

#pragma omp parallel sections
            {
#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), pLL);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), pLH);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), pHL);

#pragma omp section
                this->apply_harr_scal(N3D, (value_type)(0.5), pHH);
            }

            // ------------------------------------------
            // E1
            // ------------------------------------------

#pragma omp parallel for default(none) private(e2) shared(RO, E1, E2, N2D, pLL, pHL, pLH, pHH) 
            for (e2 = 0; e2<(long long)E2; e2++)
            {
                size_t ind3D = e2*N2D;
                for (size_t ro = 0; ro<RO; ro++)
                {
                    size_t ind = ro + (E1 - 1)* RO + ind3D;

                    T v1 = pLL[ind];
                    T v2 = pHL[ind];

                    for (size_t e1 = E1 - 1; e1>0; e1--)
                    {
                        ind = ro + e1*RO + ind3D;
                        pLL[ind] = (pLL[ind] + pLL[ind - RO]) + (pLH[ind] - pLH[ind - RO]);
                        pHL[ind] = (pHL[ind] + pHL[ind - RO]) + (pHH[ind] - pHH[ind - RO]);
                    }

                    ind = ro + ind3D;

                    pLL[ind] = (pLL[ind] + v1) + (pLH[ind] - pLH[ind + (E1 - 1)* RO]);
                    pHL[ind] = (pHL[ind] + v2) + (pHH[ind] - pHH[ind + (E1 - 1)* RO]);
                }
            }

            this->apply_harr_scal(N3D, (value_type)(0.5), pLL);
            this->apply_harr_scal(N3D, (value_type)(0.5), pHL);

            // ------------------------------------------
            // E2
            // ------------------------------------------

            long long e1;

#pragma omp parallel for default(none) private(e1) shared(RO, E1, E2, N2D, pLL,pHL, out) 
            for (e1 = 0; e1<(long long)E1; e1++)
            {
                for (size_t ro = 0; ro<RO; ro++)
                {
                    size_t ind2D = ro + e1*RO;

                    size_t ind;
                    for (size_t e2 = E2 - 1; e2>0; e2--)
                    {
                        ind = ind2D + e2*N2D;
                        out[ind] = (pLL[ind] + pLL[ind - N2D]) + (pHL[ind] - pHL[ind - N2D]);
                    }

                    ind = ind2D + (E2 - 1)*N2D;
                    out[ind2D] = (pLL[ind2D] + pLL[ind]) + (pHL[ind2D] - pHL[ind]);
                }
            }

            this->apply_harr_scal(N3D, (value_type)(0.5), out);
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDWavelet<T>::idwt3D(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class hoNDHarrWavelet<float>;
template class hoNDHarrWavelet<double>;
template class hoNDHarrWavelet< std::complex<float> >;
template class hoNDHarrWavelet< std::complex<double> >;
template class hoNDHarrWavelet< complext<float> >;
template class hoNDHarrWavelet< complext<double> >;

}
