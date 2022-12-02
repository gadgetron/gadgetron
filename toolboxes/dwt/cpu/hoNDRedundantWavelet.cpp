
#include "hoNDRedundantWavelet.h"
#include <sstream>

namespace Gadgetron{

template<typename T> 
hoNDRedundantWavelet<T>::hoNDRedundantWavelet()
{
}

template<typename T>
hoNDRedundantWavelet<T>::~hoNDRedundantWavelet()
{
}

template<typename T>
void hoNDRedundantWavelet<T>::compute_wavelet_filter(const std::vector<T>& s)
{
    try
    {
        size_t len = s.size();

        fl_d_.resize(len);
        fh_d_.resize(len);
        fl_r_.resize(len);
        fh_r_.resize(len);

        size_t n;

        // decomposition filter
        for (n = 0; n < len; n++)
        {
            fl_d_[n] = s[len - n - 1];
            fh_d_[n] = s[n];
        }

        for (n = 0; n < len; n += 2)
        {
            fh_d_[n] = -fh_d_[n];
        }

        // reconstruction filter
        for (n = 0; n < len; n++)
        {
            fl_r_[n] = s[n] / (T)2;
            fh_r_[n] = s[len - n - 1] / (T)2;
        }

        for (n = 1; n <= len; n += 2)
        {
            fh_r_[n] = -fh_r_[n];
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDWavelet<T>::compute_wavelet_filter(s) ... ");
    }
}

template<typename T>
void hoNDRedundantWavelet<T>::compute_wavelet_filter(const std::string& wav_name)
{
    try
    {
        if (wav_name == "db2")
        {
            s_.resize(4);
            s_[0] = (T)0.482962913144534;
            s_[1] = (T)0.836516303737808;
            s_[2] = (T)0.224143868042013;
            s_[3] = (T)-0.129409522551260;
        }
        else if (wav_name == "db3")
        {
            s_.resize(6);
            s_[0] = (T)0.332670552950083;
            s_[1] = (T)0.806891509311093;
            s_[2] = (T)0.459877502118491;
            s_[3] = (T)-0.135011020010255;
            s_[4] = (T)-0.085441273882027;
            s_[5] = (T)0.035226291885710;
        }
        else if (wav_name == "db4")
        {
            s_.resize(8);
            s_[0] = (T)0.230377813308896;
            s_[1] = (T)0.714846570552915;
            s_[2] = (T)0.630880767929859;
            s_[3] = (T)-0.027983769416859;
            s_[4] = (T)-0.187034811719093;
            s_[5] = (T)0.030841381835560;
            s_[6] = (T)0.032883011666885;
            s_[7] = (T)-0.010597401785069;
        }
        else if (wav_name == "db5")
        {
            s_.resize(10);
            s_[0] = (T)0.160102397974193;
            s_[1] = (T)0.603829269797189;
            s_[2] = (T)0.724308528437772;
            s_[3] = (T)0.138428145901320;
            s_[4] = (T)-0.242294887066382;
            s_[5] = (T)-0.032244869584637;
            s_[6] = (T)0.077571493840046;
            s_[7] = (T)-0.006241490212798;
            s_[8] = (T)-0.012580751999082;
            s_[9] = (T)0.003335725285474;
        }
        else
        {
            std::stringstream ss;
            ss << "Unsupported wavelet name : " << wav_name << std::endl;
            GADGET_THROW(ss.str());
        }

        this->compute_wavelet_filter(s_);
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDWavelet<T>::compute_wavelet_filter(wav_name) ... ");
    }
}

template<typename T>
void hoNDRedundantWavelet<T>::set_wavelet_filter(const std::vector<T>& fl_d, const std::vector<T>& fh_d, const std::vector<T>& fl_r, const std::vector<T>& fh_r)
{
    GADGET_CHECK_THROW(fl_d.size() == fh_d.size());
    GADGET_CHECK_THROW(fl_d.size() == fl_r.size());
    GADGET_CHECK_THROW(fl_d.size() == fh_r.size());

    fl_d_ = fl_d;
    fh_d_ = fh_d;
    fl_r_ = fl_r;
    fh_r_ = fh_r;
}

template<typename T>
void hoNDRedundantWavelet<T>::filter_d(const T* const in, size_t len_in, size_t stride_in, T* out_l, T* out_h, size_t stride_out)
{
    size_t len = fl_d_.size();

    size_t n, m;
    for (n = 0; n < len_in-len+1; n++)
    {
        T vl = 0;
        T vh = 0;
        for (m = 0; m < len; m++)
        {
            vl += in[(n + m)*stride_in] * fl_d_[len - m - 1];
            vh += in[(n + m)*stride_in] * fh_d_[len - m - 1];
        }

        out_l[n*stride_out] = vl;
        out_h[n*stride_out] = vh;
    }

    for (n = len_in - len + 1; n < len_in; n++)
    {
        T vl = 0;
        T vh = 0;
        for (m = 0; m < len; m++)
        {
            size_t k = n + m;
            if (k >= len_in) k -= len_in;
            vl += in[k*stride_in] * fl_d_[len - m - 1];
            vh += in[k*stride_in] * fh_d_[len - m - 1];
        }

        out_l[n*stride_out] = vl;
        out_h[n*stride_out] = vh;
    }
}

template<typename T>
void hoNDRedundantWavelet<T>::filter_r(const T* const in_l, const T* const in_h, size_t len_in, size_t stride_in, T* out, size_t stride_out)
{
    long long len = fl_r_.size();

    long long n, m;
    for (n = len-1; n < (long long)len_in; n++)
    {
        T v = 0;
        for (m = 0; m < len; m++)
        {
            size_t k = (n + m + 1 - len)*stride_in;
            v += (in_l[k] * fl_r_[len - m - 1]) + (in_h[k] * fh_r_[len - m - 1]);
        }

        out[n*stride_out] = v;
    }

    for (n = 0; n < len -1; n++)
    {
        T v = 0;
        for (m = 0; m < len; m++)
        {
            long long k = (n + m + 1 - len);
            if (k < 0) k += len_in;
            v += (in_l[k*stride_in] * fl_r_[len - m - 1]) + (in_h[k*stride_in] * fh_r_[len - m - 1]);
        }

        out[n*stride_out] = v;
    }
}

template<typename T>
void hoNDRedundantWavelet<T>::dwt1D(const T* const in, T* out, size_t RO, size_t level)
{
    memcpy(out, in, sizeof(T)*RO);

    std::vector<T> buf_ro(RO);

    for (size_t n = 0; n < level; n++)
    {
        T* l = out;
        T* h = l + n * RO + RO;

        this->filter_d(l, RO, 1, &buf_ro[0], h, 1);

        memcpy(out, &buf_ro[0], sizeof(T)*RO);
    }
}

template<typename T>
void hoNDRedundantWavelet<T>::idwt1D(const T* const in, T* out, size_t RO, size_t level)
{
    memcpy(out, in, sizeof(T)*RO);

    std::vector<T> buf_ro(RO);

    long long n;
    for (n = (long long)level - 1; n >= 0; n--)
    {
        T* l = out;
        const T* const h = in + n * RO + RO;

        this->filter_r(l, h, RO, 1, &buf_ro[0], 1);
        memcpy(out, &buf_ro[0], sizeof(T)*RO);
    }
}

template<typename T>
void hoNDRedundantWavelet<T>::dwt2D(const T* const in, T* out, size_t RO, size_t E1, size_t level)
{
    memcpy(out, in, sizeof(T)*RO*E1);

    std::vector<T> buf_l(E1), buf_h(E1), buf_ro(RO);

    for (size_t n = 0; n<level; n++)
    {
        T* LH = out + (3 * n + 1)*RO*E1;

        size_t ro, e1;
        // along E1
        for (ro = 0; ro < RO; ro++)
        {
            this->filter_d(out + ro, E1, RO, &buf_l[0], &buf_h[0], 1);

            for (e1 = 0; e1<E1; e1++)
            {
                out[ro + e1*RO] = buf_l[e1];
                LH[ro + e1*RO] = buf_h[e1];
            }
        }

        T* HL = LH + RO*E1;
        T* HH = HL + RO*E1;

        // along RO
        for (e1 = 0; e1<E1; e1++)
        {
            this->filter_d(out + e1*RO, RO, 1, &buf_ro[0], HL + e1*RO, 1);
            memcpy(out + e1*RO, &buf_ro[0], sizeof(T)*RO);

            this->filter_d(LH + e1*RO, RO, 1, &buf_ro[0], HH + e1*RO, 1);
            memcpy(LH + e1*RO, &buf_ro[0], sizeof(T)*RO);
        }
    }
}

template<typename T>
void hoNDRedundantWavelet<T>::idwt2D(const T* const in, T* out, size_t RO, size_t E1, size_t level)
{
    memcpy(out, in, sizeof(T)*RO*E1);

    std::vector<T> buf_ro(RO), buf_e1(E1);

    hoNDArray<T> tmp(RO*E1);
    T* pTmp = tmp.begin();

    long long n;
    for (n = (long long)level - 1; n >= 0; n--)
    {
        const T* const LH = in + (3 * n + 1)*RO*E1;
        const T* const HL = LH + RO*E1;
        const T* const HH = HL + RO*E1;

        size_t ro, e1;
        // along RO
        for (e1 = 0; e1<E1; e1++)
        {
            this->filter_r(out + e1*RO, HL + e1*RO, RO, 1, &buf_ro[0], 1);
            memcpy(out + e1*RO, &buf_ro[0], sizeof(T)*RO);

            this->filter_r(LH + e1*RO, HH + e1*RO, RO, 1, pTmp + e1*RO, 1);
        }

        // along e1
        for (ro = 0; ro<RO; ro++)
        {
            this->filter_r(out + ro, pTmp + ro, E1, RO, &buf_e1[0], 1);

            for (e1 = 0; e1<E1; e1++)
            {
                out[ro + e1*RO] = buf_e1[e1];
            }
        }
    }
}

template<typename T>
void hoNDRedundantWavelet<T>::dwt3D(const T* const in, T* out, size_t RO, size_t E1, size_t E2, size_t level)
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
#pragma omp parallel private(e1) shared(RO, E1, E2, N2D, lll, hll)
            {
                std::vector<T> buf_e2(E2), buf_l(E2), buf_h(E2);
#pragma omp for
                for (e1 = 0; e1 < (long long)E1; e1++)
                {
                    for (size_t ro = 0; ro < RO; ro++)
                    {
                        size_t ind = ro + e1*RO;

                        for (size_t e2 = 0; e2 < E2; e2++)
                        {
                            size_t ind3d = ind + e2*N2D;
                            buf_e2[e2] = lll[ind3d];
                        }

                        this->filter_d(&buf_e2[0], E2, 1, &buf_l[0], &buf_h[0], 1);

                        for (size_t e2 = 0; e2 < E2; e2++)
                        {
                            size_t ind3d = ind + e2*N2D;
                            lll[ind3d] = buf_l[e2];
                            hll[ind3d] = buf_h[e2];
                        }
                    }
                }
            }

            // ------------------------------------------
            // E1
            // ------------------------------------------

            long long e2;

#pragma omp parallel private(e2) shared(RO, E1, E2, N2D, lll, lhl, hll, hhl)
            {
                std::vector<T> buf_e1(E1), buf_l(E1), buf_h(E1);
                std::vector<T> buf_e1_2(E1), buf_l_2(E1), buf_h_2(E1);
#pragma omp for
                for (e2 = 0; e2 < (long long)E2; e2++)
                {
                    size_t ind3D = e2*N2D;
                    for (size_t ro = 0; ro < RO; ro++)
                    {
                        size_t ind;
                        for (size_t e1 = 0; e1 < E1; e1++)
                        {
                            ind = ro + ind3D + e1*RO;
                            buf_e1[e1] = lll[ind];
                            buf_e1_2[e1] = hll[ind];
                        }

                        this->filter_d(&buf_e1[0], E1, 1, &buf_l[0], &buf_h[0], 1);
                        this->filter_d(&buf_e1_2[0], E1, 1, &buf_l_2[0], &buf_h_2[0], 1);

                        for (size_t e1 = 0; e1 < E1; e1++)
                        {
                            ind = ro + ind3D + e1*RO;
                            lll[ind] = buf_l[e1];
                            lhl[ind] = buf_h[e1];

                            hll[ind] = buf_l_2[e1];
                            hhl[ind] = buf_h_2[e1];
                        }
                    }
                }
            }

            // ------------------------------------------
            // RO
            // ------------------------------------------

#pragma omp parallel private(e2) shared(RO, E1, E2, N2D, lll, hll, lhl, hhl, llh, hlh, lhh, hhh)
            {
                std::vector<T> buf_l(RO);
#pragma omp for
                for (e2 = 0; e2 < (long long)E2; e2++)
                {
                    for (size_t e1 = 0; e1 < E1; e1++)
                    {
                        size_t ind3D = e1*RO + e2*N2D;

                        this->filter_d(lll + ind3D, RO, 1, &buf_l[0], llh + ind3D, 1);
                        memcpy(lll + ind3D, &buf_l[0], sizeof(T)*RO);

                        this->filter_d(lhl + ind3D, RO, 1, &buf_l[0], lhh + ind3D, 1);
                        memcpy(lhl + ind3D, &buf_l[0], sizeof(T)*RO);

                        this->filter_d(hll + ind3D, RO, 1, &buf_l[0], hlh + ind3D, 1);
                        memcpy(hll + ind3D, &buf_l[0], sizeof(T)*RO);

                        this->filter_d(hhl + ind3D, RO, 1, &buf_l[0], hhh + ind3D, 1);
                        memcpy(hhl + ind3D, &buf_l[0], sizeof(T)*RO);
                    }
                }
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDWavelet<T>::dwt3D(...) ... ");
    }
}

template<typename T>
void hoNDRedundantWavelet<T>::idwt3D(const T* const in, T* out, size_t RO, size_t E1, size_t E2, size_t level)
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

                    this->filter_r(lll + ind3D, llh + ind3D, RO, 1, pLL + ind3D, 1);
                    this->filter_r(lhl + ind3D, lhh + ind3D, RO, 1, pLH + ind3D, 1);
                    this->filter_r(hll + ind3D, hlh + ind3D, RO, 1, pHL + ind3D, 1);
                    this->filter_r(hhl + ind3D, hhh + ind3D, RO, 1, pHH + ind3D, 1);
                }
            }

            // ------------------------------------------
            // E1
            // ------------------------------------------

#pragma omp parallel private(e2) shared(RO, E1, E2, N2D, pLL, pHL, pLH, pHH)
            {
                std::vector<T> buf_l(E1), buf_l_2(E1);
#pragma omp for
                for (e2 = 0; e2 < (long long)E2; e2++)
                {
                    size_t ind3D = e2*N2D;
                    for (size_t ro = 0; ro < RO; ro++)
                    {
                        size_t ind = ro + ind3D;

                        this->filter_r(pLL + ind, pLH + ind, E1, RO, &buf_l[0], 1);
                        this->filter_r(pHL + ind, pHH + ind, E1, RO, &buf_l_2[0], 1);

                        for (size_t e1 = 0; e1 < E1; e1++)
                        {
                            ind = ro + e1*RO + ind3D;
                            pLL[ind] = buf_l[e1];
                            pHL[ind] = buf_l_2[e1];
                        }
                    }
                }
            }

            // ------------------------------------------
            // E2
            // ------------------------------------------

            long long e1;

#pragma omp parallel for private(e1) shared(RO, E1, E2, N2D, pLL,pHL, out) 
            for (e1 = 0; e1<(long long)E1; e1++)
            {
                for (size_t ro = 0; ro<RO; ro++)
                {
                    size_t ind2D = ro + e1*RO;
                    this->filter_r(pLL + ind2D, pHL + ind2D, E2, N2D, out+ind2D, N2D);
                }
            }
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

template class hoNDRedundantWavelet<float>;
template class hoNDRedundantWavelet<double>;
template class hoNDRedundantWavelet< std::complex<float> >;
template class hoNDRedundantWavelet< std::complex<double> >;
template class hoNDRedundantWavelet< complext<float> >;
template class hoNDRedundantWavelet< complext<double> >;

}
