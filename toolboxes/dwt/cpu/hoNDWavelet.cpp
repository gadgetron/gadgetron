
#include "hoNDWavelet.h"

namespace Gadgetron{

template<typename T> 
hoNDWavelet<T>::hoNDWavelet()
{
}

template<typename T>
hoNDWavelet<T>::~hoNDWavelet()
{
}

template<typename T>
void hoNDWavelet<T>::transform(const hoNDArray<T>& in, hoNDArray<T>& out, size_t NDim, size_t level, bool forward)
{
    try
    {
        GADGET_CHECK_THROW(NDim>=1 && NDim<=3);

        std::vector<size_t> dim;
        in.get_dimensions(dim);

        std::vector<size_t> dimOut;
        size_t N = 1;
        size_t NOut = 1;

        size_t ii;
        if (forward)
        {
            dimOut.resize(dim.size() + 1, 1);

            for (ii = 0; ii < NDim; ii++)
            {
                dimOut[ii] = dim[ii];
            }

            N = in.get_size(0);
            NOut = (1 + level) * N;
            if (NDim == 1)
            {
                dimOut[NDim] = 1 + level;
            }
            else if (NDim == 2)
            {
                dimOut[NDim] = 1 + 3 * level;
                N = in.get_size(0) * in.get_size(1);
                NOut = (1 + 3 * level) * N;
            }
            else if (NDim == 3)
            {
                dimOut[NDim] = 1 + 7 * level;
                N = in.get_size(0) * in.get_size(1) * in.get_size(2);
                NOut = (1 + 7 * level) * N;
            }

            for (ii = NDim + 1; ii < dimOut.size(); ii++)
            {
                dimOut[ii] = dim[ii - 1];
            }
        }
        else
        {
            dimOut.resize(dim.size() - 1, 1);

            for (ii = 0; ii < NDim; ii++)
            {
                dimOut[ii] = dim[ii];
            }

            for (ii = NDim + 1; ii < dim.size(); ii++)
            {
                dimOut[ii - 1] = dim[ii];
            }

            NOut = in.get_size(0);
            N = (1 + level) * NOut;
            if (NDim == 1)
            {
            }
            else if (NDim == 2)
            {
                NOut = in.get_size(0) * in.get_size(1);
                N = (1 + 3 * level) * NOut;
            }
            else if (NDim == 3)
            {
                NOut = in.get_size(0) * in.get_size(1) * in.get_size(2);
                N = (1 + 7 * level) * NOut;
            }
        }

        out.create(dimOut);

        if (level == 0)
        {
            memcpy(out.begin(), in.begin(), in.get_number_of_bytes());
            return;
        }

        long long num = in.get_number_of_elements() / N;

        long long n;

        if (NDim == 1)
        {
#pragma omp parallel for default(none) private(n) shared(num, in, out, N, NOut, level, forward) if(num>16)
            for (n = 0; n < num; n++)
            {
                const T* pIn = in.begin() + n*N;
                T* pOut = out.begin() + n*NOut;

                if (forward)
                    this->dwt1D(pIn, pOut, in.get_size(0), level);
                else
                    this->idwt1D(pIn, pOut, in.get_size(0), level);
            }
        }
        else if (NDim == 2)
        {
#pragma omp parallel for default(none) private(n) shared(num, in, out, N, NOut, level, forward) if(num>8)
            for (n = 0; n < num; n++)
            {
                const T* pIn = in.begin() + n*N;
                T* pOut = out.begin() + n*NOut;

                if (forward)
                    this->dwt2D(pIn, pOut, in.get_size(0), in.get_size(1), level);
                else
                    this->idwt2D(pIn, pOut, in.get_size(0), in.get_size(1), level);
            }
        }
        else  if (NDim == 3)
        {
#pragma omp parallel for default(none) private(n) shared(num, in, out, N, NOut, level, forward) if(num>8)
            for (n = 0; n < num; n++)
            {
                const T* pIn = in.begin() + n*N;
                T* pOut = out.begin() + n*NOut;

                if (forward)
                    this->dwt3D(pIn, pOut, in.get_size(0), in.get_size(1), in.get_size(2), level);
                else
                    this->idwt3D(pIn, pOut, in.get_size(0), in.get_size(1), in.get_size(2), level);
            }
        }
    }
    catch (...)
    {
        GADGET_THROW("Errors in hoNDWavelet<T>::transform(...) ... ");
    }
}

// ------------------------------------------------------------
// Instantiation
// ------------------------------------------------------------

template class EXPORTCPUDWT hoNDWavelet<float>;
template class EXPORTCPUDWT hoNDWavelet<double>;
template class EXPORTCPUDWT hoNDWavelet< std::complex<float> >;
template class EXPORTCPUDWT hoNDWavelet< std::complex<double> >;
template class EXPORTCPUDWT hoNDWavelet< complext<float> >;
template class EXPORTCPUDWT hoNDWavelet< complext<double> >;

}
