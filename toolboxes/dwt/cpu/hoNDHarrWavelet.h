/** \file   hoNDHarrWavelet.h
    \brief  Implementation REDUDANT harr wavelet transformation
    \author Hui Xue
*/

#ifndef hoNDHarrWavelet_H
#define hoNDHarrWavelet_H

#include "hoNDWavelet.h"

namespace Gadgetron{

    /// support wavelet transform of 1D, 2D and 3D
    template <typename T> class hoNDHarrWavelet : public hoNDWavelet<T>
    {
    public:

        typedef typename realType<T>::Type value_type;
        typedef hoNDHarrWavelet<T> Self;

        hoNDHarrWavelet();
        virtual ~hoNDHarrWavelet();

    protected:

        /// implementation for 1D dwt and idwt
        /// out: [RO 1+level] array
        virtual void dwt1D(const T* const in, T* out, size_t RO, size_t level);
        /// in: [RO 1+level] array
        virtual void idwt1D(const T* const in, T* out, size_t RO, size_t level);

        /// implementation for 2D dwt and idwt
        /// out: [RO 1+3*level] array
        virtual void dwt2D(const T* const in, T* out, size_t RO, size_t E1, size_t level);
        /// in: [RO 1+3*level] array
        virtual void idwt2D(const T* const in, T* out, size_t RO, size_t E1, size_t level);

        /// implementation for 2D dwt and idwt
        /// out: [RO 1+7*level] array
        virtual void dwt3D(const T* const in, T* out, size_t RO, size_t E1, size_t E2, size_t level);
        /// in: [RO 1+7*level] array
        virtual void idwt3D(const T* const in, T* out, size_t RO, size_t E1, size_t E2, size_t level);

        /// utility functions
        template <typename T2> 
        void apply_harr_scal(size_t N, T2 a, T2* x)
        {
            for (size_t n = 0; n < N; n++)
            {
                x[n] *= a;
            }
        }

        template <typename T2>
        void apply_harr_scal(size_t N, T2 a, std::complex<T2>* x)
        {
            for (size_t n = 0; n < N; n++)
            {
                const  std::complex<T2> & c = x[n];
                const T2 re = c.real();
                const T2 im = c.imag();

                reinterpret_cast<T2(&)[2]>(x[n])[0] = re*a;
                reinterpret_cast<T2(&)[2]>(x[n])[1] = im*a;
            }
        }

        template <typename T2>
        void apply_harr_scal(size_t N, T2 a, complext<T2>* x)
        {
            for (size_t n = 0; n < N; n++)
            {
                const complext<T2> & c = x[n];
                const T2 re = c.real();
                const T2 im = c.imag();

                reinterpret_cast<T2(&)[2]>(x[n])[0] = re*a;
                reinterpret_cast<T2(&)[2]>(x[n])[1] = im*a;
            }
        }
    };
}

#endif // hoNDHarrWavelet_H
