/** \file   hoNDWavelet.h
    \brief  Implementation REDUDANT wavelet transformation
    \author Hui Xue
*/

#ifndef hoNDWavelet_H
#define hoNDWavelet_H

#include "hoNDArray.h"

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron{

    /// support wavelet transform of 1D, 2D and 3D
    template <typename T> class hoNDWavelet
    {
    public:

        typedef typename realType<T>::Type value_type;
        typedef hoNDWavelet<T> Self;

        hoNDWavelet();
        virtual ~hoNDWavelet();

        /// Calculates the forward/inverse wavelet transform
        /// in: input array [RO E1 E2 ...]
        /// out: output array of wavelet coefficients
        /// level: number of wavelet transformation levels
        /// NDim: the number of transformation dimensions
        /// if NDim==1, 1D transformation is performed on the first dimension, out will have the size [RO 1+level E1 E2 ...]
        /// if NDim==2, 2D transformation is performed on the first two dimensions, out will have the size [RO E1 1+3*level E2 ...]
        /// if NDim==3, 3D transformation is performed on the first three dimensions, out will have the size [RO E1 E2 1+7*level ...]
        /// if forward==false, the role of in and out is switched and inverse wavelet transform is performed
        virtual void transform(const hoNDArray<T>& in, hoNDArray<T>& out, size_t NDim, size_t level, bool forward);

    protected:

        /// implementation for 1D dwt and idwt
        /// out: [RO 1+level] array
        virtual void dwt1D(const T* const in, T* out, size_t RO, size_t level) = 0;
        /// in: [RO 1+level] array
        virtual void idwt1D(const T* const in, T* out, size_t RO, size_t level) = 0;

        /// implementation for 2D dwt and idwt
        /// out: [RO 1+3*level] array
        virtual void dwt2D(const T* const in, T* out, size_t RO, size_t E1, size_t level) = 0;
        /// in: [RO 1+3*level] array
        virtual void idwt2D(const T* const in, T* out, size_t RO, size_t E1, size_t level) = 0;

        /// implementation for 2D dwt and idwt
        /// out: [RO 1+7*level] array
        virtual void dwt3D(const T* const in, T* out, size_t RO, size_t E1, size_t E2, size_t level) = 0;
        /// in: [RO 1+7*level] array
        virtual void idwt3D(const T* const in, T* out, size_t RO, size_t E1, size_t E2, size_t level) = 0;
    };
}

#endif // hoNDWavelet_H
