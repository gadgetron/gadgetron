/** \file   hoNDRedundantWavelet.h
    \brief  Implementation REDUDANT wavelet transformation
    \author Hui Xue
*/

#ifndef hoNDRedundantWavelet_H
#define hoNDRedundantWavelet_H

#include "hoNDWavelet.h"

namespace Gadgetron{

    /// support wavelet transform of 1D, 2D and 3D
    template <typename T> class hoNDRedundantWavelet : public hoNDWavelet<T>
    {
    public:

        typedef typename realType<T>::Type value_type;
        typedef hoNDRedundantWavelet<T> Self;

        hoNDRedundantWavelet();
        virtual ~hoNDRedundantWavelet();

        /// these compute_wavelet_filter should be called first before calling transform
        /// transform is NOT thread-safe, since some class member variables are used as buffer for computation

        /// utility function to compute wavelet filter from commonly used wavelet scale functions
        /// wav_name : "db2", "db3", "db4", "db5"
        /// this function call will fill s_ and compute fl_d_, fh_d_, fl_r_, fh_r_
        void compute_wavelet_filter(const std::string& wav_name);
        /// compute wavelet filters from scale function
        void compute_wavelet_filter(const std::vector<T>& s);
        /// set wavelet filters
        void set_wavelet_filter(const std::vector<T>& fl_d, const std::vector<T>& fh_d, const std::vector<T>& fl_r, const std::vector<T>& fh_r);

    protected:

        /// wavelet scale function for decomposition and reconstruction
        std::vector<T> s_;

        /// wavelet high and low pass filter for decomposition and reconstruction
        std::vector<T> fl_d_;
        std::vector<T> fh_d_;
        std::vector<T> fl_r_;
        std::vector<T> fh_r_;

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

        /// perform decomposition filter
        void filter_d(const T* const in, size_t len_in, size_t stride_in, T* out_l, T* out_h, size_t stride_out);
        /// perform reconstruction filter
        void filter_r(const T* const in_l, const T* const in_h, size_t len_in, size_t stride_in, T* out, size_t stride_out);
    };
}

#endif // hoNDRedundantWavelet_H
