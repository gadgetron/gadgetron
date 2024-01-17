/** \file vector_td.h
    \brief The class vector_td defines a D-dimensional vector of type T.

    The class vector_td defines a D-dimensional vector of type T.
    It is used in the Gadgetron to represent short vectors.
    I.e. it is purposedly templetated with dimensionality D as type unsigned int instead of size_t.
    For larger vectors consider using the NDArray class instead (or a std::vector).
    The vector_td class can be used on both the cpu and gpu.
    The accompanying headers vector_td_opeators.h and vector_td_utilities.h define most of the functionality.
    Note that vector_td should not be used to represent complex numbers. For that we provide the custom class complext
   instead.
*/

#pragma once

#include "core_defines.h"

#include <algorithm>
#include <stdlib.h> // for size_t
#include <type_traits>
#ifdef max
#undef max
#endif // max

namespace Gadgetron {

    template <class T, unsigned int D> class vector_td {
    public:
        T vec[D];
        __inline__ vector_td() = default;

        template <typename... X, typename = std::enable_if_t<(sizeof...(X) > 1)>> constexpr __inline__ __host__ __device__ explicit vector_td(X... xs) : vec{ T(xs)... } {}

        __inline__ vector_td(const vector_td& other) = default;

        template <class T2> __inline__ __host__ __device__ explicit vector_td(const vector_td<T2, D>& other) {
            for (unsigned int i = 0; i < D; i++)
                vec[i] = (T)other[i];
        }

        __inline__ __host__ __device__ explicit vector_td(T x) {
            for (unsigned int i = 0; i < D; i++)
                vec[i] = x;
        }

        template <class STATIC_CONTAINER, class SFINAE = std::enable_if_t<STATIC_CONTAINER().size() == D>>
        explicit vector_td(const STATIC_CONTAINER& other) {
            std::copy(other.begin(), other.end(), this->begin());
        }

        template <class TI, typename std::enable_if<(D > 1) && std::is_convertible<TI,T>::value >::type* = nullptr>
        explicit vector_td(TI input[D]) {
            std::copy(input,input+D,vec);

        }
        __inline__ __host__ __device__ T& operator[](size_t i) {
            return vec[i];
        }

        __inline__ __host__ __device__ const T& operator[](size_t i) const {
            return vec[i];
        }

        __inline__ __host__ __device__ T* begin() {
            return vec;
        }
        __inline__ __host__ __device__ const T* begin() const {
            return vec;
        }
        __inline__ __host__ __device__ T* end() {
            return vec + D;
        }
        __inline__ __host__ __device__ const T* end() const {
            return vec + D;
        }

        static constexpr size_t size() {
            return D;
        }
    };

    template <class T, class... ARGS> auto make_vector_td(ARGS&&... args) {
        return vector_td<T, sizeof...(args)>{ std::forward<ARGS>(args)... };
    }

    //
    // Some typedefs for convenience (templated typedefs are not (yet) available in C++)
    //

    template <class REAL, unsigned int D> struct reald { typedef vector_td<REAL, D> Type; };

    template <unsigned int D> struct uintd { typedef vector_td<unsigned int, D> Type; };

    template <unsigned int D> struct uint64d { typedef vector_td<size_t, D> Type; };

    template <unsigned int D> struct intd { typedef vector_td<int, D> Type; };

    template <unsigned int D> struct int64d { typedef vector_td<long long, D> Type; };

    template <unsigned int D> struct floatd { typedef typename reald<float, D>::Type Type; };

    template <unsigned int D> struct doubled { typedef typename reald<double, D>::Type Type; };


    typedef vector_td<unsigned int, 1> uintd1;
    typedef vector_td<unsigned int, 2> uintd2;
    typedef vector_td<unsigned int, 3> uintd3;
    typedef vector_td<unsigned int, 4> uintd4;
    typedef vector_td<unsigned int, 5> uintd5;

    typedef vector_td<size_t, 1> uint64d1;
    typedef vector_td<size_t, 2> uint64d2;
    typedef vector_td<size_t, 3> uint64d3;
    typedef vector_td<size_t, 4> uint64d4;
    typedef vector_td<size_t, 5> uint64d5;

    typedef vector_td<int, 1> intd1;
    typedef vector_td<int, 2> intd2;
    typedef vector_td<int, 3> intd3;
    typedef vector_td<int, 4> intd4;
    typedef vector_td<int, 5> intd5;

    typedef vector_td<long long, 1> int64d1;
    typedef vector_td<long long, 2> int64d2;
    typedef vector_td<long long, 3> int64d3;
    typedef vector_td<long long, 4> int64d4;
    typedef vector_td<long long, 5> int64d5;

    typedef vector_td<float, 1> floatd1;
    typedef vector_td<float, 2> floatd2;
    typedef vector_td<float, 3> floatd3;
    typedef vector_td<float, 4> floatd4;
    typedef vector_td<float, 5> floatd5;

    typedef vector_td<double, 1> doubled1;
    typedef vector_td<double, 2> doubled2;
    typedef vector_td<double, 3> doubled3;
    typedef vector_td<double, 4> doubled4;
    typedef vector_td<double, 5> doubled5;
}

template <class T, unsigned int N> class std::tuple_size<Gadgetron::vector_td<T, N>> : public std::integral_constant<size_t, N> {};

template <std::size_t I, class T, unsigned int N> struct std::tuple_element<I, Gadgetron::vector_td<T, N>> {
    using type = T;
};

template <size_t I, class T, unsigned int N>
    constexpr std::enable_if_t < I<N, T&> get(Gadgetron::vector_td<T, N>& a) noexcept {
    return a[I];
}

template <size_t I, class T, unsigned int N>
    constexpr std::enable_if_t < I<N, const T&> get(const Gadgetron::vector_td<T, N>& a) noexcept {
    return a[I];
}

namespace Gadgetron {
    template <size_t I, class T, unsigned int N>
        constexpr std::enable_if_t < I<N, T&> get(Gadgetron::vector_td<T, N>& a) noexcept {
        return a[I];
    }
    template <size_t I, class T, unsigned int N>
        constexpr std::enable_if_t < I<N, const T&> get(const Gadgetron::vector_td<T, N>& a) noexcept {
        return a[I];
    }
}
