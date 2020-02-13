//
// Created by dchansen on 2/12/20.
//

#pragma once
namespace {
    using namespace Gadgetron;
    namespace gadgetron_detail {
        constexpr long long NumElementsUseThreading
            = 1024 * 64; // This is probably WAAAY too low. Any real gain for operations this small comes from OMP being
                         // more aggresive with SIMD
        template <class T, class R, class F> inline void omp_transform(const T* t, size_t N, R* r, F f) {

#pragma omp parallel for if (N > gadgetron_detail::NumElementsUseThreading)
            for (long long i = 0; i < (long long)N; i++) {
                r[i] = f(t[i]);
            }
        }
        //
        // Math internal complex types
        // this replaces std::complex<T> with complext<T>
        //
        template <class T> struct mathInternalType { typedef T type; };
        template <class T> struct mathInternalType<std::complex<T>> { typedef Gadgetron::complext<T> type; };

        // --------------------------------------------------------------------------------

        // internal low level function for element-wise addition of two arrays
        template <class T, class S, class BinaryOperator>
        inline void transform_impl(size_t sizeX, size_t sizeY, const T* x, const S* y,
            typename mathReturnType<T, S>::type* r, BinaryOperator op) {

            // cast to internal types
            const typename mathInternalType<T>::type* a
                = reinterpret_cast<const typename mathInternalType<T>::type*>(x);
            const typename mathInternalType<S>::type* b
                = reinterpret_cast<const typename mathInternalType<S>::type*>(y);
            typename mathInternalType<typename mathReturnType<T, S>::type>::type* c
                = reinterpret_cast<typename mathInternalType<typename mathReturnType<T, S>::type>::type*>(r);

            if (sizeY > sizeX) {
                throw std::runtime_error("Add cannot broadcast when the size of x is less than the size of y.");
            }

            if (sizeX == sizeY) {
                // No Broadcasting
                long long loopsize = sizeX;
                long long n;
#ifdef USE_OMP
#pragma omp simd
#endif
                for (n = 0; n < loopsize; n++) {
                    c[n] = op(a[n], b[n]);
                }
            } else {
                // Broadcasting
                long long outerloopsize = sizeX / sizeY;
                long long innerloopsize = sizeX / outerloopsize;
                if (sizeX < gadgetron_detail::NumElementsUseThreading) {
                    // No OMP at All
                    for (long long outer = 0; outer < outerloopsize; outer++) {
                        size_t offset                                = outer * innerloopsize;
                        const typename mathInternalType<T>::type* ai = &a[offset];
                        typename mathInternalType<typename mathReturnType<T, S>::type>::type* ci = &c[offset];
                        for (long long n = 0; n < innerloopsize; n++) {
                            ci[n] = op(a[n], b[n]);
                        }
                    }
                } else if (innerloopsize > gadgetron_detail::NumElementsUseThreading) {
                    // OMP in the inner loop
                    for (long long outer = 0; outer < outerloopsize; outer++) {
                        size_t offset                                = outer * innerloopsize;
                        const typename mathInternalType<T>::type* ai = &a[offset];
                        typename mathInternalType<typename mathReturnType<T, S>::type>::type* ci = &c[offset];
                        long long n;
#ifdef USE_OMP
#pragma omp simd
#endif
                        for (n = 0; n < innerloopsize; n++) {
                            ci[n] = op(ai[n], b[n]);
                        }
                    }
                } else {
                    // OMP in the outer loop
                    long long outer;
#ifdef USE_OMP
#pragma omp simd
#endif
                    for (outer = 0; outer < outerloopsize; outer++) {
                        size_t offset                                = outer * innerloopsize;
                        const typename mathInternalType<T>::type* ai = &a[offset];
                        typename mathInternalType<typename mathReturnType<T, S>::type>::type* ci = &c[offset];
                        for (long long n = 0; n < innerloopsize; n++) {
                            ci[n] = op(ai[n], b[n]);
                        }
                    }
                }
            }
        }

        template <class T, class S, class BinaryFunction>
        void transform_arrays_inplace(hoNDArray<T>& x, const hoNDArray<S>& y, BinaryFunction&& op) {
            if (!compatible_dimensions<T, S>(x, y)) {
                throw std::runtime_error("add: x and y have incompatible dimensions.");
            }
            transform_impl(x.get_number_of_elements(), y.get_number_of_elements(), x.data(), y.data(), x.data(),
                std::forward<BinaryFunction>(op));
        }

        template <class T, class S, class BinaryFunction>
        void transform_arrays(const hoNDArray<T>& x, const hoNDArray<S>& y,
            hoNDArray<typename mathReturnType<T, S>::type>& r, BinaryFunction&& op) {
            // Check the dimensions os x and y for broadcasting.
            if (!compatible_dimensions<T, S>(x, y)) {
                throw std::runtime_error("add: x and y have incompatible dimensions.");
            }

            // Resize r if necessary
            size_t sx = x.get_number_of_elements();
            size_t sy = y.get_number_of_elements();
            size_t sr = r.get_number_of_elements();
            if (sx >= sy) {
                // x is bigger than y or they have the same size
                if (sx != sr) {
                    r.create(x.get_dimensions());
                }
            } else {
                // y is bigger than x
                if (sy != sr) {
                    r.create(y.get_dimensions());
                }
            }

            transform_impl(x.get_number_of_elements(), y.get_number_of_elements(), x.begin(), y.begin(), r.begin(),
                std::forward<BinaryFunction>(op));
        }
    }
}

template <class T, class S>
void Gadgetron::add(const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T, S>::type>& r) {
    gadgetron_detail::transform_arrays(x, y, r, std::plus<>());
}

template <class T, class S>
void Gadgetron::subtract(
    const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T, S>::type>& r) {
    gadgetron_detail::transform_arrays(x, y, r, std::minus<>());
}

template <class T, class S>
void Gadgetron::multiply(
    const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T, S>::type>& r) {
    gadgetron_detail::transform_arrays(x, y, r, std::multiplies<>());
}

template <class T, class S>
void Gadgetron::divide(
    const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T, S>::type>& r) {
    gadgetron_detail::transform_arrays(x, y, r, std::divides<>());
}
template <class T, class S>
void Gadgetron::multiplyConj(
    const hoNDArray<T>& x, const hoNDArray<S>& y, hoNDArray<typename mathReturnType<T, S>::type>& r) {
    gadgetron_detail::transform_arrays(x, y, r, [](auto& a, auto& b) { return a + conj(b); });
}

template <class T, class S> Gadgetron::hoNDArray<T>& Gadgetron::operator+=(hoNDArray<T>& x, const hoNDArray<S>& y) {
    gadgetron_detail::transform_arrays_inplace(x, y, std::plus<>());
    return x;
}
template <class T, class S> Gadgetron::hoNDArray<T>& Gadgetron::operator-=(hoNDArray<T>& x, const hoNDArray<S>& y) {
    gadgetron_detail::transform_arrays_inplace(x, y, std::minus<>());
    return x;
}
template <class T, class S> Gadgetron::hoNDArray<T>& Gadgetron::operator*=(hoNDArray<T>& x, const hoNDArray<S>& y) {
    gadgetron_detail::transform_arrays_inplace(x, y, std::multiplies<>());
    return x;
}
template <class T, class S> Gadgetron::hoNDArray<T>& Gadgetron::operator/=(hoNDArray<T>& x, const hoNDArray<S>& y) {
    gadgetron_detail::transform_arrays_inplace(x, y, std::divides<>());
    return x;
}

template <class T, class S> Gadgetron::hoNDArray<T>& Gadgetron::operator+=(hoNDArray<T>& x, const S& y) {
    long long n;
    size_t N = x.get_number_of_elements();
#pragma omp simd
    for (n = 0; n < (long long)N; ++n) {
        x[n] += y;
    }
    return x;
}

// --------------------------------------------------------------------------------

template <class T, class S> Gadgetron::hoNDArray<T>& Gadgetron::operator-=(hoNDArray<T>& x, const S& y) {

    long long n;

    size_t N = x.get_number_of_elements();

#pragma omp simd
    for (n = 0; n < (long long)N; ++n) {
        x[n] -= y;
    }

    return x;
}

template <class T, class S> Gadgetron::hoNDArray<T>& Gadgetron::operator*=(hoNDArray<T>& x, const S& y) {

    long long n;

    size_t N = x.get_number_of_elements();

#pragma omp simd
    for (n = 0; n < (long long)N; ++n) {
        x[n] *= y;
    }

    return x;
}

// --------------------------------------------------------------------------------

template <class T, class S> Gadgetron::hoNDArray<T>& Gadgetron::operator/=(hoNDArray<T>& x, const S& y) {

    long long n;

    size_t N = x.get_number_of_elements();

#pragma omp simd
    for (n = 0; n < (long long)N; ++n) {
        x[n] /= y;
    }

    return x;
}
