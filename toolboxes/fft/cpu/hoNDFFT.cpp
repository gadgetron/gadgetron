
// Include for Visual studio, 'cos reasons.
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <set>

#include "hoMatrix.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_math.h"
#include "hoNDFFT.h"
#include <boost/container/flat_set.hpp>

namespace Gadgetron {

    template <typename T> hoNDFFT<T>* hoNDFFT<T>::instance() {
        if (!instance_)
            instance_ = new hoNDFFT<T>();
        return instance_;
    }

    template <class T> hoNDFFT<T>* hoNDFFT<T>::instance_ = NULL;

    namespace {

        template <class T> struct fftw_types {};

        template <> struct fftw_types<float> {
            using complex                      = fftwf_complex;
            using plan                         = fftwf_plan_s;
            static constexpr auto plan_guru    = fftwf_plan_guru64_dft;
            static constexpr auto execute_dft  = fftwf_execute_dft;
            static constexpr auto destroy_plan = fftwf_destroy_plan;
        };

        template <> struct fftw_types<double> {
            using complex                      = fftw_complex;
            using plan                         = fftw_plan_s;
            static constexpr auto plan_guru    = fftw_plan_guru64_dft;
            static constexpr auto execute_dft  = fftw_execute_dft;
            static constexpr auto destroy_plan = fftw_destroy_plan;
        };
        class FFTLock {
        protected:
            static std::mutex lock;
        };
        std::mutex FFTLock::lock;
        template <class T> class SingleFFTPlan : FFTLock {
        public:
            using FFTWComplex = typename fftw_types<T>::complex;

            SingleFFTPlan(int dimension, const hoNDArray<std::complex<T>>& input, hoNDArray<std::complex<T>>& output,
                bool forward) {
                std::lock_guard<std::mutex> guard(lock);

                const auto& dimensions = input.dimensions();
                size_t stride
                    = std::accumulate(dimensions.begin(), dimensions.begin() + dimension, 1, std::multiplies<>());

                auto fftw_dimensions = fftw_iodim64{ static_cast<ptrdiff_t>(dimensions[dimension]), static_cast<ptrdiff_t>(stride), static_cast<ptrdiff_t>(stride) };

                plan = fftw_types<T>::plan_guru(1, &fftw_dimensions, 0, nullptr, (FFTWComplex*)input.data(),
                    (FFTWComplex*)output.data(), forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
            }

            ~SingleFFTPlan() {
                std::lock_guard<std::mutex> guard(lock);
                fftw_types<T>::destroy_plan(plan);
            }

            void execute(const std::complex<T>* input, std::complex<T>* output) {
                fftw_types<T>::execute_dft(plan, (FFTWComplex*)input, (FFTWComplex*)output);
            }

        private:
            typename fftw_types<T>::plan* plan;
        };

        template <class T> class ContigousFFTPlan : FFTLock {
        public:
            using FFTWComplex = typename fftw_types<T>::complex;
            ContigousFFTPlan(
                int rank, const hoNDArray<std::complex<T>>& input, hoNDArray<std::complex<T>>& output, bool forward) {
                std::lock_guard<std::mutex> guard(lock);

                const auto& dimensions = input.dimensions();

                auto strides = std::vector<size_t>(rank + 1, 1);
                std::partial_sum(
                    dimensions.begin(), dimensions.begin() + rank, strides.begin() + 1, std::multiplies<>());

                auto fftw_dimensions = std::vector<fftw_iodim64>(rank);

                for (int i = 0; i < rank; i++) {
                    fftw_dimensions[i] = { (int64_t)dimensions[i], (int64_t)strides[i], (int64_t)strides[i] };
                }
                plan = fftw_types<T>::plan_guru(rank, fftw_dimensions.data(), 0, nullptr, (FFTWComplex*)input.data(),
                    (FFTWComplex*)output.data(), forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
            }
            ~ContigousFFTPlan() {
                std::lock_guard<std::mutex> guard(lock);
                fftw_types<T>::destroy_plan(plan);
            }

            void execute(const std::complex<T>* input, std::complex<T>* output) {
                fftw_types<T>::execute_dft(plan, (FFTWComplex*)input, (FFTWComplex*)output);
            }

        private:
            typename fftw_types<T>::plan* plan;
        };

        const int num_max_threads = omp_get_max_threads();

        int contigous_rank(const boost::container::flat_set<int>& dimensions) {
            if (!dimensions.count(0))
                return 0;

            int rank = std::distance(std::adjacent_find(dimensions.begin(), dimensions.end(),
                                         [](auto val1, auto val2) { return val1 != val2 - 1; }),
                dimensions.end());
            return rank;
        }

        template <typename T>
        static void contigous_fftn(const hoNDArray<std::complex<T>>& a, hoNDArray<std::complex<T>>& r, int rank,
            bool forward, bool normalize) {

            auto plan = ContigousFFTPlan<T>(rank, a, r, forward);
            size_t batch_size
                = std::accumulate(a.dimensions().begin(), a.dimensions().begin() + rank, 1, std::multiplies<>());
            size_t batches = a.size() / batch_size;

#pragma omp parallel for default(none) shared(plan,  a, r, batches, batch_size)
            for (long long i = 0; i < batches; i++) {

                plan.execute(a.data() + i * batch_size, r.data() + i * batch_size);
            }

            if (normalize)
                r *= T(1) / std::sqrt<T>(batch_size);
        }

        template <typename T>
        static void single_fft(int dimension, const hoNDArray<std::complex<T>>& a, hoNDArray<std::complex<T>>& r,
            bool forward, bool normalize) {
            assert(dimension >= 0);
            auto plan              = SingleFFTPlan<T>(dimension, a, r, forward);
            const auto& dimensions = a.dimensions();
            size_t inner_batches
                = std::accumulate(dimensions.begin(), dimensions.begin() + dimension, 1, std::multiplies<>());
            size_t outer_batches
                = std::accumulate(dimensions.begin() + dimension, dimensions.end(), 1, std::multiplies<>());
            size_t outer_batchsize = inner_batches * dimensions[dimension];

#pragma omp parallel for default(none) shared(plan, a, r , outer_batches, inner_batches, outer_batchsize ) collapse(2)
            for (long long outer = 0; outer < outer_batches; outer++) {
                for (long long inner = 0; inner < inner_batches; inner++) {
                    plan.execute(
                        a.data() + inner + outer * outer_batchsize, r.data() + inner + outer * outer_batchsize);
                }
            }

            if (normalize)
                r *= T(1) / std::sqrt<T>(dimensions[dimension]);
        }
    }

    static inline size_t fftshiftPivot(size_t x) {
        return (x + 1) / 2;
    }

    static inline size_t ifftshiftPivot(size_t x) {
        return x - (x + 1) / 2;
    }

    namespace {
        template <typename T>
        inline void fftshift(const std::complex<T>* a, std::complex<T>* r, size_t stride, size_t x, size_t pivot) {
            std::rotate_copy(a, a + pivot, a + x, r);
        }

        template <typename T, typename... INDICES>
        inline void fftshift(const std::complex<T>* a, std::complex<T>* r, size_t stride, size_t n, size_t pivot,
            size_t n2, INDICES... indices) {
            for (size_t i = 0; i < n; i++) {
                auto line_begin  = a + i * stride;
                size_t new_y     = i < pivot ? i + pivot : i - pivot;
                auto output_line = r + new_y * stride;
                fftshift(line_begin, output_line, stride / n2, n2, indices...);
            }
        }

        template <typename T>
        inline void fftshift(std::complex<T>* a, std::complex<T>* a2, std::vector<std::complex<T>>& buffer,
            size_t stride, size_t ny, size_t pivoty, size_t nx, size_t pivotx) {
            for (size_t iy = 0; iy < (ny + 1) / 2; iy++) {
                auto line_begin = a + iy * stride;
                auto line_pivot = line_begin + pivotx;
                auto line_end   = line_begin + stride;

                auto line_begin2 = a2 + ((iy + pivoty) % ny) * stride;
                auto line_pivot2 = line_begin2 + pivotx;
                auto line_end2   = line_begin2 + stride;
                std::rotate_copy(line_begin2, line_pivot2, line_end2, buffer.begin());
                std::rotate_copy(line_begin, line_pivot, line_end, line_begin2);
                std::copy(buffer.begin(), buffer.end(), line_begin);
            }
        }

        template <typename T, typename... INDICES>
        inline void fftshift(std::complex<T>* a, std::complex<T>* a2, std::vector<std::complex<T>>& buffer,
            size_t stride, size_t n, size_t pivot, size_t n2, INDICES... indices) {

            for (size_t i = 0; i < n; i++) {
                auto line_begin  = a + i * stride;
                size_t new_y     = (i + pivot) % n;
                auto line_begin2 = a2 + new_y * stride;
                fftshift(line_begin, line_begin2, buffer, stride / n2, n2, indices...);
            }
        }

    }

    template <typename T> static void fftshiftPivot1D(std::complex<T>* a, size_t x, size_t n, size_t pivot) {

#pragma omp parallel for shared(n, x, pivot, a) if (n > 256) default(none)
        for (long long counter = 0; counter < (long long)n; counter++) {
            std::rotate(a + counter * x, a + counter * x + pivot, a + x + counter * x);
        }
    }

    template <typename T>
    static void fftshiftPivot1D(const std::complex<T>* a, std::complex<T>* r, size_t x, size_t n, size_t pivot) {

#pragma omp parallel for shared(n, x, pivot, a, r) if (n > 256) default(none)
        for (long long counter = 0; counter < (long long)n; counter++) {
            std::rotate_copy(a + counter * x, a + counter * x + pivot, a + x + counter * x, r + counter * x);
        }
    }

    template <typename T> void hoNDFFT<T>::fftshift1D(hoNDArray<std::complex<T>>& a) {
        size_t x           = a.get_size(0);
        size_t pivot       = fftshiftPivot(x);
        size_t numOfShifts = a.get_number_of_elements() / x;
        fftshiftPivot1D(a.begin(), x, numOfShifts, pivot);
    }

    template <typename T>
    void hoNDFFT<T>::fftshift1D(const hoNDArray<std::complex<T>>& a, hoNDArray<std::complex<T>>& r) {
        if (!r.dimensions_equal(&a)) {
            r = a;
        }

        size_t x           = a.get_size(0);
        size_t pivot       = fftshiftPivot(x);
        size_t numOfShifts = a.get_number_of_elements() / x;
        fftshiftPivot1D(a.begin(), r.begin(), x, numOfShifts, pivot);
    }

    template <typename T> void hoNDFFT<T>::ifftshift1D(hoNDArray<std::complex<T>>& a) {
        size_t x           = a.get_size(0);
        size_t pivot       = ifftshiftPivot(x);
        size_t numOfShifts = a.get_number_of_elements() / x;

        fftshiftPivot1D(a.begin(), x, numOfShifts, pivot);
    }

    template <typename T>
    void hoNDFFT<T>::ifftshift1D(const hoNDArray<std::complex<T>>& a, hoNDArray<std::complex<T>>& r) {
        if (!r.dimensions_equal(&a)) {
            r = a;
        }

        size_t x           = a.get_size(0);
        size_t pivot       = ifftshiftPivot(x);
        size_t numOfShifts = a.get_number_of_elements() / x;

        fftshiftPivot1D(a.begin(), r.begin(), x, numOfShifts, pivot);
    }

    template <typename T>
    static void fftshiftPivot2D(
        const std::complex<T>* a, std::complex<T>* r, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty) {
        if (a == NULL || r == NULL)
            throw std::runtime_error("hoNDFFT::fftshiftPivot2D: void ptr provided");
        assert(a != r);

#pragma omp parallel for shared(a, r, x, y, n, pivotx, pivoty) if (n > 16) default(none)
        for (long long tt = 0; tt < (long long)n; tt++) {
            const std::complex<T>* ac = a + tt * x * y;
            std::complex<T>* rc       = r + tt * x * y;
            fftshift(ac, rc, x, y, pivoty, x, pivotx);
        }
    }

    template <typename T>
    static void fftshiftPivot2D(std::complex<T>* a, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty) {

        if (a == NULL)
            throw std::runtime_error("hoNDFFT::fftshiftPivot2D: void ptr provided");

#pragma omp parallel if (n > 16) default(shared)
        {
            std::vector<std::complex<T>> buffer(x);

#pragma omp for
            for (long long tt = 0; tt < (long long)n; tt++) {
                std::complex<T>* ac = a + tt * x * y;
                fftshift(ac, ac, buffer, x, y, pivoty, x, pivotx);
            }
        }
    }

    template <typename T>
    static inline void fftshift2D(const std::complex<T>* a, std::complex<T>* r, size_t x, size_t y, size_t n) {

        if (a == NULL || r == NULL)
            throw std::runtime_error("hoNDFFT::fftshift2D: void ptr provided");

        size_t pivotx = fftshiftPivot(x);
        size_t pivoty = fftshiftPivot(y);

        fftshiftPivot2D(a, r, x, y, n, pivotx, pivoty);
    }

    template <typename T> inline void hoNDFFT<T>::fftshift2D(hoNDArray<ComplexType>& a) {
        size_t n = a.get_number_of_elements() / (a.get_size(0) * a.get_size(1));
        fftshiftPivot2D(
            a.data(), a.get_size(0), a.get_size(1), n, fftshiftPivot(a.get_size(0)), fftshiftPivot(a.get_size(1)));
    }

    template <typename T>
    inline void hoNDFFT<T>::fftshift2D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        if (!r.dimensions_equal(&a)) {
            r = a;
        }

        size_t n = a.get_number_of_elements() / (a.get_size(0) * a.get_size(1));
        fftshiftPivot2D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), n, fftshiftPivot(a.get_size(0)),
            fftshiftPivot(a.get_size(1)));
    }

    template <typename T> inline void hoNDFFT<T>::ifftshift2D(hoNDArray<ComplexType>& a) {
        size_t n = a.get_number_of_elements() / (a.get_size(0) * a.get_size(1));
        fftshiftPivot2D(
            a.begin(), a.get_size(0), a.get_size(1), n, ifftshiftPivot(a.get_size(0)), ifftshiftPivot(a.get_size(1)));
    }

    template <typename T>
    inline void hoNDFFT<T>::ifftshift2D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        if (!r.dimensions_equal(&a)) {
            r = a;
        }
        size_t n = a.get_number_of_elements() / (a.get_size(0) * a.get_size(1));
        fftshiftPivot2D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), n, ifftshiftPivot(a.get_size(0)),
            ifftshiftPivot(a.get_size(1)));
    }

    template <typename T>
    void fftshiftPivot3D(const std::complex<T>* a, std::complex<T>* r, size_t x, size_t y, size_t z, size_t n,
        size_t pivotx, size_t pivoty, size_t pivotz) {

        if (a == NULL || r == NULL)
            throw std::runtime_error("hoNDFFT::fftshift2D: void ptr provided");

        long long tt;

#pragma omp parallel for private(tt) shared(a, r, x, y, z, n, pivotx, pivoty, pivotz) if (n > 16) default(none)
        for (tt = 0; tt < (long long)n; tt++) {
            const std::complex<T>* ac = a + tt * x * y * z;
            std::complex<T>* rc       = r + tt * x * y * z;
            fftshift(ac, rc, x * y, z, pivotz, y, pivoty, x, pivotx);
        }
    }

    template <typename T>
    void fftshiftPivot3D(
        std::complex<T>* a, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty, size_t pivotz) {

        if (a == NULL)
            throw std::runtime_error("hoNDFFT::fftshiftPivot3D: void ptr provided");

        long long tt;

#pragma omp parallel private(tt)  if (n > 16) default(shared)
        {
            std::vector<std::complex<T>> buffer(x);

#pragma omp for
            for (tt = 0; tt < (long long)n; tt++) {
                std::complex<T>* ac = a + tt * x * y * z;
                fftshift(ac, ac, buffer, x * y, z, pivotz, y, pivoty, x, pivotx);
            }
        }
    }

    template <typename T>
    inline void fftshift3D(const std::complex<T>* a, std::complex<T>* r, size_t x, size_t y, size_t z, size_t n) {

        if (a == NULL || r == NULL)
            throw std::runtime_error("hoNDFFT::fftshift3D: void ptr provided");

        size_t pivotx = fftshiftPivot(x);
        size_t pivoty = fftshiftPivot(y);
        size_t pivotz = fftshiftPivot(z);

        fftshiftPivot3D(a, r, x, y, z, n, pivotx, pivoty, pivotz);
    }

    template <typename T>
    inline void ifftshift3D(const std::complex<T>* a, std::complex<T>* r, size_t x, size_t y, size_t z, size_t n) {

        if (a == NULL || r == NULL)
            throw std::runtime_error("hoNDFFT::ifftshift3D: void ptr provided");

        size_t pivotx = ifftshiftPivot(x);
        size_t pivoty = ifftshiftPivot(y);
        size_t pivotz = ifftshiftPivot(z);

        fftshiftPivot3D(a, r, x, y, z, n, pivotx, pivoty, pivotz);
    }

    template <typename T> inline void fftshift3D(std::complex<T>* a, size_t x, size_t y, size_t z, size_t n) {
        if (a == NULL)
            throw std::runtime_error("hoNDFFT::fftshift3D: void ptr provided");
        size_t pivotx = fftshiftPivot(x);
        size_t pivoty = fftshiftPivot(y);
        size_t pivotz = fftshiftPivot(z);
        fftshiftPivot3D(a, x, y, z, n, pivotx, pivoty, pivotz);
    }

    template <typename T> inline void ifftshift3D(std::complex<T>* a, size_t x, size_t y, size_t z, size_t n) {
        if (a == NULL)
            throw std::runtime_error("hoNDFFT::ifftshift3D: void ptr provided");

        size_t pivotx = ifftshiftPivot(x);
        size_t pivoty = ifftshiftPivot(y);
        size_t pivotz = ifftshiftPivot(z);

        fftshiftPivot3D(a, x, y, z, n, pivotx, pivoty, pivotz);
    }

    template <typename T> inline void hoNDFFT<T>::fftshift3D(hoNDArray<ComplexType>& a) {
        size_t n = a.get_number_of_elements() / (a.get_size(0) * a.get_size(1) * a.get_size(2));

        fftshiftPivot3D(a.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n, fftshiftPivot(a.get_size(0)),
            fftshiftPivot(a.get_size(1)), fftshiftPivot(a.get_size(2)));
    }

    template <typename T>
    inline void hoNDFFT<T>::fftshift3D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        if (!r.dimensions_equal(&a)) {
            r = a;
        }

        size_t n = a.get_number_of_elements() / (a.get_size(0) * a.get_size(1) * a.get_size(2));
        fftshiftPivot3D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n,
            fftshiftPivot(a.get_size(0)), fftshiftPivot(a.get_size(1)), fftshiftPivot(a.get_size(2)));
    }

    template <typename T> inline void hoNDFFT<T>::ifftshift3D(hoNDArray<ComplexType>& a) {
        size_t n = a.get_number_of_elements() / (a.get_size(0) * a.get_size(1) * a.get_size(2));
        fftshiftPivot3D(a.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n, ifftshiftPivot(a.get_size(0)),
            ifftshiftPivot(a.get_size(1)), ifftshiftPivot(a.get_size(2)));
    }

    template <typename T>
    inline void hoNDFFT<T>::ifftshift3D(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        if (!r.dimensions_equal(&a)) {
            r = a;
        }
        size_t n = a.get_number_of_elements() / (a.get_size(0) * a.get_size(1) * a.get_size(2));
        fftshiftPivot3D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n,
            ifftshiftPivot(a.get_size(0)), ifftshiftPivot(a.get_size(1)), ifftshiftPivot(a.get_size(2)));
    }

    // -----------------------------------------------------------------------------------------

    template <typename T> inline void hoNDFFT<T>::fft1(hoNDArray<ComplexType>& a) {
        contigous_fftn(a, a, 1, true, true);
    }

    template <typename T> inline void hoNDFFT<T>::ifft1(hoNDArray<ComplexType>& a) {
        contigous_fftn(a, a, 1, false, true);
    }

    template <typename T> inline void hoNDFFT<T>::fft1(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        if (!r.dimensions_equal(&a)) {
            r.create(a.dimensions());
        }

        contigous_fftn(a, r, 1, true, true);
    }

    template <typename T> inline void hoNDFFT<T>::ifft1(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        if (!r.dimensions_equal(&a)) {
            r.create(a.dimensions());
        }
        contigous_fftn(a, r, 1, false, true);
    }

    template <typename T> inline void hoNDFFT<T>::fft1c(hoNDArray<ComplexType>& a) {
        ifftshift1D(a);
        fft1(a);
        fftshift1D(a);
    }

    template <typename T> inline void hoNDFFT<T>::ifft1c(hoNDArray<ComplexType>& a) {
        ifftshift1D(a);
        ifft1(a);
        fftshift1D(a);
    }

    template <typename T> inline void hoNDFFT<T>::fft1c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        ifftshift1D(a, r);
        fft1(r);
        fftshift1D(r);
    }

    template <typename T> inline void hoNDFFT<T>::ifft1c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        ifftshift1D(a, r);
        ifft1(r);
        fftshift1D(r);
    }

    template <typename T>
    inline void hoNDFFT<T>::fft1c(
        const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf) {
        ifftshift1D(a, r);
        fft1(r, buf);
        fftshift1D(buf, r);
    }

    template <typename T>
    inline void hoNDFFT<T>::ifft1c(
        const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf) {
        ifftshift1D(a, r);
        ifft1(r, buf);
        fftshift1D(buf, r);
    }

    // -----------------------------------------------------------------------------------------

    template <typename T> inline void hoNDFFT<T>::fft2(hoNDArray<ComplexType>& a) {
        contigous_fftn(a, a, 2, true, true);
    }

    template <typename T> inline void hoNDFFT<T>::ifft2(hoNDArray<ComplexType>& a) {
        contigous_fftn(a, a, 2, false, true);
    }

    template <typename T> inline void hoNDFFT<T>::fft2(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        if (!r.dimensions_equal(&a)) {
            r.create(a.dimensions());
        }

        contigous_fftn(a, r, 2, true, true);
    }

    template <typename T> inline void hoNDFFT<T>::ifft2(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        if (!r.dimensions_equal(&a)) {
            r.create(a.dimensions());
        }

        contigous_fftn(a, r, 2, false, true);
    }

    template <typename T> inline void hoNDFFT<T>::fft2c(hoNDArray<ComplexType>& a) {
        ifftshift2D(a);
        fft2(a);
        fftshift2D(a);
    }

    template <typename T> inline void hoNDFFT<T>::ifft2c(hoNDArray<ComplexType>& a) {
        ifftshift2D(a);
        ifft2(a);
        fftshift2D(a);
    }

    template <typename T> inline void hoNDFFT<T>::fft2c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        ifftshift2D(a, r);
        fft2(r);
        fftshift2D(r);
    }

    template <typename T> inline void hoNDFFT<T>::ifft2c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        ifftshift2D(a, r);
        ifft2(r);
        fftshift2D(r);
    }

    template <typename T>
    inline void hoNDFFT<T>::fft2c(
        const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf) {
        ifftshift2D(a, r);
        fft2(r, buf);
        fftshift2D(buf, r);
    }

    template <typename T>
    inline void hoNDFFT<T>::ifft2c(
        const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf) {
        ifftshift2D(a, r);
        ifft2(r, buf);
        fftshift2D(buf, r);
    }

    // -----------------------------------------------------------------------------------------

    template <typename T> inline void hoNDFFT<T>::fft3(hoNDArray<ComplexType>& a) {
        contigous_fftn(a, a, 3, true, true);
    }

    template <typename T> inline void hoNDFFT<T>::ifft3(hoNDArray<ComplexType>& a) {
        contigous_fftn(a, a, 3, false, true);
    }

    template <typename T> inline void hoNDFFT<T>::fft3(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        /*r = a;
        return fft3(r);*/
        if (!r.dimensions_equal(&a)) {
            r.create(a.dimensions());
        }

        contigous_fftn(a, r, 3, true, true);
    }

    template <typename T> inline void hoNDFFT<T>::ifft3(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        /*r = a;
        return ifft3(r);*/
        if (!r.dimensions_equal(&a)) {
            r.create(a.dimensions());
        }

        contigous_fftn(a, r, 3, true, true);
    }

    template <typename T> inline void hoNDFFT<T>::fft3c(hoNDArray<ComplexType>& a) {
        ifftshift3D(a);
        fft3(a);
        fftshift3D(a);
    }

    template <typename T> inline void hoNDFFT<T>::ifft3c(hoNDArray<ComplexType>& a) {
        ifftshift3D(a);
        ifft3(a);
        fftshift3D(a);
    }

    template <typename T> inline void hoNDFFT<T>::fft3c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        ifftshift3D(a, r);
        fft3(r);
        fftshift3D(r);
    }

    template <typename T> inline void hoNDFFT<T>::ifft3c(const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r) {
        ifftshift3D(a, r);
        ifft3(r);
        fftshift3D(r);
    }

    template <typename T>
    inline void hoNDFFT<T>::fft3c(
        const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf) {
        ifftshift3D(a, r);
        fft3(r, buf);
        fftshift3D(buf, r);
    }

    template <typename T>
    inline void hoNDFFT<T>::ifft3c(
        const hoNDArray<ComplexType>& a, hoNDArray<ComplexType>& r, hoNDArray<ComplexType>& buf) {
        ifftshift3D(a, r);
        ifft3(r, buf);
        fftshift3D(buf, r);
    }

    template <typename T> void fft1(hoNDArray<std::complex<T>>& a, bool forward) {
        hoNDArray<std::complex<T>> res(a);
        fft1(res, a, forward);
    }

    template <typename T> void fft2(hoNDArray<std::complex<T>>& a, bool forward) {
        hoNDArray<std::complex<T>> res(a);
        fft2(res, a, forward);
    }

    template <typename T> void fft3(hoNDArray<std::complex<T>>& a, bool forward) {
        hoNDArray<std::complex<T>> res(a);
        fft3(res, a, forward);
    }

    template <typename T> void fft1(hoNDArray<std::complex<T>>& a, hoNDArray<std::complex<T>>& r, bool forward) {
        contigous_fftn(a, r, 1, forward, true);
    }

    template <typename T> void fft2(hoNDArray<std::complex<T>>& a, hoNDArray<std::complex<T>>& r, bool forward) {
        contigous_fftn(a, r, 2, forward, true);
    }

    template <typename T> void fft3(hoNDArray<std::complex<T>>& a, hoNDArray<std::complex<T>>& r, bool forward) {
        contigous_fftn(a, r, 3, forward, true);
    }

    template <typename T>
    void fftn(const hoNDArray<std::complex<T>>& a, hoNDArray<std::complex<T>>& r, const std::vector<size_t>& dimensions,
        bool forward) {

        auto dimensions_set = boost::container::flat_set<int>(dimensions.begin(), dimensions.end());
        int c_rank          = contigous_rank(dimensions_set);
        if (c_rank)
            contigous_fftn(a, r, c_rank, forward, false);

        for (auto it = dimensions_set.lower_bound(c_rank); it != dimensions_set.end(); ++it) {
            single_fft(*it,a, r, forward, false);
        }

        size_t batch_size = 1;
        for (auto d : dimensions)
            batch_size *= a.get_size(d);

        T fftRatio = T(1.0 / std::sqrt(T(batch_size)));
        r *= fftRatio;
    }

    // TODO: implement more optimized threading strategy
    inline int get_num_threads_fft1(size_t n0, size_t num) {
        if (num_max_threads == 1)
            return 1;

        if (n0 * num > 1024 * 128) {
            return num_max_threads;
        } else if (n0 * num > 512 * 128) {
            return ((num_max_threads > 8) ? 8 : num_max_threads);
        } else if (n0 * num > 256 * 128) {
            return ((num_max_threads > 4) ? 4 : num_max_threads);
        } else if (n0 * num > 128 * 128) {
            return 2;
        }

        return 1;
    }

    inline int get_num_threads_fft2(size_t n0, size_t n1, size_t num) {
        if (num_max_threads == 1)
            return 1;

        if (n0 * n1 * num > 128 * 128 * 64) {
            return num_max_threads;
        } else if (n0 * n1 * num > 128 * 128 * 32) {
            return ((num_max_threads > 8) ? 8 : num_max_threads);
        } else if (n0 * n1 * num > 128 * 128 * 16) {
            return ((num_max_threads > 4) ? 4 : num_max_threads);
        } else if (n0 * n1 * num > 128 * 128 * 8) {
            return 2;
        }

        return 1;
    }

    inline int get_num_threads_fft3(size_t n0, size_t n1, size_t n2, size_t num) {
        if (num_max_threads == 1)
            return 1;

        if (num >= num_max_threads) {
            return num_max_threads;
        }

        return 1;
    }

    int get_num_threads_fftn(const std::vector<int>& dimensions, size_t num) {

        switch (dimensions.size()) {
        case 1:
            return get_num_threads_fft1(dimensions[0], num);
        case 2:
            return get_num_threads_fft2(dimensions[0], dimensions[1], num);
        case 3:
            return get_num_threads_fft3(dimensions[0], dimensions[1], dimensions[2], num);
        default:
            return std::min<long long>(num_max_threads, num);
        }
    }

    template <typename T> void hoNDFFT<T>::fft(hoNDArray<ComplexType>* input, unsigned int dim_to_transform) {
        single_fft(dim_to_transform, *input, *input, true, true);
    }
    template <typename T> void hoNDFFT<T>::ifft(hoNDArray<ComplexType>* input, unsigned int dim_to_transform) {
        single_fft(dim_to_transform, *input, *input, false, true);
    }
    template <typename T> void hoNDFFT<T>::fft(hoNDArray<ComplexType>* input) {
        contigous_fftn(*input, *input, input->get_number_of_dimensions(), true, true);
    }
    template <typename T> void hoNDFFT<T>::ifft(hoNDArray<ComplexType>* input) {
        contigous_fftn(*input, *input, input->get_number_of_dimensions(), true, true);
    }
    template <typename T> void hoNDFFT<T>::fft(hoNDArray<complext<T>>* input, unsigned int dim_to_transform) {
        fft(reinterpret_cast<hoNDArray<ComplexType>*>(input), dim_to_transform);
    }
    template <typename T> void hoNDFFT<T>::ifft(hoNDArray<complext<T>>* input, unsigned int dim_to_transform) {
        ifft(reinterpret_cast<hoNDArray<ComplexType>*>(input), dim_to_transform);
    }
    template <typename T> void hoNDFFT<T>::fft(hoNDArray<complext<T>>* input) {
        fft(reinterpret_cast<hoNDArray<ComplexType>*>(input));
    }
    template <typename T> void hoNDFFT<T>::ifft(hoNDArray<complext<T>>* input) {
        ifft(reinterpret_cast<hoNDArray<ComplexType>*>(input));
    }


    template <class ComplexType, class ENABLER>
    void FFT::fft(hoNDArray<ComplexType>& data, std::vector<size_t> dimensions) {
        std::sort(dimensions.begin(), dimensions.end());
        if (std::adjacent_find(dimensions.begin(), dimensions.end()) != dimensions.end()) {
            throw std::runtime_error("Duplicate dimensions in list to be transformed");
        }
        fftn(data,data,dimensions,true);
    }

    template <class ComplexType, class ENABLER> void FFT::fft(hoNDArray<ComplexType>& data, size_t dimension) {
        single_fft(dimension,data,data,true,true);
    }
    template <class ComplexType, class ENABLER>
    void FFT::ifft(hoNDArray<ComplexType>& data, std::vector<size_t> dimensions) {
        std::sort(dimensions.begin(), dimensions.end());
        if (std::adjacent_find(dimensions.begin(), dimensions.end()) != dimensions.end()) {
            throw std::runtime_error("Duplicate dimensions in list to be transformed");
        }
        fftn(data,data,dimensions,false);
    }

    template <class ComplexType, class ENABLER> void FFT::ifft(hoNDArray<ComplexType>& data, size_t dimension) {
        single_fft(dimension,data,data,false,true);
    }

    // -----------------------------------------------------------------------------------------

    //
    // Instantiation
    //

    template class EXPORTCPUFFT hoNDFFT<float>;
    template class EXPORTCPUFFT hoNDFFT<double>;


    template void FFT::fft<std::complex<float>>(hoNDArray<std::complex<float>>& data, std::vector<size_t> dimensions);
    template void FFT::fft<std::complex<double>>(hoNDArray<std::complex<double>>& data, std::vector<size_t> dimensions);
    template void FFT::fft<std::complex<float>>(hoNDArray<std::complex<float>>& data, size_t dimensions);
    template void FFT::fft<std::complex<double>>(hoNDArray<std::complex<double>>& data, size_t dimensions);

    template void FFT::ifft<std::complex<float>>(hoNDArray<std::complex<float>>& data, std::vector<size_t> dimensions);
    template void FFT::ifft<std::complex<double>>(hoNDArray<std::complex<double>>& data, std::vector<size_t> dimensions);
    template void FFT::ifft<std::complex<float>>(hoNDArray<std::complex<float>>& data, size_t dimensions);
    template void FFT::ifft<std::complex<double>>(hoNDArray<std::complex<double>>& data, size_t dimensions);
}
