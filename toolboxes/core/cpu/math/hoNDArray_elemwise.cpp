#include "hoNDArray_elemwise.h"
#include "complext.h"
#include "cpp_blas.h"
#include "hoArmadillo.h"
#include "hoNDArray_reductions.h"

#ifdef USE_OMP
#include <omp.h>
#endif

#ifndef lapack_int
#define lapack_int int
#endif // lapack_int

#ifndef lapack_complex_float
#define lapack_complex_float std::complex<float>
#endif // lapack_complex_float

#ifndef lapack_complex_double
#define lapack_complex_double std::complex<double>
#endif // #ifndef lapack_complex_double

#define NumElementsUseThreading 64 * 1024

namespace {

}

namespace Gadgetron {

    template <typename T> void conjugate(const hoNDArray<T>& x, hoNDArray<T>& r) {
        if (r.get_number_of_elements() != x.get_number_of_elements()) {
            r.create(x.dimensions());
        }
        Gadgetron::transform(x,r,[](auto val){return conj(val);});
    }

    template  void conjugate(
        const hoNDArray<std::complex<float>>& x, hoNDArray<std::complex<float>>& r);
    template  void conjugate(
        const hoNDArray<std::complex<double>>& x, hoNDArray<std::complex<double>>& r);

    template  void conjugate(const hoNDArray<complext<float>>& x, hoNDArray<complext<float>>& r);
    template  void conjugate(const hoNDArray<complext<double>>& x, hoNDArray<complext<double>>& r);
    // --------------------------------------------------------------------------------



    template <typename T> void addEpsilon(hoNDArray<T>& x) {
        constexpr auto eps = std::numeric_limits<realType_t<T>>::epsilon();
        using std::abs;
        transform(x,x,[&](auto val){return (abs(val) < eps ) ? val+eps : val; });
    }

    template  void addEpsilon(hoNDArray<float>& x);
    template  void addEpsilon(hoNDArray<double>& x);
    template  void addEpsilon(hoNDArray<std::complex<float>>& x);
    template  void addEpsilon(hoNDArray<std::complex<double>>& x);
    template  void addEpsilon(hoNDArray<complext<float>>& x);
    template  void addEpsilon(hoNDArray<complext<double>>& x);

    // --------------------------------------------------------------------------------

    template <typename T> void argument(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r) {
        if (r.get_number_of_elements() != x.get_number_of_elements()) {
            r.create(x.dimensions());
        }
        using std::arg;
        transform(x,r,[](auto val){return arg(val);});
    }

    template  void argument(const hoNDArray<std::complex<float>>& x, hoNDArray<float>& r);
    template  void argument(const hoNDArray<std::complex<double>>& x, hoNDArray<double>& r);

    template <class T> hoNDArray<realType_t<T>> argument(const hoNDArray<T>& x) {
        using std::arg;
        return Gadgetron::transform(x,[](const auto& v){return arg(v);});
    }

    template  hoNDArray<float> argument(const hoNDArray<std::complex<float>>& x);
    template  hoNDArray<double> argument(const hoNDArray<std::complex<double>>& x);

    // --------------------------------------------------------------------------------

    template <typename T> void inv(const hoNDArray<T>& x, hoNDArray<T>& r) {
        if (!r.dimensions_equal(&x)) {
            r = x;
        }
        transform(x,r,[](auto val){return T(1)/val;});
    }

    template  void inv(const hoNDArray<float>& x, hoNDArray<float>& r);
    template  void inv(const hoNDArray<double>& x, hoNDArray<double>& r);
    template  void inv(const hoNDArray<std::complex<float>>& x, hoNDArray<std::complex<float>>& r);
    template  void inv(const hoNDArray<std::complex<double>>& x, hoNDArray<std::complex<double>>& r);
    template  void inv(const hoNDArray<complext<float>>& x, hoNDArray<complext<float>>& r);
    template  void inv(const hoNDArray<complext<double>>& x, hoNDArray<complext<double>>& r);

    // --------------------------------------------------------------------------------

    template <typename T, typename R> void abs(const hoNDArray<T>& x, hoNDArray<R>& r) {
        if (r.get_number_of_elements() != x.get_number_of_elements()) {
            r.create(x.dimensions());
        }
        transform(x,r,[](auto val){return abs(val);});
    }

    template  void abs(const hoNDArray<float>& x, hoNDArray<float>& r);
    template  void abs(const hoNDArray<double>& x, hoNDArray<double>& r);
    template  void abs(const hoNDArray<std::complex<float>>& x, hoNDArray<float>& r);
    template  void abs(const hoNDArray<std::complex<double>>& x, hoNDArray<double>& r);
    template  void abs(const hoNDArray<complext<float>>& x, hoNDArray<float>& r);
    template  void abs(const hoNDArray<complext<double>>& x, hoNDArray<double>& r);
    template  void abs(const hoNDArray<std::complex<float>>& x, hoNDArray<std::complex<float>>& r);
    template  void abs(const hoNDArray<std::complex<double>>& x, hoNDArray<std::complex<double>>& r);
    template  void abs(const hoNDArray<complext<float>>& x, hoNDArray<complext<float>>& r);
    template  void abs(const hoNDArray<complext<double>>& x, hoNDArray<complext<double>>& r);

    template <class T> hoNDArray<realType_t<T>> abs(const hoNDArray<T>& x) {
        using std::abs;
        return transform(x,[](auto val){return abs(val);}) ;
    }

    template  hoNDArray<float> abs(const hoNDArray<float>& x);
    template  hoNDArray<double> abs(const hoNDArray<double>& x);
    template  hoNDArray<float> abs(const hoNDArray<std::complex<float>>& x);
    template  hoNDArray<double> abs(const hoNDArray<std::complex<double>>& x);
    template  hoNDArray<float> abs(const hoNDArray<complext<float>>& x);
    template  hoNDArray<double> abs(const hoNDArray<complext<double>>& x);

    template <class T> boost::shared_ptr<hoNDArray<typename realType<T>::Type>> abs(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::abs(): Invalid input array");

        boost::shared_ptr<hoNDArray<typename realType<T>::Type>> result(new hoNDArray<typename realType<T>::Type>());
        result->create(x->dimensions());
        abs(*x, *result);
        return result;
    }

    template <class T> void abs_inplace(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::abs_inplace(): Invalid input array");

        abs(*x, *x);
    }

    template <class T> boost::shared_ptr<hoNDArray<typename realType<T>::Type>> abs_square(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::abs_square(): Invalid input array");

        boost::shared_ptr<hoNDArray<typename realType<T>::Type>> result(new hoNDArray<typename realType<T>::Type>());
        result->create(x->dimensions());
        abs(*x, *result);
        multiply(*result, *result, *result);
        return result;
    }

    // --------------------------------------------------------------------------------

    template <typename T> void sqrt(const hoNDArray<T>& x, hoNDArray<T>& r) {
        using std::sqrt;
        if (r.get_number_of_elements() != x.get_number_of_elements()) {
            r.create(x.dimensions());
        }

        size_t N    = x.get_number_of_elements();
        const T* pX = x.begin();
        T* pR       = r.begin();

        long long n;
        for (n = 0; n < (long long)N; n++) {
            pR[n] = sqrt(pX[n]);
        }
    }

    template  void sqrt(const hoNDArray<float>& x, hoNDArray<float>& r);
    template  void sqrt(const hoNDArray<double>& x, hoNDArray<double>& r);
    template  void sqrt(const hoNDArray<std::complex<float>>& x, hoNDArray<std::complex<float>>& r);
    template  void sqrt(const hoNDArray<std::complex<double>>& x, hoNDArray<std::complex<double>>& r);

    template  void sqrt(const hoNDArray<complext<float>>& x, hoNDArray<complext<float>>& r);
    template  void sqrt(const hoNDArray<complext<double>>& x, hoNDArray<complext<double>>& r);

    template <class T> boost::shared_ptr<hoNDArray<T>> sqrt(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::sqrt(): Invalid input array");

        boost::shared_ptr<hoNDArray<T>> result(new hoNDArray<T>());
        result->create(x->dimensions());
        sqrt(*x, *result);
        return result;
    }

    template <class T> void sqrt_inplace(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::sqrt_inplace(): Invalid input array");

        sqrt(*x, *x);
    }

    // --------------------------------------------------------------------------------

    template <class T> boost::shared_ptr<hoNDArray<T>> square(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::square(): Invalid input array");

        boost::shared_ptr<hoNDArray<T>> result(new hoNDArray<T>());
        result->create(x->dimensions());
        /*arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = arma::square(as_arma_col(x));*/
        multiply(*x, *x, *result);
        return result;
    }

    template <class T> void square_inplace(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::square_inplace(): Invalid input array");

        /*arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
        aRes = arma::square(aRes);*/

        multiply(*x, *x, *x);
    }

    // --------------------------------------------------------------------------------

    template <class T> boost::shared_ptr<hoNDArray<T>> reciprocal(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::reciprocal(): Invalid input array");

        /*arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
        ones.ones();*/
        boost::shared_ptr<hoNDArray<T>> result(new hoNDArray<T>());
        result->create(x->dimensions());
        /*arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = ones/as_arma_col(x);*/
        inv(*x, *result);
        return result;
    }

    template <class T> void reciprocal_inplace(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::reciprocal_inplace(): Invalid input array");

        /*arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
        arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
        ones.ones();
        aRes = ones/aRes;*/

        inv(*x, *x);
    }

    template <class T> boost::shared_ptr<hoNDArray<T>> reciprocal_sqrt(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::reciprocal_sqrt(): Invalid input array");

        /*arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
        ones.ones();*/
        boost::shared_ptr<hoNDArray<T>> result(new hoNDArray<T>());
        result->create(x->dimensions());
        /*arma::Col<typename stdType<T>::Type> aRes = as_arma_col(result.get());
        aRes = ones/arma::sqrt(as_arma_col(x));*/

        sqrt(*x, *result);
        inv(*result, *result);
        return result;
    }

    template <class T> void reciprocal_sqrt_inplace(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::reciprocal_sqrt_inplace(): Invalid input array");

        /*arma::Col<typename stdType<T>::Type> ones(x->get_number_of_elements());
        ones.ones();
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(x);
        aRes = ones/arma::sqrt(aRes);*/

        sqrt(*x, *x);
        inv(*x, *x);
    }

    // --------------------------------------------------------------------------------
    namespace {

    template<class T>
    std::complex<T> sgn(std::complex<T> x) {
      if (std::norm(x) <= T(0)) return std::complex<T>(0);
      return (x / std::abs(x));
    }

    }

    template <class T> boost::shared_ptr<hoNDArray<T>> sgn(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::sgn(): Invalid input array");

        boost::shared_ptr<hoNDArray<T>> res(new hoNDArray<T>());
        res->create(x->dimensions());
#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (long long i = 0; i < (long long)res->get_number_of_elements(); i++) {
            res->get_data_ptr()[i] = sgn(x->get_data_ptr()[i]);
        }
        return res;
    }

    template <class T> void sgn_inplace(hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::sgn_inplace(): Invalid input array");

        for (long long i = 0; i < (long long)x->get_number_of_elements(); i++)
            x->get_data_ptr()[i] = sgn(x->get_data_ptr()[i]);
    }


    template<class T> hoNDArray<realType_t<T>> real(const hoNDArray<T>& x){
        return Gadgetron::transform(x,[](const auto& cplx){return real(cplx);});
    }

    template hoNDArray<float> real(const hoNDArray<std::complex<float>>&);
    template hoNDArray<double> real(const hoNDArray<std::complex<double>>&);
    template hoNDArray<float> real(const hoNDArray<complext<float>>&);
    template hoNDArray<double> real(const hoNDArray<complext<double>>&);
    template hoNDArray<float> real(const hoNDArray<float>&);
    template hoNDArray<double> real(const hoNDArray<double>&);
    // --------------------------------------------------------------------------------

    template <class T> boost::shared_ptr<hoNDArray<typename realType<T>::Type>> real(const hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::real(): Invalid input array");

        boost::shared_ptr<hoNDArray<typename realType<T>::Type>> result(new hoNDArray<typename realType<T>::Type>());
        result->create(x->dimensions());
        arma::Col<typename realType<T>::Type> aRes = as_arma_col(*result);
        aRes                                       = arma::real(as_arma_col(*x));
        return result;
    }

    // --------------------------------------------------------------------------------

    template <class T> boost::shared_ptr<hoNDArray<typename realType<T>::Type>> imag(const hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::imag(): Invalid input array");

        boost::shared_ptr<hoNDArray<typename realType<T>::Type>> result(new hoNDArray<typename realType<T>::Type>());
        result->create(x->dimensions());
        arma::Col<typename realType<T>::Type> aRes = as_arma_col(*result);
        aRes                                       = arma::imag(as_arma_col(*x));
        return result;
    }

    template<class T> hoNDArray<realType_t<T>> imag(const hoNDArray<T>& x){
        return Gadgetron::transform(x,[](const auto& cplx){return imag(cplx);});
    }

    template hoNDArray<float> imag(const hoNDArray<std::complex<float>>&);
    template hoNDArray<double> imag(const hoNDArray<std::complex<double>>&);
    template hoNDArray<float> imag(const hoNDArray<complext<float>>&);
    template hoNDArray<double> imag(const hoNDArray<complext<double>>&);
    template hoNDArray<float> imag(const hoNDArray<float>&);
    // --------------------------------------------------------------------------------

    template <class T> boost::shared_ptr<hoNDArray<T>> conj(const hoNDArray<T>* x) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::conj(): Invalid input array");

        boost::shared_ptr<hoNDArray<T>> result(new hoNDArray<T>());
        result->create(x->dimensions());
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(*result);
        aRes                                      = arma::conj(as_arma_col(*x));
        return result;
    }

    // --------------------------------------------------------------------------------

    template <class T> boost::shared_ptr<hoNDArray<T>> real_to_complex(const hoNDArray<typename realType<T>::Type>* x) {
        if (x == 0x0)
            BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_to_complex(): Invalid input array"));

        boost::shared_ptr<hoNDArray<T>> result(new hoNDArray<T>());
        result->create(x->dimensions());
        arma::Col<typename stdType<T>::Type> aRes = as_arma_col(*result);
        aRes                                      = arma::Col<typename stdType<T>::Type>(
            as_arma_col(*x), arma::Col<typename realType<T>::Type>(x->get_number_of_elements()).zeros());
        return result;
    }

    template <class T>
    boost::shared_ptr<hoNDArray<T>> real_imag_to_complex(
        hoNDArray<typename realType<T>::Type>* real, hoNDArray<typename realType<T>::Type>* imag) {
        if (real == 0x0 || imag == 0x0)
            BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_imag_to_complex(): Invalid input array"));

        if (real->get_number_of_elements() != imag->get_number_of_elements())
            BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::real_imag_to_complex(): Invalid input array"));

        boost::shared_ptr<hoNDArray<T>> result(new hoNDArray<T>());
        result->create(real->dimensions());

        T* pRes = result->begin();

        size_t N = real->get_number_of_elements();
        for (size_t n = 0; n < N; n++) {
            pRes[n] = T(real->at(n), imag->at(n));
        }

        return result;
    }

    // --------------------------------------------------------------------------------

    template <class T>
    void real_imag_to_complex(const hoNDArray<typename realType<T>::Type>& real,
        const hoNDArray<typename realType<T>::Type>& imag, hoNDArray<T>& cplx) {
        try {
            GADGET_CHECK_THROW(real.dimensions_equal(&imag));

            if (!cplx.dimensions_equal(&real)) {
                cplx.create(real.dimensions());
            }

            T* pRes                                 = cplx.begin();
            const typename realType<T>::Type* pReal = real.begin();
            const typename realType<T>::Type* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
#pragma omp parallel for private(n) shared(N, pRes, pReal, pImag)
            for (n = 0; n < (long long)N; n++) {
                pRes[n] = T(pReal[n], pImag[n]);
            }
        } catch (...) {
            GADGET_THROW("Errors in real_imag_to_complex(...) ... ");
        }
    }

    template  void real_imag_to_complex(
        const hoNDArray<float>& real, const hoNDArray<float>& imag, hoNDArray<std::complex<float>>& cplx);
    template  void real_imag_to_complex(
        const hoNDArray<double>& real, const hoNDArray<double>& imag, hoNDArray<std::complex<double>>& cplx);

    // --------------------------------------------------------------------------------

    template <class T>
    void complex_to_real_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real,
        hoNDArray<typename realType<T>::Type>& imag) {
        try {
            if (!real.dimensions_equal(&cplx)) {
                real.create(cplx.dimensions());
            }

            if (!imag.dimensions_equal(&cplx)) {
                imag.create(cplx.dimensions());
            }

            const T* pRes                     = cplx.begin();
            typename realType<T>::Type* pReal = real.begin();
            typename realType<T>::Type* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
#pragma omp parallel for default(none) private(n) shared(N, pRes, pReal, pImag)
            for (n = 0; n < (long long)N; n++) {
                pReal[n] = pRes[n].real();
                pImag[n] = pRes[n].imag();
            }
        } catch (...) {
            GADGET_THROW("Errors in complex_to_real_imag(...) ... ");
        }
    }

    template  void complex_to_real_imag(
        const hoNDArray<std::complex<float>>& cplx, hoNDArray<float>& real, hoNDArray<float>& imag);
    template  void complex_to_real_imag(
        const hoNDArray<std::complex<double>>& cplx, hoNDArray<double>& real, hoNDArray<double>& imag);

    void complex_to_real_imag(const hoNDArray<float>& cplx, hoNDArray<float>& real, hoNDArray<float>& imag) {
        try {
            if (!real.dimensions_equal(&cplx)) {
                real.create(cplx.dimensions());
            }

            if (!imag.dimensions_equal(&cplx)) {
                imag.create(cplx.dimensions());
            }

            const float* pRes = cplx.begin();
            float* pReal      = real.begin();
            float* pImag      = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
#pragma omp parallel for default(none) private(n) shared(N, pRes, pReal, pImag)
            for (n = 0; n < (long long)N; n++) {
                pReal[n] = pRes[n];
                pImag[n] = 0;
            }
        } catch (...) {
            GADGET_THROW("Errors in complex_to_real_imag(...) ... ");
        }
    }

    void complex_to_real_imag(const hoNDArray<double>& cplx, hoNDArray<double>& real, hoNDArray<double>& imag) {
        try {
            if (!real.dimensions_equal(&cplx)) {
                real.create(cplx.dimensions());
            }

            if (!imag.dimensions_equal(&cplx)) {
                imag.create(cplx.dimensions());
            }

            const double* pRes = cplx.begin();
            double* pReal      = real.begin();
            double* pImag      = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
#pragma omp parallel for default(none) private(n) shared(N, pRes, pReal, pImag)
            for (n = 0; n < (long long)N; n++) {
                pReal[n] = pRes[n];
                pImag[n] = 0;
            }
        } catch (...) {
            GADGET_THROW("Errors in complex_to_real_imag(...) ... ");
        }
    }

    // --------------------------------------------------------------------------------

    template <class T> void complex_to_real(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real) {
        try {
            if (!real.dimensions_equal(&cplx)) {
                real.create(cplx.dimensions());
            }

            const T* pRes                     = cplx.begin();
            typename realType<T>::Type* pReal = real.begin();

            size_t N = real.get_number_of_elements();

            long long n;
#pragma omp parallel for default(none) private(n) shared(N, pRes, pReal)
            for (n = 0; n < (long long)N; n++) {
                pReal[n] = pRes[n].real();
            }
        } catch (...) {
            GADGET_THROW("Errors in complex_to_real(...) ... ");
        }
    }

    template  void complex_to_real(const hoNDArray<std::complex<float>>& cplx, hoNDArray<float>& real);
    template  void complex_to_real(
        const hoNDArray<std::complex<double>>& cplx, hoNDArray<double>& real);
    template  void complex_to_real(const hoNDArray<complext<float>>& cplx, hoNDArray<float>& real);
    template  void complex_to_real(const hoNDArray<complext<double>>& cplx, hoNDArray<double>& real);

    template <class T> void complex_to_real(const hoNDArray<T>& cplx, hoNDArray<T>& real) {
        try {
            if (!real.dimensions_equal(&cplx)) {
                real.create(cplx.dimensions());
            }

            const T* pRes = cplx.begin();
            T* pReal      = real.begin();

            size_t N = real.get_number_of_elements();

            long long n;
#pragma omp parallel for private(n) shared(N, pRes, pReal)
            for (n = 0; n < (long long)N; n++) {
                pReal[n] = T(pRes[n].real(), 0);
            }
        } catch (...) {
            GADGET_THROW("Errors in complex_to_real(...) ... ");
        }
    }

    template  void complex_to_real(
        const hoNDArray<std::complex<float>>& cplx, hoNDArray<std::complex<float>>& real);
    template  void complex_to_real(
        const hoNDArray<std::complex<double>>& cplx, hoNDArray<std::complex<double>>& real);
    template  void complex_to_real(
        const hoNDArray<complext<float>>& cplx, hoNDArray<complext<float>>& real);
    template  void complex_to_real(
        const hoNDArray<complext<double>>& cplx, hoNDArray<complext<double>>& real);

    template <class T> void complex_to_real(hoNDArray<T>& cplx) {
        try {
            T* pRes = cplx.begin();

            size_t N = cplx.get_number_of_elements();

            long long n;
#pragma omp parallel for private(n) shared(N, pRes)
            for (n = 0; n < (long long)N; n++) {
                pRes[n] = T(pRes[n].real(), 0);
            }
        } catch (...) {
            GADGET_THROW("Errors in complex_to_real(...) ... ");
        }
    }

    template  void complex_to_real(hoNDArray<std::complex<float>>& cplx);
    template  void complex_to_real(hoNDArray<std::complex<double>>& cplx);
    template  void complex_to_real(hoNDArray<complext<float>>& cplx);
    template  void complex_to_real(hoNDArray<complext<double>>& cplx);

    // --------------------------------------------------------------------------------

    template <class T> void complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& imag) {
        try {
            if (!imag.dimensions_equal(&cplx)) {
                imag.create(cplx.dimensions());
            }

            const T* pRes                     = cplx.begin();
            typename realType<T>::Type* pImag = imag.begin();

            size_t N = imag.get_number_of_elements();

            long long n;
#pragma omp parallel for default(none) private(n) shared(N, pRes, pImag)
            for (n = 0; n < (long long)N; n++) {
                pImag[n] = pRes[n].imag();
            }
        } catch (...) {
            GADGET_THROW("Errors in complex_to_imag(...) ... ");
        }
    }

    template  void complex_to_imag(const hoNDArray<std::complex<float>>& cplx, hoNDArray<float>& imag);
    template  void complex_to_imag(
        const hoNDArray<std::complex<double>>& cplx, hoNDArray<double>& imag);

    template <class T> void complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<T>& imag) {
        try {
            if (!imag.dimensions_equal(&cplx)) {
                imag.create(cplx.dimensions());
            }

            const T* pRes = cplx.begin();
            T* pImag      = imag.begin();

            size_t N = imag.get_number_of_elements();

            long long n;
#pragma omp parallel for private(n) shared(N, pRes, pImag)
            for (n = 0; n < (long long)N; n++) {
                pImag[n] = T(0, pRes[n].imag());
            }
        } catch (...) {
            GADGET_THROW("Errors in complex_to_imag(...) ... ");
        }
    }

    template  void complex_to_imag(
        const hoNDArray<std::complex<float>>& cplx, hoNDArray<std::complex<float>>& imag);
    template  void complex_to_imag(
        const hoNDArray<std::complex<double>>& cplx, hoNDArray<std::complex<double>>& imag);

    template <class T> void complex_to_imag(hoNDArray<T>& cplx) {
        try {
            T* pRes = cplx.begin();

            size_t N = cplx.get_number_of_elements();

            long long n;
#pragma omp parallel for private(n) shared(N, pRes)
            for (n = 0; n < (long long)N; n++) {
                pRes[n] = T(pRes[n].real(), 0);
            }
        } catch (...) {
            GADGET_THROW("Errors in complex_to_imag(...) ... ");
        }
    }

    template  void complex_to_imag(hoNDArray<std::complex<float>>& cplx);
    template  void complex_to_imag(hoNDArray<std::complex<double>>& cplx);

    // --------------------------------------------------------------------------------

    template <class T> void real_to_complex(const hoNDArray<typename realType<T>::Type>& real, hoNDArray<T>& cplx) {
        try {
            if (!cplx.dimensions_equal(&real)) {
                cplx.create(real.dimensions());
            }

            const typename realType<T>::Type* pReal = real.begin();
            T* pRes                                 = cplx.begin();

            size_t N = real.get_number_of_elements();

            long long n;
#pragma omp parallel for private(n) shared(N, pRes, pReal)
            for (n = 0; n < (long long)N; n++) {
                pRes[n] = T(pReal[n], 0);
            }
        } catch (...) {
            GADGET_THROW("Errors in real_to_complex(...) ... ");
        }
    }

    template  void real_to_complex(const hoNDArray<float>& real, hoNDArray<std::complex<float>>& cplx);
    template  void real_to_complex(
        const hoNDArray<double>& real, hoNDArray<std::complex<double>>& cplx);

    // --------------------------------------------------------------------------------

    template <typename T> void fill(hoNDArray<T>* x, T val) {
        size_t N = x->get_number_of_elements();
        T* pX    = x->begin();

        long long n;
#pragma omp parallel for default(none) private(n) shared(N, pX, val) if (N > NumElementsUseThreading)
        for (n = 0; n < (long long)N; n++) {
            pX[n] = val;
        }
    }

    template  void fill(hoNDArray<int>* x, int val);
    template  void fill(hoNDArray<float>* x, float val);
    template  void fill(hoNDArray<double>* x, double val);
    template  void fill(hoNDArray<bool>* x, bool val);
    template  void fill(hoNDArray<std::complex<float>>* x, std::complex<float> val);
    template  void fill(hoNDArray<std::complex<double>>* x, std::complex<double> val);
    template  void fill(hoNDArray<complext<float>>* x, complext<float> val);
    template  void fill(hoNDArray<complext<double>>* x, complext<double> val);

    // --------------------------------------------------------------------------------

    template <typename T> void fill(hoNDArray<T>& x, T val) {
        Gadgetron::fill(&x, val);
    }

    template  void fill(hoNDArray<int>& x, int val);
    template  void fill(hoNDArray<float>& x, float val);
    template  void fill(hoNDArray<double>& x, double val);
    template  void fill(hoNDArray<std::complex<float>>& x, std::complex<float> val);
    template  void fill(hoNDArray<std::complex<double>>& x, std::complex<double> val);
    template  void fill(hoNDArray<complext<float>>& x, complext<float> val);
    template  void fill(hoNDArray<complext<double>>& x, complext<double> val);

    // --------------------------------------------------------------------------------

    //
    // TODO:
    // The clamp functions could (probably) be implemented much like we use Thrust for the device versions
    // - i.e. using Armadillo's transform on the array.
    // However this requires a newer version of Armadillo as current Linux distributions provide...
    //

    template <typename T>
    struct hoNDA_clamp //: public thrust::unary_function<T,T>
    {
        hoNDA_clamp(T _min, T _max, T _min_val, T _max_val)
            : min(_min), max(_max), min_val(_min_val), max_val(_max_val) {}
        T operator()(const T& x) const {
            if (x < min)
                return min_val;
            else if (x >= max)
                return max_val;
            else
                return x;
        }
        T min, max;
        T min_val, max_val;
    };

    template <typename T>
    struct hoNDA_clamp<std::complex<T>> //: public thrust::unary_function< std::complex<T>, std::complex<T> >
    {
        hoNDA_clamp(T _min, T _max, std::complex<T> _min_val, std::complex<T> _max_val)
            : min(_min), max(_max), min_val(_min_val), max_val(_max_val) {}
        std::complex<T> operator()(const std::complex<T>& x) const {
            if (real(x) < min)
                return min_val;
            else if (real(x) >= max)
                return max_val;
            else
                return std::complex<T>(real(x));
        }
        T min, max;
        std::complex<T> min_val, max_val;
    };

    template <typename T>
    struct hoNDA_clamp<complext<T>> //: public thrust::unary_function< complext<T>, complext<T> >
    {
        hoNDA_clamp(T _min, T _max, complext<T> _min_val, complext<T> _max_val)
            : min(_min), max(_max), min_val(_min_val), max_val(_max_val) {}
        complext<T> operator()(const complext<T>& x) const {
            if (real(x) < min)
                return min_val;
            else if (real(x) >= max)
                return max_val;
            else
                return complext<T>(real(x));
        }
        T min, max;
        complext<T> min_val, max_val;
    };

    template <class T>
    void clamp(hoNDArray<T>* x, typename realType<T>::Type min, typename realType<T>::Type max, T min_val, T max_val) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::clamp(): Invalid input array");

        hoNDA_clamp<T> functor(min, max, min_val, max_val);
        std::transform(x->begin(), x->end(), x->begin(), functor);
    }

    template <class T> void clamp(hoNDArray<T>* x, typename realType<T>::Type min, typename realType<T>::Type max) {
        clamp(x, min, max, T(min), T(max));
    }

    template <typename T>
    struct hoNDA_clamp_min //: public thrust::unary_function<T,T>
    {
        hoNDA_clamp_min(T _min) : min(_min) {}
        T operator()(const T& x) const {
            if (x < min)
                return min;
            else
                return x;
        }
        T min;
    };

    template <typename T>
    struct hoNDA_clamp_min<std::complex<T>> //: public thrust::unary_function< std::complex<T>, std::complex<T> >
    {
        hoNDA_clamp_min(T _min) : min(_min) {}
        std::complex<T> operator()(const std::complex<T>& x) const {
            if (real(x) < min)
                return std::complex<T>(min);
            else
                return std::complex<T>(real(x));
        }
        T min;
    };

    template <typename T>
    struct hoNDA_clamp_min<complext<T>> //: public thrust::unary_function< complext<T>, complext<T> >
    {
        hoNDA_clamp_min(T _min) : min(_min) {}
        complext<T> operator()(const complext<T>& x) const {
            if (real(x) < min)
                return complext<T>(min);
            else
                return complext<T>(real(x));
        }
        T min;
    };

    template <class T> void clamp_min(hoNDArray<T>* x, typename realType<T>::Type min) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::clamp_min(): Invalid input array");

        hoNDA_clamp_min<T> functor(min);
        std::transform(x->begin(), x->end(), x->begin(), functor);
    }

    template <typename T>
    struct hoNDA_clamp_max //: public thrust::unary_function<T,T>
    {
        hoNDA_clamp_max(T _max) : max(_max) {}
        T operator()(const T& x) const {
            if (x > max)
                return max;
            else
                return x;
        }
        T max;
    };

    template <typename T>
    struct hoNDA_clamp_max<std::complex<T>> //: public thrust::unary_function< std::complex<T>, std::complex<T> >
    {
        hoNDA_clamp_max(T _max) : max(_max) {}
        std::complex<T> operator()(const std::complex<T>& x) const {
            if (real(x) > max)
                return std::complex<T>(max);
            else
                return std::complex<T>(real(x));
        }
        T max;
    };

    template <typename T>
    struct hoNDA_clamp_max<complext<T>> //: public thrust::unary_function< complext<T>, complext<T> >
    {
        hoNDA_clamp_max(T _max) : max(_max) {}
        complext<T> operator()(const complext<T>& x) const {
            if (real(x) > max)
                return complext<T>(max);
            else
                return complext<T>(real(x));
        }
        T max;
    };

    template <class T> void clamp_max(hoNDArray<T>* x, typename realType<T>::Type max) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::clamp_max(): Invalid input array");

        hoNDA_clamp_max<T> functor(max);
        std::transform(x->begin(), x->end(), x->begin(), functor);
    }

    template <class T> void normalize(hoNDArray<T>* x, typename realType<T>::Type val) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::normalize(): Invalid input array");

        size_t max_idx                   = amax(x);
        T max_val_before                 = x->get_data_ptr()[max_idx];
        typename realType<T>::Type scale = val / abs(max_val_before);
        *x *= scale;
    }

    template <class T> void shrink1(hoNDArray<T>* x, typename realType<T>::Type gamma, hoNDArray<T>* out) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::shrink1(): Invalid input array");

        T* outPtr = (out == 0x0) ? x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (long long i = 0; i < (long long)x->get_number_of_elements(); i++) {
            T prev                             = x->get_data_ptr()[i];
            typename realType<T>::Type absPrev = abs(prev);
            T sgnPrev                          = (absPrev <= typename realType<T>::Type(0)) ? T(0) : prev / absPrev;
            outPtr[i]                          = sgnPrev * std::max(absPrev - gamma, typename realType<T>::Type(0));
        }
    }

    template <class T>
    void pshrink(hoNDArray<T>* x, typename realType<T>::Type gamma, typename realType<T>::Type p, hoNDArray<T>* out) {
        if (x == 0x0)
            throw std::runtime_error("Gadgetron::pshrink(): Invalid input array");

        T* outPtr = (out == 0x0) ? x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (long long i = 0; i < (long long)x->get_number_of_elements(); i++) {
            T prev                             = x->get_data_ptr()[i];
            typename realType<T>::Type absPrev = abs(prev);
            T sgnPrev                          = (absPrev <= typename realType<T>::Type(0)) ? T(0) : prev / absPrev;
            outPtr[i] = sgnPrev * std::max(absPrev - gamma * std::pow(absPrev, p - 1), typename realType<T>::Type(0));
        }
    }

    template <class T>
    void shrinkd(hoNDArray<T>* _x, hoNDArray<typename realType<T>::Type>* _s, typename realType<T>::Type gamma,
        hoNDArray<T>* out) {
        if (_x == 0x0 || _s == 0)
            throw std::runtime_error("Gadgetron::shrinkd(): Invalid input array");

        T* outPtr = (out == 0x0) ? _x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (long long i = 0; i < (long long)_x->get_number_of_elements(); i++) {
            T x                          = _x->get_data_ptr()[i];
            typename realType<T>::Type s = _s->get_data_ptr()[i];
            if (s > gamma)
                outPtr[i] = x / s * (s - gamma);
            else
                outPtr[i] = 0;
        }
    }

    template <class T>
    void pshrinkd(hoNDArray<T>* _x, hoNDArray<typename realType<T>::Type>* _s, typename realType<T>::Type gamma,
        typename realType<T>::Type p, hoNDArray<T>* out) {
        if (_x == 0x0)
            throw std::runtime_error("Gadgetron::pshrinkd(): Invalid input array");

        T* outPtr = (out == 0x0) ? _x->get_data_ptr() : out->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
        for (long long i = 0; i < (long long)_x->get_number_of_elements(); i++) {
            T x                          = _x->get_data_ptr()[i];
            typename realType<T>::Type s = _s->get_data_ptr()[i];
            outPtr[i] = x / s * std::max(s - gamma * std::pow(s, p - 1), typename realType<T>::Type(0));
        }
    }

    // --------------------------------------------------------------------------------
    // template<class T> void axpy(T a, const hoNDArray<T> *x, hoNDArray<T> *y)
    //{
    //    BLAS::axpy(x->get_number_of_elements(), a, x->get_data_ptr(), 1, y->get_data_ptr(), 1);
    //}

    // template <typename T> void axpy(T a, const hoNDArray<T>& x, hoNDArray<T>& y)
    //{
    //    axpy(a, &x, &y);
    //}

    // template <typename R, typename T> void scal(R a, hoNDArray<T>& x)
    //{
    //    BLAS::scal(x.get_number_of_elements(), a, x.get_data_ptr(), 1);
    //}

    // --------------------------------------------------------------------------------

    template <typename T>
    void conv2(size_t RO, size_t E1, size_t num, const T* x, size_t kRO, size_t kE1, const T* y, T* z) {
        try {
            long long halfKRO = (long long)(kRO / 2);
            long long halfKE1 = (long long)(kE1 / 2);

            hoNDArray<T> flipY(2 * halfKRO + 1, 2 * halfKE1 + 1);
            T* pKer = flipY.begin();

            long long n;
            long long ro, e1;

            // flip the kernel
            for (e1 = 0; e1 < (long long)kE1; e1++) {
                long long flip_e1 = 2 * halfKE1 - e1;

                for (ro = 0; ro < (long long)kRO; ro++) {
                    long long flip_ro = 2 * halfKRO - ro;

                    flipY(flip_ro, flip_e1) = y[ro + e1 * kRO];
                }
            }

// perform the convolution
#pragma omp parallel for default(none) private(n, ro, e1) shared(num, x, RO, E1, z, halfKRO, halfKE1, pKer)
            for (n = 0; n < (long long)num; n++) {
                const T* pX = x + n * RO * E1;
                T* pZ       = z + n * RO * E1;

                long long kro, ke1, dro, de1;

                for (e1 = 0; e1 < (long long)E1; e1++) {
                    for (ro = 0; ro < (long long)RO; ro++) {
                        pZ[ro + e1 * RO] = 0;
                        for (ke1 = -halfKE1; ke1 <= halfKE1; ke1++) {
                            de1 = ke1 + e1;
                            if (de1 < 0) {
                                de1 += E1;
                            } else if (de1 >= (long long)E1) {
                                de1 -= E1;
                            }

                            for (kro = -halfKRO; kro <= halfKRO; kro++) {
                                dro = kro + ro;
                                if (dro < 0) {
                                    dro += RO;
                                } else if (dro >= (long long)RO) {
                                    dro -= RO;
                                }

                                pZ[ro + e1 * RO]
                                    += pKer[kro + halfKRO + (ke1 + halfKE1) * (2 * halfKRO + 1)] * pX[dro + de1 * RO];
                            }
                        }
                    }
                }
            }
        } catch (...) {
            GADGET_THROW("Errors happened in conv2(size_t RO, size_t E1, size_t num, const T* x, size_t kRO, size_t "
                         "kE1, const T* y, T* z) ... ");
        }
    }

    template <typename T> void conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) {
        try {
            if (!z.dimensions_equal(&x)) {
                z = x;
            }

            long long RO  = (long long)x.get_size(0);
            long long E1  = (long long)x.get_size(1);
            long long num = ((long long)x.get_number_of_elements()) / (RO * E1);

            long long kRO = (long long)y.get_size(0);
            long long kE1 = (long long)y.get_size(1);

            conv2(RO, E1, num, x.begin(), kRO, kE1, y.begin(), z.begin());
        } catch (...) {
            GADGET_THROW(
                "Errors happened in conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) ... ");
        }
    }

    template  void conv2(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z);
    template  void conv2(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    template  void conv2(const hoNDArray<std::complex<float>>& x,
        const hoNDArray<std::complex<float>>& y, hoNDArray<std::complex<float>>& z);
    template  void conv2(const hoNDArray<std::complex<double>>& x,
        const hoNDArray<std::complex<double>>& y, hoNDArray<std::complex<double>>& z);

    template  void conv2(
        const hoNDArray<complext<float>>& x, const hoNDArray<complext<float>>& y, hoNDArray<complext<float>>& z);
    template  void conv2(
        const hoNDArray<complext<double>>& x, const hoNDArray<complext<double>>& y, hoNDArray<complext<double>>& z);
    // --------------------------------------------------------------------------------

    template <typename T>
    void conv3(
        size_t RO, size_t E1, size_t E2, size_t num, const T* x, size_t kRO, size_t kE1, size_t kE2, const T* y, T* z) {
        try {
            long long halfKRO = (long long)(kRO / 2);
            long long halfKE1 = (long long)(kE1 / 2);
            long long halfKE2 = (long long)(kE2 / 2);

            hoNDArray<T> flipY(2 * halfKRO + 1, 2 * halfKE1 + 1, 2 * halfKE2 + 1);
            T* pKer = flipY.begin();

            long long n, e2;
            long long ro, e1;

            // flip the kernel
            for (e2 = 0; e2 < (long long)kE2; e2++) {
                long long flip_e2 = 2 * halfKE2 - e2;

                for (e1 = 0; e1 < (long long)kE1; e1++) {
                    long long flip_e1 = 2 * halfKE1 - e1;

                    for (ro = 0; ro < (long long)kRO; ro++) {
                        long long flip_ro = 2 * halfKRO - ro;

                        flipY(flip_ro, flip_e1, flip_e2) = y[ro + e1 * kRO + e2 * kRO * kE1];
                    }
                }
            }

// perform the convolution
#pragma omp parallel for default(none) private(n)                                                                      \
    shared(num, x, RO, E1, E2, z, halfKRO, halfKE1, halfKE2, pKer, e1, ro) if (num > 8)
            for (n = 0; n < (long long)num; n++) {
                const T* pX = x + n * RO * E1 * E2;
                T* pZ       = z + n * RO * E1 * E2;

                long long kro, ke1, ke2, dro, de1, de2;

#pragma omp parallel for default(none) private(ro, e1, e2, kro, ke1, ke2, dro, de1, de2)                               \
    shared(pX, RO, E1, E2, pZ, halfKRO, halfKE1, halfKE2, pKer)
                for (e2 = 0; e2 < (long long)E2; e2++) {
                    for (e1 = 0; e1 < (long long)E1; e1++) {
                        for (ro = 0; ro < (long long)RO; ro++) {
                            pZ[ro + e1 * RO + e2 * RO * E1] = 0;
                            for (ke2 = -halfKE2; ke2 <= halfKE2; ke2++) {
                                de2 = ke2 + e2;
                                if (de2 < 0) {
                                    de2 += E2;
                                } else if (de2 >= (long long)E2) {
                                    de2 -= E2;
                                }

                                for (ke1 = -halfKE1; ke1 <= halfKE1; ke1++) {
                                    de1 = ke1 + e1;
                                    if (de1 < 0) {
                                        de1 += E1;
                                    } else if (de1 >= (long long)E1) {
                                        de1 -= E1;
                                    }

                                    for (kro = -halfKRO; kro <= halfKRO; kro++) {
                                        dro = kro + ro;
                                        if (dro < 0) {
                                            dro += RO;
                                        } else if (dro >= (long long)RO) {
                                            dro -= RO;
                                        }

                                        pZ[ro + e1 * RO + e2 * RO * E1]
                                            += pKer[kro + halfKRO + (ke1 + halfKE1) * (2 * halfKRO + 1)
                                                    + (ke2 + halfKE2) * (2 * halfKRO + 1) * (2 * halfKE1 + 1)]
                                               * pX[dro + de1 * RO + de2 * RO * E1];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } catch (...) {
            GADGET_THROW("Errors happened in conv3(size_t RO, size_t E1, size_t E2, size_t num, const T* x, size_t "
                         "kRO, size_t kE1, size_t kE2, const T* y, T* z) ... ");
        }
    }

    template <typename T> void conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) {
        try {
            if (!z.dimensions_equal(&x)) {
                z = x;
            }

            long long RO  = (long long)x.get_size(0);
            long long E1  = (long long)x.get_size(1);
            long long E2  = (long long)x.get_size(2);
            long long num = ((long long)x.get_number_of_elements()) / (RO * E1 * E2);

            long long kRO = (long long)y.get_size(0);
            long long kE1 = (long long)y.get_size(1);
            long long kE2 = (long long)y.get_size(2);

            conv3(RO, E1, E2, num, x.begin(), kRO, kE1, kE2, y.begin(), z.begin());
        } catch (...) {
            GADGET_THROW(
                "Errors happened in conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) ... ");
        }
    }

    template  void conv3(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z);
    template  void conv3(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    template  void conv3(const hoNDArray<std::complex<float>>& x,
        const hoNDArray<std::complex<float>>& y, hoNDArray<std::complex<float>>& z);
    template  void conv3(const hoNDArray<std::complex<double>>& x,
        const hoNDArray<std::complex<double>>& y, hoNDArray<std::complex<double>>& z);

    template  void conv3(
        const hoNDArray<complext<float>>& x, const hoNDArray<complext<float>>& y, hoNDArray<complext<float>>& z);
    template  void conv3(
        const hoNDArray<complext<double>>& x, const hoNDArray<complext<double>>& y, hoNDArray<complext<double>>& z);
    // --------------------------------------------------------------------------------

    template <typename T> void sum_over_dimension(const hoNDArray<T>& x, hoNDArray<T>& r, size_t dim) {
        try {
            size_t D = x.get_number_of_dimensions();
            if (dim >= D) {
                r = x;
                return;
            }

            auto dimX = x.dimensions();

            auto dimR      = dimX;
            dimR[dim] = 1;

            if (!r.dimensions_equal(&dimR)) {
                r.create(dimR);
            }

            if (dim == 0) {
                size_t X   = x.get_size(0);
                size_t num = x.get_number_of_elements() / X;

                const T* pX = x.begin();
                T* pR       = r.begin();

                long long n;

                for (n = 0; n < (long long)num; n++) {
                    T xsum = pX[n * X];
                    for (size_t ro = 1; ro < X; ro++) {
                        xsum += pX[n * X + ro];
                    }

                    pR[n] = xsum;
                }
            } else {
                size_t strideX = x.get_size(0);
                for (size_t d = 1; d <= dim; d++) {
                    strideX *= x.get_size(d);
                }

                size_t strideR = strideX / x.get_size(dim);
                size_t num     = x.get_number_of_elements() / strideX;
                size_t nDim    = x.get_size(dim);

                const T* pX = x.begin();
                T* pR       = r.begin();

                if (nDim == 1) {
                    memcpy(pR, pX, x.get_number_of_bytes());
                    return;
                }

                long long n;

                for (n = 0; n < (long long)num; n++) {
                    const T* pX_curr = pX + n * strideX;
                    T* pR_curr       = pR + n * strideR;

                    memcpy(pR_curr, pX_curr, sizeof(T) * strideR);

                    size_t p, c;
                    for (p = 1; p < nDim; p++) {
                        for (c = 0; c < strideR; c++) {
                            pR_curr[c] += pX_curr[p * strideR + c];
                        }
                    }
                }
            }
        } catch (...) {
            GADGET_THROW(
                "Errors happened in sum_over_dimension(const hoNDArray<T>& x, hoNDArray<T>& y, size_t dim) ... ");
        }
    }

    template  void sum_over_dimension(const hoNDArray<float>& x, hoNDArray<float>& y, size_t dim);
    template  void sum_over_dimension(const hoNDArray<double>& x, hoNDArray<double>& y, size_t dim);
    template  void sum_over_dimension(
        const hoNDArray<std::complex<float>>& x, hoNDArray<std::complex<float>>& y, size_t dim);
    template  void sum_over_dimension(
        const hoNDArray<std::complex<double>>& x, hoNDArray<std::complex<double>>& y, size_t dim);

    template  void sum_over_dimension(
        const hoNDArray<complext<float>>& x, hoNDArray<complext<float>>& y, size_t dim);
    template  void sum_over_dimension(
        const hoNDArray<complext<double>>& x, hoNDArray<complext<double>>& y, size_t dim);
    // --------------------------------------------------------------------------------


    hoNDArray<bool>& operator&=(hoNDArray<bool>& x, const hoNDArray<bool>& y) {
        if (compatible_dimensions<bool, bool>(x, y)) {
            const size_t elementsX = x.get_number_of_elements();
            const size_t elementsY = y.get_number_of_elements();
            bool* x_ptr            = x.get_data_ptr();
            const bool* y_ptr      = y.get_data_ptr();
            for (size_t i = 0; i < elementsX; i++)
                x_ptr[i] &= y_ptr[i % elementsY];
            return x;
        } else {
            throw std::runtime_error("&= incompatible dimensions.");
        }
    }
    hoNDArray<bool>& operator|=(hoNDArray<bool>& x, const hoNDArray<bool>& y) {
        if (compatible_dimensions<bool, bool>(x, y)) {
            const size_t elementsX = x.get_number_of_elements();
            const size_t elementsY = y.get_number_of_elements();
            bool* x_ptr            = x.get_data_ptr();
            const bool* y_ptr      = y.get_data_ptr();
            for (size_t i = 0; i < elementsX; i++)
                x_ptr[i] |= y_ptr[i % elementsY];
            return x;
        } else {
            throw std::runtime_error("|= incompatible dimensions.");
        }
    }
    // --------------------------------------------------------------------------------

    //
    // Instantiation
    //

    template  boost::shared_ptr<hoNDArray<float>> abs<float>(hoNDArray<float>*);
    template  void abs_inplace<float>(hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float>> abs_square<float>(hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float>> sqrt<float>(hoNDArray<float>*);
    template  void sqrt_inplace<float>(hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float>> square<float>(hoNDArray<float>*);
    template  void square_inplace<float>(hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float>> reciprocal<float>(hoNDArray<float>*);
    template  void reciprocal_inplace<float>(hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float>> reciprocal_sqrt<float>(hoNDArray<float>*);
    template  void reciprocal_sqrt_inplace<float>(hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float>> sgn<float>(hoNDArray<float>*);
    template  void sgn_inplace<float>(hoNDArray<float>*);
    template  void clamp<float>(hoNDArray<float>*, float, float);
    template  void clamp<float>(hoNDArray<float>*, float, float, float, float);
    template  void clamp_min<float>(hoNDArray<float>*, float);
    template  void clamp_max<float>(hoNDArray<float>*, float);
    template  void normalize<float>(hoNDArray<float>*, float);
    template  void shrink1<float>(hoNDArray<float>*, float, hoNDArray<float>*);
    template  void pshrink<float>(hoNDArray<float>*, float, float, hoNDArray<float>*);
    template  void shrinkd<float>(hoNDArray<float>*, hoNDArray<float>*, float, hoNDArray<float>*);
    template  void pshrinkd<float>(
        hoNDArray<float>*, hoNDArray<float>*, float, float, hoNDArray<float>*);

    template  boost::shared_ptr<hoNDArray<double>> abs<double>(hoNDArray<double>*);
    template  void abs_inplace<double>(hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double>> abs_square<double>(hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double>> sqrt<double>(hoNDArray<double>*);
    template  void sqrt_inplace<double>(hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double>> square<double>(hoNDArray<double>*);
    template  void square_inplace<double>(hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double>> reciprocal<double>(hoNDArray<double>*);
    template  void reciprocal_inplace<double>(hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double>> reciprocal_sqrt<double>(hoNDArray<double>*);
    template  void reciprocal_sqrt_inplace<double>(hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double>> sgn<double>(hoNDArray<double>*);
    template  void sgn_inplace<double>(hoNDArray<double>*);
    template  void clamp<double>(hoNDArray<double>*, double, double);
    template  void clamp_min<double>(hoNDArray<double>*, double);
    template  void clamp_max<double>(hoNDArray<double>*, double);
    template  void normalize<double>(hoNDArray<double>*, double);
    template  void shrink1<double>(hoNDArray<double>*, double, hoNDArray<double>*);
    template  void pshrink<double>(hoNDArray<double>*, double, double, hoNDArray<double>*);
    template  void shrinkd<double>(hoNDArray<double>*, hoNDArray<double>*, double, hoNDArray<double>*);
    template  void pshrinkd<double>(
        hoNDArray<double>*, hoNDArray<double>*, double, double, hoNDArray<double>*);

    template  boost::shared_ptr<hoNDArray<float>> abs<std::complex<float>>(
        hoNDArray<std::complex<float>>*);
    template  boost::shared_ptr<hoNDArray<float>> abs_square<std::complex<float>>(
        hoNDArray<std::complex<float>>*);
    template  boost::shared_ptr<hoNDArray<std::complex<float>>> sqrt<std::complex<float>>(
        hoNDArray<std::complex<float>>*);
    template  void sqrt_inplace<std::complex<float>>(hoNDArray<std::complex<float>>*);
    template  void sgn_inplace<std::complex<float>>(hoNDArray<std::complex<float>>*);
template  boost::shared_ptr<hoNDArray<std::complex<float>>> square<std::complex<float>>(
        hoNDArray<std::complex<float>>*);
    template  void square_inplace<std::complex<float>>(hoNDArray<std::complex<float>>*);
    template  boost::shared_ptr<hoNDArray<std::complex<float>>> reciprocal<std::complex<float>>(
        hoNDArray<std::complex<float>>*);
    template  void reciprocal_inplace<std::complex<float>>(hoNDArray<std::complex<float>>*);
    template  boost::shared_ptr<hoNDArray<std::complex<float>>> reciprocal_sqrt<std::complex<float>>(
        hoNDArray<std::complex<float>>*);
    template  void reciprocal_sqrt_inplace<std::complex<float>>(hoNDArray<std::complex<float>>*);

    template  void clamp<std::complex<float>>(hoNDArray<std::complex<float>>*, float, float);
    template  void clamp_min<std::complex<float>>(hoNDArray<std::complex<float>>*, float);
    template  void clamp_max<std::complex<float>>(hoNDArray<std::complex<float>>*, float);
    template  void normalize<std::complex<float>>(hoNDArray<std::complex<float>>*, float);
    template  void shrink1<std::complex<float>>(
        hoNDArray<std::complex<float>>*, float, hoNDArray<std::complex<float>>*);
    template  void pshrink<std::complex<float>>(
        hoNDArray<std::complex<float>>*, float, float, hoNDArray<std::complex<float>>*);
    template  void shrinkd<std::complex<float>>(
        hoNDArray<std::complex<float>>*, hoNDArray<float>*, float, hoNDArray<std::complex<float>>*);
    template  void pshrinkd<std::complex<float>>(
        hoNDArray<std::complex<float>>*, hoNDArray<float>*, float, float, hoNDArray<std::complex<float>>*);

    template  boost::shared_ptr<hoNDArray<double>> abs<std::complex<double>>(
        hoNDArray<std::complex<double>>*);
    template  boost::shared_ptr<hoNDArray<double>> abs_square<std::complex<double>>(
        hoNDArray<std::complex<double>>*);
    template  boost::shared_ptr<hoNDArray<std::complex<double>>> sqrt<std::complex<double>>(
        hoNDArray<std::complex<double>>*);
    template  void sqrt_inplace<std::complex<double>>(hoNDArray<std::complex<double>>*);
    template  void sgn_inplace<std::complex<double>>(hoNDArray<std::complex<double>>*);
template  boost::shared_ptr<hoNDArray<std::complex<double>>> square<std::complex<double>>(
        hoNDArray<std::complex<double>>*);
    template  void square_inplace<std::complex<double>>(hoNDArray<std::complex<double>>*);
    template  boost::shared_ptr<hoNDArray<std::complex<double>>> reciprocal<std::complex<double>>(
        hoNDArray<std::complex<double>>*);
    template  void reciprocal_inplace<std::complex<double>>(hoNDArray<std::complex<double>>*);
    template  boost::shared_ptr<hoNDArray<std::complex<double>>> reciprocal_sqrt<std::complex<double>>(
        hoNDArray<std::complex<double>>*);
    template  void reciprocal_sqrt_inplace<std::complex<double>>(hoNDArray<std::complex<double>>*);
    template  void clamp<std::complex<double>>(hoNDArray<std::complex<double>>*, double, double);
    template  void clamp_min<std::complex<double>>(hoNDArray<std::complex<double>>*, double);
    template  void clamp_max<std::complex<double>>(hoNDArray<std::complex<double>>*, double);
    template  void normalize<std::complex<double>>(hoNDArray<std::complex<double>>*, double);
    template  void shrink1<std::complex<double>>(
        hoNDArray<std::complex<double>>*, double, hoNDArray<std::complex<double>>*);
    template  void pshrink<std::complex<double>>(
        hoNDArray<std::complex<double>>*, double, double, hoNDArray<std::complex<double>>*);
    template  void shrinkd<std::complex<double>>(
        hoNDArray<std::complex<double>>*, hoNDArray<double>*, double, hoNDArray<std::complex<double>>*);
    template  void pshrinkd<std::complex<double>>(
        hoNDArray<std::complex<double>>*, hoNDArray<double>*, double, double, hoNDArray<std::complex<double>>*);

    template  boost::shared_ptr<hoNDArray<float>> abs<complext<float>>(hoNDArray<complext<float>>*);
    template  boost::shared_ptr<hoNDArray<float>> abs_square<complext<float>>(
        hoNDArray<complext<float>>*);
    template  boost::shared_ptr<hoNDArray<complext<float>>> sqrt<complext<float>>(
        hoNDArray<complext<float>>*);
    template  void sqrt_inplace<complext<float>>(hoNDArray<complext<float>>*);
    template  boost::shared_ptr<hoNDArray<complext<float>>> square<complext<float>>(
        hoNDArray<complext<float>>*);
    template  void square_inplace<complext<float>>(hoNDArray<complext<float>>*);
    template  boost::shared_ptr<hoNDArray<complext<float>>> reciprocal<complext<float>>(
        hoNDArray<complext<float>>*);
    template  void reciprocal_inplace<complext<float>>(hoNDArray<complext<float>>*);
    template  boost::shared_ptr<hoNDArray<complext<float>>> reciprocal_sqrt<complext<float>>(
        hoNDArray<complext<float>>*);
    template  void reciprocal_sqrt_inplace<complext<float>>(hoNDArray<complext<float>>*);
    template  void clamp<complext<float>>(hoNDArray<complext<float>>*, float, float);
    template  void clamp_min<complext<float>>(hoNDArray<complext<float>>*, float);
    template  void clamp_max<complext<float>>(hoNDArray<complext<float>>*, float);
    template  void normalize<complext<float>>(hoNDArray<complext<float>>*, float);
    template  void shrink1<complext<float>>(
        hoNDArray<complext<float>>*, float, hoNDArray<complext<float>>*);
    template  void pshrink<complext<float>>(
        hoNDArray<complext<float>>*, float, float, hoNDArray<complext<float>>*);
    template  void shrinkd<complext<float>>(
        hoNDArray<complext<float>>*, hoNDArray<float>*, float, hoNDArray<complext<float>>*);
    template  void pshrinkd<complext<float>>(
        hoNDArray<complext<float>>*, hoNDArray<float>*, float, float, hoNDArray<complext<float>>*);

    template  boost::shared_ptr<hoNDArray<double>> abs<complext<double>>(hoNDArray<complext<double>>*);
    template  boost::shared_ptr<hoNDArray<double>> abs_square<complext<double>>(
        hoNDArray<complext<double>>*);
    template  boost::shared_ptr<hoNDArray<complext<double>>> sqrt<complext<double>>(
        hoNDArray<complext<double>>*);
    template  void sqrt_inplace<complext<double>>(hoNDArray<complext<double>>*);
    template  boost::shared_ptr<hoNDArray<complext<double>>> square<complext<double>>(
        hoNDArray<complext<double>>*);
    template  void square_inplace<complext<double>>(hoNDArray<complext<double>>*);
    template  boost::shared_ptr<hoNDArray<complext<double>>> reciprocal<complext<double>>(
        hoNDArray<complext<double>>*);
    template  void reciprocal_inplace<complext<double>>(hoNDArray<complext<double>>*);
    template  boost::shared_ptr<hoNDArray<complext<double>>> reciprocal_sqrt<complext<double>>(
        hoNDArray<complext<double>>*);
    template  void reciprocal_sqrt_inplace<complext<double>>(hoNDArray<complext<double>>*);
    template  void clamp<complext<double>>(hoNDArray<complext<double>>*, double, double);
    template  void clamp_min<complext<double>>(hoNDArray<complext<double>>*, double);
    template  void clamp_max<complext<double>>(hoNDArray<complext<double>>*, double);
    template  void normalize<complext<double>>(hoNDArray<complext<double>>*, double);
    template  void shrink1<complext<double>>(
        hoNDArray<complext<double>>*, double, hoNDArray<complext<double>>*);
    template  void pshrink<complext<double>>(
        hoNDArray<complext<double>>*, double, double, hoNDArray<complext<double>>*);
    template  void shrinkd<complext<double>>(
        hoNDArray<complext<double>>*, hoNDArray<double>*, double, hoNDArray<complext<double>>*);
    template  void pshrinkd<complext<double>>(
        hoNDArray<complext<double>>*, hoNDArray<double>*, double, double, hoNDArray<complext<double>>*);

    template  boost::shared_ptr<hoNDArray<std::complex<float>>> real_to_complex<std::complex<float>>(
        const hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<std::complex<float>>>
    real_imag_to_complex<std::complex<float>>(hoNDArray<float>*, hoNDArray<float>*);

    template  boost::shared_ptr<hoNDArray<float_complext>> real_to_complex<float_complext>(
        const hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float_complext>> real_imag_to_complex<float_complext>(
        hoNDArray<float>*, hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float>> real<float>(const hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float>> real<std::complex<float>>(
        const hoNDArray<std::complex<float>>*);
    template  boost::shared_ptr<hoNDArray<float>> real<float_complext>(
        const hoNDArray<float_complext>*);
    template  boost::shared_ptr<hoNDArray<float>> imag<float>(const hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<float>> imag<std::complex<float>>(
        const hoNDArray<std::complex<float>>*);
    template  boost::shared_ptr<hoNDArray<float>> imag<float_complext>(
        const hoNDArray<float_complext>*);
    template  boost::shared_ptr<hoNDArray<float>> conj<float>(const hoNDArray<float>*);
    template  boost::shared_ptr<hoNDArray<std::complex<float>>> conj<std::complex<float>>(
        const hoNDArray<std::complex<float>>*);
    template  boost::shared_ptr<hoNDArray<float_complext>> conj<float_complext>(
        const hoNDArray<float_complext>*);

    template  boost::shared_ptr<hoNDArray<std::complex<double>>> real_to_complex<std::complex<double>>(
        const hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<std::complex<double>>>
    real_imag_to_complex<std::complex<double>>(hoNDArray<double>*, hoNDArray<double>*);

    template  boost::shared_ptr<hoNDArray<double_complext>> real_to_complex<double_complext>(
        const hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double_complext>> real_imag_to_complex<double_complext>(
        hoNDArray<double>*, hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double>> real<double>(const hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double>> real<std::complex<double>>(
        const hoNDArray<std::complex<double>>*);
    template  boost::shared_ptr<hoNDArray<double>> real<double_complext>(
        const hoNDArray<double_complext>*);
    template  boost::shared_ptr<hoNDArray<double>> imag<std::complex<double>>(
        const hoNDArray<std::complex<double>>*);
    template  boost::shared_ptr<hoNDArray<double>> imag<double>(const hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<double>> imag<double_complext>(
        const hoNDArray<double_complext>*);
    template  boost::shared_ptr<hoNDArray<double>> conj<double>(const hoNDArray<double>*);
    template  boost::shared_ptr<hoNDArray<std::complex<double>>> conj<std::complex<double>>(
        const hoNDArray<std::complex<double>>*);
    template  boost::shared_ptr<hoNDArray<double_complext>> conj<double_complext>(
        const hoNDArray<double_complext>*);



}
