#include "hoNDArray_math_util.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------------------

    template<typename T> void fill( hoNDArray<T>* x, T val)
    {
        size_t N = x->get_number_of_elements();
        T* pX = x->begin();
        Gadgetron::math::fill(N, pX, val);
    }

    template EXPORTCPUCOREMATH void fill( hoNDArray<float>* x, float val);
    template EXPORTCPUCOREMATH void fill( hoNDArray<double>* x, double val);
    template EXPORTCPUCOREMATH void fill( hoNDArray<GT_Complex8>* x, GT_Complex8 val);
    template EXPORTCPUCOREMATH void fill( hoNDArray<GT_Complex16>* x, GT_Complex16 val);

    // --------------------------------------------------------------------------------

    template<typename T> void fill( hoNDArray<T>& x, T val )
    {
        size_t N = x.get_number_of_elements();
        T* pX = x.begin();
        Gadgetron::math::fill(N, pX, val);
    }

    template EXPORTCPUCOREMATH void fill( hoNDArray<float>& x, float val);
    template EXPORTCPUCOREMATH void fill( hoNDArray<double>& x, double val);
    template EXPORTCPUCOREMATH void fill( hoNDArray<GT_Complex8>& x, GT_Complex8 val);
    template EXPORTCPUCOREMATH void fill( hoNDArray<GT_Complex16>& x, GT_Complex16 val);

    // --------------------------------------------------------------------------------

    template<class T> 
    bool real_imag_to_complex(const hoNDArray<typename realType<T>::Type>& real, const hoNDArray<typename realType<T>::Type>& imag, hoNDArray<T>& cplx)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(real.dimensions_equal(&imag));

            if ( !cplx.dimensions_equal(&real) )
            {
                cplx.create(real.get_dimensions());
            }

            T* pRes = cplx.begin();
            const typename realType<T>::Type* pReal = real.begin();
            const typename realType<T>::Type* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for private(n) shared(N, pRes, pReal, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pRes[n] = T(pReal[n], pImag[n]);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in real_imag_to_complex(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool real_imag_to_complex(const hoNDArray<float>& real, const hoNDArray<float>& imag, hoNDArray<GT_Complex8>& cplx);
    template EXPORTCPUCOREMATH bool real_imag_to_complex(const hoNDArray<double>& real, const hoNDArray<double>& imag, hoNDArray<GT_Complex16>& cplx);

    // --------------------------------------------------------------------------------

    template<class T> 
    bool complex_to_real_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real, hoNDArray<typename realType<T>::Type>& imag)
    {
        try
        {
            if ( !real.dimensions_equal(&cplx) )
            {
                real.create(cplx.get_dimensions());
            }

            if ( !imag.dimensions_equal(&cplx) )
            {
                imag.create(cplx.get_dimensions());
            }

            const T* pRes = cplx.begin();
            typename realType<T>::Type* pReal = real.begin();
            typename realType<T>::Type* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pReal, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pReal[n] = pRes[n].real();
                pImag[n] = pRes[n].imag();
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real_imag(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray<GT_Complex8>& cplx, hoNDArray<float>& real, hoNDArray<float>& imag);
    template EXPORTCPUCOREMATH bool complex_to_real_imag(const hoNDArray<GT_Complex16>& cplx, hoNDArray<double>& real, hoNDArray<double>& imag);

    template <> EXPORTCPUCOREMATH
    bool complex_to_real_imag(const hoNDArray<float>& cplx, hoNDArray<float>& real, hoNDArray<float>& imag)
    {
        try
        {
            if ( !real.dimensions_equal(&cplx) )
            {
                real.create(cplx.get_dimensions());
            }

            if ( !imag.dimensions_equal(&cplx) )
            {
                imag.create(cplx.get_dimensions());
            }

            const float* pRes = cplx.begin();
            float* pReal = real.begin();
            float* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pReal, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pReal[n] = pRes[n];
                pImag[n] = 0;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real_imag(...) ... ");
            return false;
        }

        return true;
    }

    template<> EXPORTCPUCOREMATH 
    bool complex_to_real_imag(const hoNDArray<double>& cplx, hoNDArray<double>& real, hoNDArray<double>& imag)
    {
        try
        {
            if ( !real.dimensions_equal(&cplx) )
            {
                real.create(cplx.get_dimensions());
            }

            if ( !imag.dimensions_equal(&cplx) )
            {
                imag.create(cplx.get_dimensions());
            }

            const double* pRes = cplx.begin();
            double* pReal = real.begin();
            double* pImag = imag.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pReal, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pReal[n] = pRes[n];
                pImag[n] = 0;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real_imag(...) ... ");
            return false;
        }

        return true;
    }

    // --------------------------------------------------------------------------------

    template<class T> 
    bool complex_to_real(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& real)
    {
        try
        {
            if ( !real.dimensions_equal(&cplx) )
            {
                real.create(cplx.get_dimensions());
            }

            const T* pRes = cplx.begin();
            typename realType<T>::Type* pReal = real.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pReal)
            for ( n=0; n<(long long)N; n++ )
            {
                pReal[n] = pRes[n].real();
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_real(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_real(const hoNDArray<GT_Complex8>& cplx, hoNDArray<float>& real);
    template EXPORTCPUCOREMATH bool complex_to_real(const hoNDArray<GT_Complex16>& cplx, hoNDArray<double>& real);

    // --------------------------------------------------------------------------------

    template<class T> 
    bool complex_to_imag(const hoNDArray<T>& cplx, hoNDArray<typename realType<T>::Type>& imag)
    {
        try
        {
            if ( !imag.dimensions_equal(&cplx) )
            {
                imag.create(cplx.get_dimensions());
            }

            const T* pRes = cplx.begin();
            typename realType<T>::Type* pImag = imag.begin();

            size_t N = imag.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pImag)
            for ( n=0; n<(long long)N; n++ )
            {
                pImag[n] = pRes[n].imag();
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in complex_to_imag(...) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool complex_to_imag(const hoNDArray<GT_Complex8>& cplx, hoNDArray<float>& imag);
    template EXPORTCPUCOREMATH bool complex_to_imag(const hoNDArray<GT_Complex16>& cplx, hoNDArray<double>& imag);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool absolute(const hoNDArray< std::complex<T> >& x, hoNDArray< std::complex<T> >& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        hoNDArray<T> rTmp;
        rTmp.create(x.get_dimensions());

        Gadgetron::absolute(x, rTmp);

        r.copyFrom(rTmp);

        return true;
    }

    template EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool add(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        Gadgetron::math::add(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool add(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool add(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool add(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool add(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool subtract(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        Gadgetron::math::subtract(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool subtract(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool subtract(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool subtract(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool subtract(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool multiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        Gadgetron::math::multiply(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool multiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool multiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool multiply(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool multiply(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool divide(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        Gadgetron::math::divide(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool divide(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool divide(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool divide(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool divide(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool sqrt(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        Gadgetron::math::sqrt(x.get_number_of_elements(), x.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool sqrt(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool sqrt(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool sqrt(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool sqrt(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool minAbsolute(const hoNDArray<T>& x, T& r, size_t& ind)
    {
        size_t N = x.get_number_of_elements();
        if ( N == 0 ) return true;
        Gadgetron::math::minAbsolute(N, x.begin(), r, ind);

        return true;
    }

    template EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<float>& x, float& r, size_t& ind);
    template EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    template EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind);
    template EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool maxAbsolute(const hoNDArray<T>& x, T& r, size_t& ind)
    {
        size_t N = x.get_number_of_elements();
        if ( N == 0 ) return true;
        Gadgetron::math::maxAbsolute(N, x.begin(), r, ind);

        return true;
    }

    template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<float>& x, float& r, size_t& ind);
    template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<double>& x, double& r, size_t& ind);
    template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind);
    template EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool multiplyConj(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        Gadgetron::math::multiplyConj(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool conjugate(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        Gadgetron::math::conjugate(x.get_number_of_elements(), x.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool conjugate(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool conjugate(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool addEpsilon(hoNDArray<T>& x)
    {
        Gadgetron::math::addEpsilon(x.get_number_of_elements(), x.begin());
        return true;
    }

    template EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<float>& x);
    template EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<double>& x);
    template EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<GT_Complex8>& x);
    template EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<GT_Complex16>& x);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool norm2(const hoNDArray<T>& x, typename realType<T>::Type& r)
    {
        Gadgetron::math::norm2(x.get_number_of_elements(), x.begin(), r);
        return true;
    }

    template EXPORTCPUCOREMATH bool norm2(const hoNDArray<float>& x, float& r);
    template EXPORTCPUCOREMATH bool norm2(const hoNDArray<double>& x, double& r);
    template EXPORTCPUCOREMATH bool norm2(const hoNDArray<GT_Complex8>& x, float& r);
    template EXPORTCPUCOREMATH bool norm2(const hoNDArray<GT_Complex16>& x, double& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool norm1(const hoNDArray<T>& x, typename realType<T>::Type& r)
    {
        Gadgetron::math::norm1(x.get_number_of_elements(), x.begin(), r);
        return true;
    }

    template EXPORTCPUCOREMATH bool norm1(const hoNDArray<float>& x, float& r);
    template EXPORTCPUCOREMATH bool norm1(const hoNDArray<double>& x, double& r);
    template EXPORTCPUCOREMATH bool norm1(const hoNDArray<GT_Complex8>& x, float& r);
    template EXPORTCPUCOREMATH bool norm1(const hoNDArray<GT_Complex16>& x, double& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool dotc(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        Gadgetron::math::dotc(x.get_number_of_elements(), x.begin(), y.begin(), r);
        return true;
    }

    template EXPORTCPUCOREMATH bool dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, GT_Complex8& r);
    template EXPORTCPUCOREMATH bool dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, GT_Complex16& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool absolute(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        Gadgetron::math::absolute(x.get_number_of_elements(), x.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool absolute(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool absolute(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool argument(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        Gadgetron::math::argument(x.get_number_of_elements(), x.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool argument(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool argument(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool inv(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        if ( !r.dimensions_equal(&x) )
        {
            r = x;
        }

        Gadgetron::math::inv(x.get_number_of_elements(), x.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool inv(const hoNDArray<float>& x, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool inv(const hoNDArray<double>& x, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool inv(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool inv(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template<typename T> 
    bool conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            long long RO = (long long) x.get_size(0);
            long long E1 = (long long) x.get_size(1);
            long long num = ((long long) x.get_number_of_elements()) / (RO*E1);

            long long kRO = (long long) y.get_size(0);
            long long kE1 = (long long) y.get_size(1);

            Gadgetron::math::conv2(RO, E1, num, x.begin(), kRO, kE1, y.begin(), z.begin());
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool conv2(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z);
    template EXPORTCPUCOREMATH bool conv2(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    template EXPORTCPUCOREMATH bool conv2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& z);
    template EXPORTCPUCOREMATH bool conv2(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& z);

    // --------------------------------------------------------------------------------

    template<typename T> 
    bool conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            long long RO = (long long) x.get_size(0);
            long long E1 = (long long) x.get_size(1);
            long long E2 = (long long) x.get_size(2);
            long long num = ((long long)x.get_number_of_elements()) / (RO*E1*E2);

            long long kRO = (long long) y.get_size(0);
            long long kE1 = (long long) y.get_size(1);
            long long kE2 = (long long) y.get_size(2);

            Gadgetron::math::conv3(RO, E1, E2, num, x.begin(), kRO, kE1, kE2, y.begin(), z.begin());
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) ... ");
            return false;
        }

        return true;
    }

    template EXPORTCPUCOREMATH bool conv3(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& z);
    template EXPORTCPUCOREMATH bool conv3(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& z);
    template EXPORTCPUCOREMATH bool conv3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& z);
    template EXPORTCPUCOREMATH bool conv3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& z);

    // ----------------------------------------------------

    template <typename T> 
    T dotc(const hoNDArray<T>& x, const hoNDArray<T>& y)
    {
        T r = 0;

        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        dotc(x, y, r);

        return r;
    }

    template EXPORTCPUCOREMATH GT_Complex8 dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y);
    template EXPORTCPUCOREMATH GT_Complex16 dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y);

    // --------------------------------------------------------------------------------

    template <typename T> 
    T dotu(const hoNDArray<T>& x, const hoNDArray<T>& y)
    {
        T r = 0;

        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        Gadgetron::math::dotu(x.get_number_of_elements(), x.begin(), y.begin(), r);

        return r;
    }

    template EXPORTCPUCOREMATH float dotu(const hoNDArray<float>& x, const hoNDArray<float>& y);
    template EXPORTCPUCOREMATH double dotu(const hoNDArray<double>& x, const hoNDArray<double>& y);
    template EXPORTCPUCOREMATH GT_Complex8 dotu(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y);
    template EXPORTCPUCOREMATH GT_Complex16 dotu(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool axpy(T a, const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        GADGET_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

        if ( r.get_number_of_elements() != x.get_number_of_elements() )
        {
            r = y;
        }
        else
        {
            if ( &r != &y )
            {
                memcpy(r.begin(), y.begin(), r.get_number_of_bytes());
            }
        }

        Gadgetron::math::axpy(a, x.get_number_of_elements(), x.begin(), y.begin(), r.begin());

        return true;
    }

    template EXPORTCPUCOREMATH bool axpy(float a, const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r);
    template EXPORTCPUCOREMATH bool axpy(double a, const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r);
    template EXPORTCPUCOREMATH bool axpy(GT_Complex8 a, const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r);
    template EXPORTCPUCOREMATH bool axpy(GT_Complex16 a, const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool scal(T a, hoNDArray<T>& x)
    {
        Gadgetron::math::scal(x.get_number_of_elements(), a, x.begin());
        return true;
    }

    template EXPORTCPUCOREMATH bool scal(float a, hoNDArray<float>& x);
    template EXPORTCPUCOREMATH bool scal(double a, hoNDArray<double>& x);
    template EXPORTCPUCOREMATH bool scal(GT_Complex8 a, hoNDArray<GT_Complex8>& x);
    template EXPORTCPUCOREMATH bool scal(GT_Complex16 a, hoNDArray<GT_Complex16>& x);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool scal(T a, hoNDArray< std::complex<T> >& x)
    {
        Gadgetron::math::scal(x.get_number_of_elements(), a, x.begin());
        return true;
    }

    template EXPORTCPUCOREMATH bool scal(float a, hoNDArray<GT_Complex8>& x);
    template EXPORTCPUCOREMATH bool scal(double a, hoNDArray<GT_Complex16>& x);

    // --------------------------------------------------------------------------------

    template <typename T> 
    bool sort(const hoNDArray<T>& x, hoNDArray<T>& r, bool isascending)
    {
        if ( &r != &x )
        {
            if ( r.get_number_of_elements()!=x.get_number_of_elements())
            {
                r = x;
            }
            else
            {
                memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
            }
        }

        Gadgetron::math::sort(x.get_number_of_elements(), x.begin(), r.begin(), isascending);

        return true;
    }

    template EXPORTCPUCOREMATH bool sort(const hoNDArray<float>& x, hoNDArray<float>& r, bool isascending);
    template EXPORTCPUCOREMATH bool sort(const hoNDArray<double>& x, hoNDArray<double>& r, bool isascending);

    // --------------------------------------------------------------------------------
}
