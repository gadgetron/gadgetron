/** \file   hoNDArray_math_util.hxx
    \brief  Implementation of some hoNDArray math functions
*/

#pragma once

namespace Gadgetron
{
    template<typename T> inline void fill( size_t N, T* pX, T val)
    {
        long long n;
        #pragma omp parallel for default(none) private(n) shared(N, pX, val)
        for ( n=0; n<(long long)N; n++ )
        {
            pX[n] = val;
        }
    }

	template<typename T> inline void fill( hoNDArray<T>* x, T val)
    {
        size_t N = x->get_number_of_elements();
        T* pX = x->begin();
        Gadgetron::fill(N, pX, val);
    }

    template<typename T> inline void fill( hoNDArray<T>& x, T val )
    {
        size_t N = x.get_number_of_elements();
        T* pX = x.begin();
        Gadgetron::fill(N, pX, val);
    }

    template<typename T, unsigned int D> inline void fill( hoNDImage<T, D>* x, T val )
    {
        size_t N = x->get_number_of_elements();
        T* pX = x->begin();
        Gadgetron::fill(N, pX, val);
    }

    template<typename T, unsigned int D> inline void fill( hoNDImage<T, D>& x, T val )
    {
        size_t N = x.get_number_of_elements();
        T* pX = x.begin();
        Gadgetron::fill(N, pX, val);
    }

    template<typename T> inline void clear( hoNDArray<T>* x )
    {
        if ( x->get_number_of_elements() > 0 )
        {
            memset( x->get_data_ptr(), 0, x->get_number_of_elements()*sizeof(T));
        }
    }

    template<typename T> inline void clear( hoNDArray<T>& x )
    {
        if ( x.get_number_of_elements() > 0 )
        {
            memset( x.get_data_ptr(), 0, x.get_number_of_elements()*sizeof(T));
        }
    }

    template<typename T, unsigned int D> inline void clear( hoNDImage<T, D>* x )
    {
        if ( x->get_number_of_elements() > 0 )
        {
            memset( x->get_data_ptr(), 0, x->get_number_of_elements()*sizeof(T));
        }
    }

    template<typename T, unsigned int D> inline void clear( hoNDImage<T, D>& x )
    {
        if ( x.get_number_of_elements() > 0 )
        {
            memset( x.get_data_ptr(), 0, x.get_number_of_elements()*sizeof(T));
        }
    }

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

    template<> 
    inline bool complex_to_real_imag(const hoNDArray<float>& cplx, hoNDArray<float>& real, hoNDArray<float>& imag)
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

    template<> 
    inline bool complex_to_real_imag(const hoNDArray<double>& cplx, hoNDArray<double>& real, hoNDArray<double>& imag)
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

    template <typename T> 
    inline bool absolute(const hoNDArray< std::complex<T> >& x, hoNDArray< std::complex<T> >& r)
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

#ifndef USE_MKL

    template <typename T> 
    bool add(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
            if ( r.get_number_of_elements()!=x.get_number_of_elements())
            {
                r = x;
            }

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            if ( pR == pX )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<N; n++ )
                {
                    pR[n] += pY[n];
                }
            }
            else if ( pR == pY )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<N; n++ )
                {
                    pR[n] += pX[n];
                }
            }
            else
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<N; n++ )
                {
                    pR[n] = pX[n] + pY[n];
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in add(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool add(size_t N, const T* x, const T* y, T* r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && r!=NULL);

            long long n;

            const T* pX = x;
            const T* pY = y;
            T* pR = r;

            if ( pR == pX )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] += pY[n];
                }
            }
            else if ( pR == pY )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] += pX[n];
                }
            }
            else
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] = pX[n] + pY[n];
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in add(size_t N, const T* x, const T* y, T* r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool subtract(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
            if ( r.get_number_of_elements()!=x.get_number_of_elements())
            {
                r = x;
            }

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            if ( pR == pX )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] -= pY[n];
                }
            }
            else
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] = pX[n] - pY[n];
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in subtract(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool subtract(size_t N, const T* x, const T* y, T* r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && r!=NULL);

            long long n;

            const T* pX = x;
            const T* pY = y;
            T* pR = r;

            if ( pR == pX )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] -= pY[n];
                }
            }
            else if ( pR == pY )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] -= pX[n];
                }
            }
            else
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<N; n++ )
                {
                    pR[n] = pX[n] - pY[n];
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in subtract(size_t N, const T* x, const T* y, T* r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool multiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
            if ( r.get_number_of_elements()!=x.get_number_of_elements())
            {
                r = x;
            }

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            if ( pR == pX )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] *= pY[n];
                }
            }
            else if ( pR == pY )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] *= pX[n];
                }
            }
            else
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] = pX[n] * pY[n];
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in multiply(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool multiply(size_t N, const T* x, const T* y, T* r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && r!=NULL);

            long long n;

            const T* pX = x;
            const T* pY = y;
            T* pR = r;

            if ( pR == pX )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] *= pY[n];
                }
            }
            else if ( pR == pY )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] *= pX[n];
                }
            }
            else
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<(long long)N; n++ )
                {
                    pR[n] = pX[n] * pY[n];
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in multiply(size_t N, const T* x, const T* y, T* r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool divide(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
            if ( r.get_number_of_elements()!=x.get_number_of_elements())
            {
                r = x;
            }

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            if ( pR == pX )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<N; n++ )
                {
                    pR[n] /= pY[n];
                }
            }
            else
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<N; n++ )
                {
                    pR[n] = pX[n]/pY[n];
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in divide(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool sqrt(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            if ( r.get_number_of_elements()!=x.get_number_of_elements())
            {
                r.create(x.get_dimensions());
            }

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            T* pR = r.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pR)
            for ( n=0; n<N; n++ )
            {
                pR[n] = std::sqrt(pX[n]);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in sqrt(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool minAbsolute(const hoNDArray<T>& x, T& r, size_t& ind)
    {
        try
        {
            long long N = (long long)x.get_number_of_elements();
            if ( N == 0 ) return true;

            long long n;

            const T* pX = x.begin();

            typename realType<T>::Type v = std::abs(pX[0]);
            typename realType<T>::Type v2;

            ind = 0;

            for ( n=1; n<N; n++ )
            {
                v2 = std::abs(pX[n]);
                if ( v2 < v )
                {
                    v = v2;
                    ind = n;
                }
            }

            r = pX[ind];
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in minAbsolute(const hoNDArray<T>& x, T& r, size_t& ind) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool maxAbsolute(const hoNDArray<T>& x, T& r, size_t& ind)
    {
        try
        {
            long long N = (long long)x.get_number_of_elements();
            if ( N == 0 ) return true;

            long long n;

            const T* pX = x.begin();

            typename realType<T>::Type v = std::abs(pX[0]);
            typename realType<T>::Type v2;

            ind = 0;

            for ( n=1; n<N; n++ )
            {
                v2 = std::abs(pX[n]);
                if ( v2 > v )
                {
                    v = v2;
                    ind = n;
                }
            }

            r = pX[ind];
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in maxAbsolute(const hoNDArray<T>& x, T& r, size_t& ind) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool multiplyConj(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
            if ( r.get_number_of_elements()!=x.get_number_of_elements())
            {
                r = x;
            }

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
            for ( n=0; n<N; n++ )
            {
                pR[n] = pX[n] * std::conj(pY[n]);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in multiplyConj(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool conjugate(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            if ( r.get_number_of_elements()!=x.get_number_of_elements())
            {
                r.create(x.get_dimensions());
            }

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            T* pR = r.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pR)
            for ( n=0; n<N; n++ )
            {
                pR[n] = std::conj(pX[n]);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in conjugate(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool addEpsilon(hoNDArray<T>& x)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            T* pX = x.begin();

            typename realType<T>::Type eps = std::numeric_limits<typename realType<T>::Type>::epsilon();

            long long i;

            #pragma omp parallel for default(none) private(i) shared(n, pX, eps)
            for (i=0; i<(long long)n; i++ )
            {
                if ( std::abs(pX[i]) < eps )
                {
                    pX[i] += eps;
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in addEpsilon(hoNDArray<T>& x) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool norm2(const hoNDArray<T>& x, typename realType<T>::Type& r)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            const T* pX = x.begin();

            typename realType<T>::Type sqrNormSum(0), v;

            long long i;

            #pragma omp parallel for default(none) private(i, v) shared(n, pX) reduction(+:sqrNormSum)
            for (i=0; i<(long long)n; i++ )
            {
                v = std::abs(pX[n]);
                sqrNormSum = sqrNormSum + v*v;
            }

            r = std::sqrt(sqrNormSum);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in norm2(const hoNDArray<T>& x, typename realType<T>::Type& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool norm1(const hoNDArray<T>& x, typename realType<T>::Type& r)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            const T* pX = x.begin();

            typename realType<T>::Type norm1Sum(0), v;

            long long i;

            #pragma omp parallel for default(none) private(i, v) shared(n, pX) reduction(+:norm1Sum)
            for (i=0; i<(long long)n; i++ )
            {
                norm1Sum = norm1Sum + std::abs(pX[n]);
            }

            r = norm1Sum;
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in norm1(const hoNDArray<T>& x, typename realType<T>::Type& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool dotc(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();
            r = 0;

            T res(0);

            // #pragma omp parallel for default(none) private(n) shared(N, pX, pY) reduction(+:res)
            for ( n=0; n<N; n++ )
            {
                res += std::conj(pX[n]) *pY[n];
            }

            r = res;
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in dotc(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool absolute(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        long long N = (long long)x.get_number_of_elements();
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r)
        for ( n=0; n<N; n++ )
        {
            r(n) = std::abs( x(n) );
        }

        return true;
    }

    template <typename T> 
    bool argument(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        long long N = (long long)x.get_number_of_elements();
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r)
        for ( n=0; n<N; n++ )
        {
            r(n) = std::arg( x(n) );
        }

        return true;
    }

    template <typename T> 
    bool inv(const hoNDArray<T>& x, hoNDArray<T>& r)
    {
        try
        {
            if ( !r.dimensions_equal(&x) )
            {
                r = x;
            }

            const T* pX = x.begin();
            T* pR = r.begin();

            T v(1.0);
            long long n = x.get_number_of_elements();
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(n, pX, pR, v)
            for ( ii=0; ii<n; ii++ )
            {
                pR[ii] = v/pX[ii];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<T>& x, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t num = x.get_number_of_elements() / (RO*E1);

            size_t kRO = y.get_size(0);
            size_t kE1 = y.get_size(1);

            size_t halfKRO = kRO/2;
            size_t halfKE1 = kE1/2;

            hoNDArray<T> flipY(2*halfKRO+1, 2*halfKE1+1);
            Gadgetron::clear(flipY);

            T* pKer = flipY.begin();

            long long n;
            size_t ro, e1;

            // flip the kernel
            for ( e1=0; e1<kE1; e1++ )
            {
                size_t flip_e1 = 2*halfKE1 - e1;

                for ( ro=0; ro<kRO; ro++ )
                {
                    size_t flip_ro = 2*halfKRO - ro;

                    flipY(flip_ro, flip_e1) = y(ro, e1);
                }
            }

            // perform the convolution
            #pragma omp parallel for default(none) private(n, ro, e1) shared(num, x, RO, E1, z, halfKRO, halfKE1, pKer)
            for ( n=0; n<num; n++ )
            {
                const T* pX = x.begin() + n*RO*E1;
                T* pZ = z.begin() + n*RO*E1;

                size_t kro, ke1, dro, de1;

                for ( e1=0; e1<E1; e1++ )
                {
                    for ( ro=0; ro<RO; ro++ )
                    {
                        pZ[ro + e1*RO] = 0;
                        for ( ke1=-halfKE1; ke1<=halfKE1; ke1++ )
                        {
                            de1 = ke1 + e1;
                            if ( de1 < 0 )
                            {
                                de1 += E1;
                            }
                            else if ( de1 >= E1 )
                            {
                                de1 -= E1;
                            }

                            for ( kro=-halfKRO; kro<=halfKRO; kro++ )
                            {
                                dro = kro + ro;
                                if ( dro < 0 )
                                {
                                    dro += RO;
                                }
                                else if ( dro >= RO )
                                {
                                    dro -= RO;
                                }

                                pZ[ro + e1*RO] += pKer[ kro+halfKRO + (ke1+halfKE1) * (2*halfKRO+1) ] * pX[dro + de1*RO];
                            }
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) ... ");
            return false;
        }

        return true;
    }

    template<typename T> 
    bool conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z)
    {
        try
        {
            if ( !z.dimensions_equal(&x) )
            {
                z = x;
            }

            size_t RO = x.get_size(0);
            size_t E1 = x.get_size(1);
            size_t E2 = x.get_size(2);
            size_t num = x.get_number_of_elements() / (RO*E1*E2);

            size_t kRO = y.get_size(0);
            size_t kE1 = y.get_size(1);
            size_t kE2 = y.get_size(2);

            size_t halfKRO = kRO/2;
            size_t halfKE1 = kE1/2;
            size_t halfKE2 = kE2/2;

            hoNDArray<T> flipY(2*halfKRO+1, 2*halfKE1+1, 2*halfKE2+1);
            Gadgetron::clear(flipY);

            T* pKer = flipY.begin();

            long long n, e2;
            size_t ro, e1;

            // flip the kernel
            for ( e2=0; e2<kE2; e2++ )
            {
                size_t flip_e2 = 2*halfKE2 - e2;

                for ( e1=0; e1<kE1; e1++ )
                {
                    size_t flip_e1 = 2*halfKE1 - e1;

                    for ( ro=0; ro<kRO; ro++ )
                    {
                        size_t flip_ro = 2*halfKRO - ro;

                        flipY(flip_ro, flip_e1, flip_e2) = y(ro, e1, e2);
                    }
                }
            }

            // perform the convolution
            #pragma omp parallel for default(none) private(n, ro, e1, e2) shared(num, x, RO, E1, E2, z, halfKRO, halfKE1, halfKE2, pKer) if ( num > 8 )
            for ( n=0; n<num; n++ )
            {
                const T* pX = x.begin() + n*RO*E1*E2;
                T* pZ = z.begin() + n*RO*E1*E2;

                size_t kro, ke1, ke2, dro, de1, de2;

                #pragma omp parallel for default(none) private(ro, e1, e2, kro, ke1, ke2, dro, de1, de2) shared(pX, RO, E1, E2, pZ, halfKRO, halfKE1, halfKE2, pKer)
                for ( e2=0; e2<(long long)E2; e2++ )
                {
                    for ( e1=0; e1<E1; e1++ )
                    {
                        for ( ro=0; ro<RO; ro++ )
                        {
                            pZ[ro + e1*RO + e2*RO*E1] = 0;
                            for ( ke2=-halfKE2; ke2<=halfKE2; ke2++ )
                            {
                                de2 = ke2 + e2;
                                if ( de2 < 0 )
                                {
                                    de2 += E2;
                                }
                                else if ( de2 >= E2 )
                                {
                                    de2 -= E2;
                                }

                                for ( ke1=-halfKE1; ke1<=halfKE1; ke1++ )
                                {
                                    de1 = ke1 + e1;
                                    if ( de1 < 0 )
                                    {
                                        de1 += E1;
                                    }
                                    else if ( de1 >= E1 )
                                    {
                                        de1 -= E1;
                                    }

                                    for ( kro=-halfKRO; kro<=halfKRO; kro++ )
                                    {
                                        dro = kro + ro;
                                        if ( dro < 0 )
                                        {
                                            dro += RO;
                                        }
                                        else if ( dro >= RO )
                                        {
                                            dro -= RO;
                                        }

                                        pZ[ro + e1*RO + e2*RO*E1] += pKer[ kro+halfKRO + (ke1+halfKE1)*(2*halfKRO+1) + (ke2+halfKE2)*(2*halfKRO+1)*(2*halfKE1+1) ] * pX[dro + de1*RO + de2*RO*E1];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) ... ");
            return false;
        }

        return true;
    }

#endif // USE_MKL

#ifndef USE_MKL

    // ----------------------------------------------------

    template <typename T> 
    T dotc(const hoNDArray<T>& x, const hoNDArray<T>& y)
    {
        T r = 0;

        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pY, r) reductions(+:r)
            for ( n=0; n<N; n++ )
            {
                r = r + std::conj(pX[n]) *pY[n];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in dotc(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r) ... ");
            return 0;
        }

        return r;
    }

    template <typename T> 
    T dotu(const hoNDArray<T>& x, const hoNDArray<T>& y)
    {
        T r = 0;

        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pY, r) reductions(+:r)
            for ( n=0; n<N; n++ )
            {
                r = r + pX[n]*pY[n];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in dotu(const hoNDArray<T>& x, const hoNDArray<T>& y, T& r) ... ");
            return 0;
        }

        return r;
    }

    template <typename T> 
    bool axpy(T a, const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r)
    {
        try
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

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();
            T* pR = r.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR, a)
            for ( n=0; n<N; n++ )
            {
                pR[n] = a*pX[n] + pY[n];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in axpy(T a, const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool scal(T a, hoNDArray<T>& x)
    {
        try
        {
            long long N = (long long)x.get_number_of_elements();
            long long n;

            T* pX = x.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, a)
            for ( n=0; n<N; n++ )
            {
                pX[n] *= a;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(T a, hoNDArray<T>& x) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool scal(T a, hoNDArray< std::complex<T> >& x)
    {
        try
        {
            long long N = (long long)x.get_number_of_elements();
            long long n;

            std::complex<T>* pX = x.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, a)
            for ( n=0; n<N; n++ )
            {
                pX[n] *= a;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(T a, hoNDArray< std::complex<T> >& x) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool scal(T a, T*x, long long N)
    {
        try
        {
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, a)
            for ( n=0; n<N; n++ )
            {
                x[n] *= a;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(T a, T*x, long long N) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    bool scal(T a, std::complex<T>*x, long long N)
    {
        try
        {
            long long N = (long long)x.get_number_of_elements();
            long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, a)
            for ( n=0; n<N; n++ )
            {
                x[n] *= a;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(T a, std::complex<T>*x, long long N) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool scal(T a, hoNDImage<T, D>& x)
    {
        try
        {
            long long N = (long long)x.get_number_of_elements();
            long long n;

            T* pX = x.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, a)
            for ( n=0; n<N; n++ )
            {
                pX[n] *= a;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(T a, hoNDImage<T, D>& x) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool scal(T a, hoNDImage< std::complex<T>, D>& x)
    {
        try
        {
            long long N = (long long)x.get_number_of_elements();
            long long n;

            std::complex<T>* pX = x.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, a)
            for ( n=0; n<N; n++ )
            {
                pX[n] *= a;
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(T a, hoNDImage< std::complex<T>, D>& x) ... ");
            return false;
        }

        return true;
    }

    template <typename T> 
    struct hoNDArrayCompAscending
    {
        bool operator() (T a, T b) { return (a>=b); }
    };

    template <typename T> 
    struct hoNDArrayCompDescending
    {
        bool operator() (T a, T b) { return (a<b); }
    };

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

        if ( isascending )
        {
            hoNDArrayCompAscending<T> obj;
            std::sort(r.begin(), r.end(), obj);
        }
        else
        {
            hoNDArrayCompDescending<T> obj;
            std::sort(r.begin(), r.end(), obj);
        }

        return true;
    }

#endif // USE_MKL

#ifdef USE_MKL

    // ----------------------------------------------------------------------------------------
    // float
    // ----------------------------------------------------------------------------------------

    EXPORTCPUCOREMATH inline bool add(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsAdd(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool subtract(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsSub(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool multiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsMul(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool divide(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsDiv(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool absolute(const hoNDArray<float>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsAbs(x.get_number_of_elements(), x.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool argument(const hoNDArray<float>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        memset(r.begin(), 0, r.get_number_of_bytes());

        return true;
    }

    EXPORTCPUCOREMATH inline bool sqrt(const hoNDArray<float>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vsSqrt(x.get_number_of_elements(), x.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool minAbsolute(const hoNDArray<float>& x, float& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(isamin(&n, x.begin(), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool maxAbsolute(const hoNDArray<float>& x, float& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(isamax(&n, x.begin(), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool addEpsilon(hoNDArray<float>& x)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            float* pX = x.begin();

            long long i;

            #pragma omp parallel for default(none) private(i) shared(n, pX)
            for (i=0; i<(long long)n; i++ )
            {
                if ( GT_ABS(pX[i]) < FLT_EPSILON )
                {
                    pX[i] += GT_SGN(pX[i])*FLT_EPSILON;
                }
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool norm2(const hoNDArray<float>& x, float& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = snrm2(&n, x.begin(), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool norm1(const hoNDArray<float>& x, float& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = sasum(&n, x.begin(), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool inv(const hoNDArray<float>& x, hoNDArray<float>& r)
    {
        try
        {
            if ( !r.dimensions_equal(&x) )
            {
                r = x;
            }

            long long n = x.get_number_of_elements();
            vsInv(n, x.begin(), r.begin());
            GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<float>& x, hoNDArray<float>& r) ... ");
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------------
    // double
    // ----------------------------------------------------------------------------------------

    EXPORTCPUCOREMATH inline bool add(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdAdd(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool subtract(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdSub(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool multiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdMul(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool divide(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdDiv(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool absolute(const hoNDArray<double>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdAbs(x.get_number_of_elements(), x.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool argument(const hoNDArray<double>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        memset(r.begin(), 0, r.get_number_of_bytes());

        return true;
    }

    EXPORTCPUCOREMATH inline bool sqrt(const hoNDArray<double>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vdSqrt(x.get_number_of_elements(), x.begin(), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool minAbsolute(const hoNDArray<double>& x, double& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(idamin(&n, x.begin(), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool maxAbsolute(const hoNDArray<double>& x, double& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(idamax(&n, x.begin(), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool addEpsilon(hoNDArray<double>& x)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            double* pX = x.begin();

            long long i;

            #pragma omp parallel for default(none) private(i) shared(n, pX)
            for (i=0; i<(long long)n; i++ )
            {
                if ( GT_ABS(pX[i]) < DBL_EPSILON )
                {
                    pX[i] += GT_SGN(pX[i])*DBL_EPSILON;
                }
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool norm2(const hoNDArray<double>& x, double& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = dnrm2(&n, x.begin(), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool norm1(const hoNDArray<double>& x, double& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = dasum(&n, x.begin(), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool inv(const hoNDArray<double>& x, hoNDArray<double>& r)
    {
        try
        {
            if ( !r.dimensions_equal(&x) )
            {
                r = x;
            }

            long long n = x.get_number_of_elements();
            vdInv(n, x.begin(), r.begin());
            GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<double>& x, hoNDArray<double>& r) ... ");
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------------
    // GT_Complex8
    // ----------------------------------------------------------------------------------------

    EXPORTCPUCOREMATH inline bool add(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        return true;
    }

    EXPORTCPUCOREMATH inline bool add(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && z!=NULL);
        vcAdd(N, reinterpret_cast<const MKL_Complex8*>(x), reinterpret_cast<const MKL_Complex8*>(y), reinterpret_cast<MKL_Complex8*>(r));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        return true;
    }

    EXPORTCPUCOREMATH inline bool subtract(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcSub(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        return true;
    }

    EXPORTCPUCOREMATH inline bool subtract(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && z!=NULL);
        vcSub(N, reinterpret_cast<const MKL_Complex8*>(x), reinterpret_cast<const MKL_Complex8*>(y), reinterpret_cast<MKL_Complex8*>(r));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        return true;
    }

    EXPORTCPUCOREMATH inline bool multiply(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool multiply(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && z!=NULL);
        vcMul(N, reinterpret_cast<const MKL_Complex8*>(x), reinterpret_cast<const MKL_Complex8*>(y), reinterpret_cast<MKL_Complex8*>(r));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        return true;
    }

    EXPORTCPUCOREMATH inline bool divide(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcDiv(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool sqrt(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vcSqrt(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool minAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(icamin(&n, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool maxAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(icamax(&n, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool multiplyConj(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vcMulByConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool conjugate(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vcConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool addEpsilon(hoNDArray<GT_Complex8>& x)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            GT_Complex8* pX = x.begin();

            long long i;

            #pragma omp parallel for default(none) private(i) shared(n, pX)
            for (i=0; i<(long long)n; i++ )
            {
                if ( std::abs(pX[i]) < FLT_EPSILON )
                {
                    pX[i] += FLT_EPSILON;
                }
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool norm2(const hoNDArray<GT_Complex8>& x, float& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = scnrm2(&n, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool norm1(const hoNDArray<GT_Complex8>& x, float& r)
    {
        try
        {
            hoNDArray<float> a;
            GADGET_CHECK_RETURN_FALSE(absolute(x, a));
            GADGET_CHECK_RETURN_FALSE(norm1(a, r));
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, GT_Complex8& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            MKL_INT N = x.get_number_of_elements();
            MKL_INT incx(1), incy(1);
            cdotc(reinterpret_cast<MKL_Complex8*>(&r), &N, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx, 
                    reinterpret_cast<const MKL_Complex8*>(y.begin()), &incy);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool absolute(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vcAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool argument(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vcArg(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool inv(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    {
        try
        {
            if ( !r.dimensions_equal(&x) )
            {
                r = x;
            }

            const GT_Complex8* pX = x.begin();
            GT_Complex8* pR = r.begin();

            GT_Complex8 v(1.0);
            long long n = x.get_number_of_elements();
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(n, pX, pR, v)
            for ( ii=0; ii<n; ii++ )
            {
                pR[ii] = v/pX[ii];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r) ... ");
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------------
    // GT_Complex16
    // ----------------------------------------------------------------------------------------

    EXPORTCPUCOREMATH inline bool add(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool add(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && r!=NULL);
        vzAdd(N, reinterpret_cast<const MKL_Complex16*>(x), reinterpret_cast<const MKL_Complex16*>(y), reinterpret_cast<MKL_Complex16*>(r));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        return true;
    }

    EXPORTCPUCOREMATH inline bool subtract(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzSub(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool subtract(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && r!=NULL);
        vzSub(N, reinterpret_cast<const MKL_Complex16*>(x), reinterpret_cast<const MKL_Complex16*>(y), reinterpret_cast<MKL_Complex16*>(r));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        return true;
    }

    EXPORTCPUCOREMATH inline bool multiply(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool multiply(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && r!=NULL);
        vzMul(N, reinterpret_cast<const MKL_Complex16*>(x), reinterpret_cast<const MKL_Complex16*>(y), reinterpret_cast<MKL_Complex16*>(r));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        return true;
    }

    EXPORTCPUCOREMATH inline bool divide(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzDiv(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vzAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        hoNDArray<double> rTmp;
        rTmp.create(x.get_dimensions());

        vzAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), rTmp.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        //GADGET_CHECK_RETURN_FALSE(r.copyFrom(rTmp));
        r.copyFrom(rTmp);

        return true;
    }

    EXPORTCPUCOREMATH inline bool sqrt(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vzSqrt(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool minAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(izamin(&n, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool maxAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind)
    {
        try
        {
            MKL_INT n = x.get_number_of_elements();
            MKL_INT incx = 1;
            ind = (size_t)(izamax(&n, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx));
            r = x.at(ind);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool multiplyConj(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r = x;
        }

        vzMulByConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    bool argument(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vzArg(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), r.begin());
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool conjugate(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        vzConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        return true;
    }

    EXPORTCPUCOREMATH inline bool addEpsilon(hoNDArray<GT_Complex16>& x)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            GT_Complex16* pX = x.begin();

            long long i;

            #pragma omp parallel for default(none) private(i) shared(n, pX)
            for (i=0; i<(long long)n; i++ )
            {
                if ( std::abs(pX[i]) < DBL_EPSILON )
                {
                    pX[i] += DBL_EPSILON;
                }
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool norm2(const hoNDArray<GT_Complex16>& x, double& r)
    {
        try
        {
            MKL_INT incx = 1;
            MKL_INT n = x.get_number_of_elements();
            r = dznrm2(&n, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool norm1(const hoNDArray<GT_Complex16>& x, double& r)
    {
        try
        {
            hoNDArray<double> a;
            GADGET_CHECK_RETURN_FALSE(absolute(x, a));
            GADGET_CHECK_RETURN_FALSE(norm1(a, r));
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, GT_Complex16& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            MKL_INT N = x.get_number_of_elements();
            MKL_INT incx(1), incy(1);
            zdotc(reinterpret_cast<MKL_Complex16*>(&r), &N, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx, 
                    reinterpret_cast<const MKL_Complex16*>(y.begin()), &incy);
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool inv(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    {
        try
        {
            if ( !r.dimensions_equal(&x) )
            {
                r = x;
            }

            const GT_Complex16* pX = x.begin();
            GT_Complex16* pR = r.begin();

            GT_Complex16 v(1.0);
            long long n = x.get_number_of_elements();
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(n, pX, pR, v)
            for ( ii=0; ii<n; ii++ )
            {
                pR[ii] = v/pX[ii];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r) ... ");
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------------
    // other functions
    // ----------------------------------------------------------------------------------------

    EXPORTCPUCOREMATH inline GT_Complex8 dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y)
    {
        if ( x.get_number_of_elements() != y.get_number_of_elements() )
        {
            GADGET_ERROR_MSG("dotc(x, y), inputs have differnet length ...");
            return 0.0;
        }

        MKL_INT N = x.get_number_of_elements();
        MKL_INT incx(1), incy(1);
        GT_Complex8 r;
        cdotc(reinterpret_cast<MKL_Complex8*>(&r), &N, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx, reinterpret_cast<const MKL_Complex8*>(y.begin()), &incy);
        return r;
    }

    EXPORTCPUCOREMATH inline GT_Complex16 dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y)
    {
        if ( x.get_number_of_elements() != y.get_number_of_elements() )
        {
            GADGET_ERROR_MSG("dotc(x, y), inputs have differnet length ...");
            return 0;
        }

        MKL_INT N = x.get_number_of_elements();
        MKL_INT incx(1), incy(1);
        GT_Complex16 r;
        zdotc(reinterpret_cast<MKL_Complex16*>(&r), &N, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx, reinterpret_cast<const MKL_Complex16*>(y.begin()), &incy);
        return r;
    }

    EXPORTCPUCOREMATH inline GT_Complex8 dotu(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y)
    {
        if ( x.get_number_of_elements() != y.get_number_of_elements() )
        {
            GADGET_ERROR_MSG("dotu(x, y), inputs have differnet length ...");
            return 0;
        }

        MKL_INT N = x.get_number_of_elements();
        MKL_INT incx(1), incy(1);
        GT_Complex8 r;
        cdotu(reinterpret_cast<MKL_Complex8*>(&r), &N, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx, reinterpret_cast<const MKL_Complex8*>(y.begin()), &incy);
        return r;
    }

    EXPORTCPUCOREMATH inline GT_Complex16 dotu(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y)
    {
        if ( x.get_number_of_elements() != y.get_number_of_elements() )
        {
            GADGET_ERROR_MSG("dotu(x, y), inputs have differnet length ...");
            return 0;
        }

        MKL_INT N = x.get_number_of_elements();
        MKL_INT incx(1), incy(1);
        GT_Complex16 r;
        zdotu(reinterpret_cast<MKL_Complex16*>(&r), &N, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx, reinterpret_cast<const MKL_Complex16*>(y.begin()), &incy);
        return r;
    }

    EXPORTCPUCOREMATH inline bool axpy(float a, const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    {
        try
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

            MKL_INT N = (MKL_INT)(x.get_number_of_elements());
            const MKL_INT incX(1), incY(1);

            cblas_saxpy (N, a, x.begin(), incX, r.begin(), incY);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in axpy(float a, const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool axpy(double a, const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    {
        try
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

            MKL_INT N = (MKL_INT)(x.get_number_of_elements());
            const MKL_INT incX(1), incY(1);

            cblas_daxpy (N, a, x.begin(), incX, r.begin(), incY);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in axpy(double a, const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool axpy(const GT_Complex8& a, const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    {
        try
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

            MKL_INT N = (MKL_INT)(x.get_number_of_elements());
            const MKL_INT incX(1), incY(1);

            cblas_caxpy (N, &a, x.begin(), incX, r.begin(), incY);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in axpy(const GT_Complex8& a, const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool axpy(const GT_Complex16& a, const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    {
        try
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

            MKL_INT N = (MKL_INT)(x.get_number_of_elements());
            const MKL_INT incX(1), incY(1);

            cblas_zaxpy (N, &a, x.begin(), incX, r.begin(), incY);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in axpy(const GT_Complex16& a, const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(float a, hoNDArray<float>& x)
    {
        try
        {
            cblas_sscal ((MKL_INT)(x.get_number_of_elements()), a, x.begin(), 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(float a, hoNDArray<float>& x) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(double a, hoNDArray<double>& x)
    {
        try
        {
            cblas_dscal ((MKL_INT)(x.get_number_of_elements()), a, x.begin(), 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(double a, hoNDArray<double>& x) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(float a, hoNDArray<GT_Complex8>& x)
    {
        try
        {
            GT_Complex8 alpha = GT_Complex8(a);
            cblas_cscal (x.get_number_of_elements(), &alpha, x.begin(), 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(float a, hoNDArray<GT_Complex8>& x) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(double a, hoNDArray<GT_Complex16>& x)
    {
        try
        {
            GT_Complex16 alpha = GT_Complex16(a);
            cblas_zscal (x.get_number_of_elements(), &alpha, x.begin(), 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(double a, hoNDArray<GT_Complex16>& x) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(GT_Complex8 a, hoNDArray<GT_Complex8>& x)
    {
        try
        {
            cblas_cscal (x.get_number_of_elements(), &a, x.begin(), 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(GT_Complex8 a, hoNDArray<GT_Complex8>& x) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(GT_Complex16 a, hoNDArray<GT_Complex16>& x)
    {
        try
        {
            cblas_zscal (x.get_number_of_elements(), &a, x.begin(), 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(GT_Complex16 a, hoNDArray<GT_Complex16>& x) ... ");
            return false;
        }

        return true;
    }

    // -----------------------

    EXPORTCPUCOREMATH inline bool scal(float a, float*x, long long N)
    {
        try
        {
            cblas_sscal ((MKL_INT)(N), a, x, 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(float a, float*x, long long N) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(double a, double*x, long long N)
    {
        try
        {
            cblas_dscal ((MKL_INT)(N), a, x, 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(double a, double*x, long long N) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(float a, GT_Complex8*x, long long N)
    {
        try
        {
            GT_Complex8 alpha = GT_Complex8(a);
            cblas_cscal (N, &alpha, x, 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(float a, GT_Complex8*x, long long N) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(double a, GT_Complex16*x, long long N)
    {
        try
        {
            GT_Complex16 alpha = GT_Complex16(a);
            cblas_zscal (N, &alpha, x, 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(double a, GT_Complex16*x, long long N) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(GT_Complex8 a, GT_Complex8*x, long long N)
    {
        try
        {
            cblas_cscal (N, &a, x, 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(GT_Complex8 a, GT_Complex8*x, long long N) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool scal(GT_Complex16 a, GT_Complex16*x, long long N)
    {
        try
        {
            cblas_zscal (N, &a, x, 1);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(GT_Complex16 a, GT_Complex16*x, long long N) ... ");
            return false;
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool sort(const hoNDArray<float>& x, hoNDArray<float>& r, bool isascending)
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

        if ( isascending )
        {
            GADGET_CHECK_RETURN_FALSE(LAPACKE_slasrt('I', r.get_number_of_elements(), r.begin())==0);
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(LAPACKE_slasrt('D', r.get_number_of_elements(), r.begin())==0);
        }

        return true;
    }

    EXPORTCPUCOREMATH inline bool sort(const hoNDArray<double>& x, hoNDArray<double>& r, bool isascending)
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

        if ( isascending )
        {
            GADGET_CHECK_RETURN_FALSE(LAPACKE_dlasrt('I', r.get_number_of_elements(), r.begin())==0);
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(LAPACKE_dlasrt('D', r.get_number_of_elements(), r.begin())==0);
        }

        return true;
    }

#endif // USE_MKL
}
