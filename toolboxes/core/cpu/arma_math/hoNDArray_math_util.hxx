/** \file   hoNDArray_math_util.hxx
    \brief  Implementation of some hoNDArray math functions
*/

#pragma once

#include "hoNDArray_elemwise.h"
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

            /*long long N = (long long)x.get_number_of_elements();
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
            }*/
            Gadgetron::math::add(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
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

            Gadgetron::math::add(N, x, y, r);
            
            /*long long n;

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
            }*/
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

            /*long long N = (long long)x.get_number_of_elements();
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
            }*/
            
            Gadgetron::math::subtract(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
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

            /*long long n;

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
            }*/

            Gadgetron::math::subtract(N, x, y, r);
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

            /*long long N = (long long)x.get_number_of_elements();
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
            }*/

            Gadgetron::math::multiply(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
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

            /*long long n;

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
            }*/

            Gadgetron::math::multiply(N, x, y, r);
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

            /*long long N = (long long)x.get_number_of_elements();
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
            }*/

            Gadgetron::math::divide(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
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

            /*long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            T* pR = r.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pR)
            for ( n=0; n<N; n++ )
            {
                pR[n] = std::sqrt(pX[n]);
            }*/
            Gadgetron::math::sqrt(x.get_number_of_elements(), x.begin(), r.begin());
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
            size_t N = x.get_number_of_elements();
            if ( N == 0 ) return true;

            /*long long n;

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

            r = pX[ind];*/

            Gadgetron::math::minAbsolute(N, x.begin(), r, ind);
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
            size_t N = x.get_number_of_elements();
            if ( N == 0 ) return true;
/*
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
*/
            Gadgetron::math::maxAbsolute(N, x.begin(), r, ind);
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

            //long long N = (long long)x.get_number_of_elements();
            //long long n;

            //const T* pX = x.begin();
            //const T* pY = y.begin();
            //T* pR = r.begin();

            //#pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
            //for ( n=0; n<N; n++ )
            //{
            //    pR[n] = pX[n] * std::conj(pY[n]);
            //}

            Gadgetron::math::multiplyConj(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
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

            /*long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            T* pR = r.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pR)
            for ( n=0; n<N; n++ )
            {
                pR[n] = std::conj(pX[n]);
            }*/

            Gadgetron::math::conjugate(x.get_number_of_elements(), x.begin(), r.begin());
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
            //size_t n = x.get_number_of_elements();
            //T* pX = x.begin();

            //typename realType<T>::Type eps = std::numeric_limits<typename realType<T>::Type>::epsilon();

            //long long i;

            //#pragma omp parallel for default(none) private(i) shared(n, pX, eps)
            //for (i=0; i<(long long)n; i++ )
            //{
            //    if ( std::abs(pX[i]) < eps )
            //    {
            //        pX[i] += eps;
            //    }
            //}

            Gadgetron::math::addEpsilon(x.get_number_of_elements(), x.begin());
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
            //size_t n = x.get_number_of_elements();
            //const T* pX = x.begin();

            //typename realType<T>::Type sqrNormSum(0), v;

            //long long i;

            //#pragma omp parallel for default(none) private(i, v) shared(n, pX) reduction(+:sqrNormSum)
            //for (i=0; i<(long long)n; i++ )
            //{
            //    v = std::abs(pX[i]);
            //    sqrNormSum = sqrNormSum + v*v;
            //}

            //r = std::sqrt(sqrNormSum);

            Gadgetron::math::norm2(x.get_number_of_elements(), x.begin(), r);
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
            /*size_t n = x.get_number_of_elements();
            const T* pX = x.begin();

            typename realType<T>::Type norm1Sum(0), v;

            long long i;

            #pragma omp parallel for default(none) private(i, v) shared(n, pX) reduction(+:norm1Sum)
            for (i=0; i<(long long)n; i++ )
            {
                norm1Sum = norm1Sum + std::abs(pX[i]);
            }

            r = norm1Sum;*/

           Gadgetron::math::norm1(x.get_number_of_elements(), x.begin(), r);
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

            //long long N = (long long)x.get_number_of_elements();
            //long long n;

            //const T* pX = x.begin();
            //const T* pY = y.begin();
            //r = 0;

            //T res(0);

            //// #pragma omp parallel for default(none) private(n) shared(N, pX, pY) reduction(+:res)
            //for ( n=0; n<N; n++ )
            //{
            //    res += std::conj(pX[n]) *pY[n];
            //}

            //r = res;

            Gadgetron::math::dotc(x.get_number_of_elements(), x.begin(), y.begin(), r);
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

        /*long long N = (long long)x.get_number_of_elements();
        long long n;

        #pragma omp parallel for default(none) private(n) shared(N, x, r)
        for ( n=0; n<N; n++ )
        {
            r(n) = std::abs( x(n) );
        }*/

        Gadgetron::math::absolute(x.get_number_of_elements(), x.begin(), r.begin());

        return true;
    }

    template <typename T> 
    bool argument(const hoNDArray<T>& x, hoNDArray<typename realType<T>::Type>& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        //long long N = (long long)x.get_number_of_elements();
        //long long n;

        //#pragma omp parallel for default(none) private(n) shared(N, x, r)
        //for ( n=0; n<N; n++ )
        //{
        //    r(n) = std::arg( x(n) );
        //}

        Gadgetron::math::argument(x.get_number_of_elements(), x.begin(), r.begin());

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

            /*const T* pX = x.begin();
            T* pR = r.begin();

            T v(1.0);
            long long n = x.get_number_of_elements();
            long long ii;

            #pragma omp parallel for default(none) private(ii) shared(n, pX, pR, v)
            for ( ii=0; ii<n; ii++ )
            {
                pR[ii] = v/pX[ii];
            }*/

            Gadgetron::math::inv(x.get_number_of_elements(), x.begin(), r.begin());
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

            long long RO = (long long) x.get_size(0);
            long long E1 = (long long) x.get_size(1);
            long long num = ((long long) x.get_number_of_elements()) / (RO*E1);

            long long kRO = (long long) y.get_size(0);
            long long kE1 = (long long) y.get_size(1);

            Gadgetron::math::conv2(RO, E1, num, x.begin(), kRO, kE1, y.begin(), z.begin());

            //long long halfKRO = kRO/2;
            //long long halfKE1 = kE1/2;

            //hoNDArray<T> flipY(2*halfKRO+1, 2*halfKE1+1);
            //Gadgetron::clear(flipY);

            //T* pKer = flipY.begin();

            //long long n;
            //long long ro, e1;

            //// flip the kernel
            //for ( e1=0; e1<kE1; e1++ )
            //{
            //    long long flip_e1 = 2*halfKE1 - e1;

            //    for ( ro=0; ro<kRO; ro++ )
            //    {
            //        long long flip_ro = 2*halfKRO - ro;

            //        flipY(flip_ro, flip_e1) = y(ro, e1);
            //    }
            //}

            //// perform the convolution
            //#pragma omp parallel for default(none) private(n, ro, e1) shared(num, x, RO, E1, z, halfKRO, halfKE1, pKer)
            //for ( n=0; n<num; n++ )
            //{
            //    const T* pX = x.begin() + n*RO*E1;
            //    T* pZ = z.begin() + n*RO*E1;

            //    long long kro, ke1, dro, de1;

            //    for ( e1=0; e1<E1; e1++ )
            //    {
            //        for ( ro=0; ro<RO; ro++ )
            //        {
            //            pZ[ro + e1*RO] = 0;
            //            for ( ke1=-halfKE1; ke1<=halfKE1; ke1++ )
            //            {
            //                de1 = ke1 + e1;
            //                if ( de1 < 0 )
            //                {
            //                    de1 += E1;
            //                }
            //                else if ( de1 >= E1 )
            //                {
            //                    de1 -= E1;
            //                }

            //                for ( kro=-halfKRO; kro<=halfKRO; kro++ )
            //                {
            //                    dro = kro + ro;
            //                    if ( dro < 0 )
            //                    {
            //                        dro += RO;
            //                    }
            //                    else if ( dro >= RO )
            //                    {
            //                        dro -= RO;
            //                    }

            //                    pZ[ro + e1*RO] += pKer[ kro+halfKRO + (ke1+halfKE1) * (2*halfKRO+1) ] * pX[dro + de1*RO];
            //                }
            //            }
            //        }
            //    }
            //}
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

            long long RO = (long long) x.get_size(0);
            long long E1 = (long long) x.get_size(1);
            long long E2 = (long long) x.get_size(2);
            long long num = ((long long)x.get_number_of_elements()) / (RO*E1*E2);

            long long kRO = (long long) y.get_size(0);
            long long kE1 = (long long) y.get_size(1);
            long long kE2 = (long long) y.get_size(2);

            Gadgetron::math::conv3(RO, E1, E2, num, x.begin(), kRO, kE1, kE2, y.begin(), z.begin());

            //long long halfKRO = kRO/2;
            //long long halfKE1 = kE1/2;
            //long long halfKE2 = kE2/2;

            //hoNDArray<T> flipY(2*halfKRO+1, 2*halfKE1+1, 2*halfKE2+1);
            //Gadgetron::clear(flipY);

            //T* pKer = flipY.begin();

            //long long n, e2;
            //long long ro, e1;

            //// flip the kernel
            //for ( e2=0; e2<kE2; e2++ )
            //{
            //    long long flip_e2 = 2*halfKE2 - e2;

            //    for ( e1=0; e1<kE1; e1++ )
            //    {
            //        long long flip_e1 = 2*halfKE1 - e1;

            //        for ( ro=0; ro<kRO; ro++ )
            //        {
            //            long long flip_ro = 2*halfKRO - ro;

            //            flipY(flip_ro, flip_e1, flip_e2) = y(ro, e1, e2);
            //        }
            //    }
            //}

            //// perform the convolution
            //#pragma omp parallel for default(none) private(n, ro, e1, e2) shared(num, x, RO, E1, E2, z, halfKRO, halfKE1, halfKE2, pKer) if ( num > 8 )
            //for ( n=0; n<num; n++ )
            //{
            //    const T* pX = x.begin() + n*RO*E1*E2;
            //    T* pZ = z.begin() + n*RO*E1*E2;

            //    long long kro, ke1, ke2, dro, de1, de2;

            //    #pragma omp parallel for default(none) private(ro, e1, e2, kro, ke1, ke2, dro, de1, de2) shared(pX, RO, E1, E2, pZ, halfKRO, halfKE1, halfKE2, pKer)
            //    for ( e2=0; e2<E2; e2++ )
            //    {
            //        for ( e1=0; e1<E1; e1++ )
            //        {
            //            for ( ro=0; ro<RO; ro++ )
            //            {
            //                pZ[ro + e1*RO + e2*RO*E1] = 0;
            //                for ( ke2=-halfKE2; ke2<=halfKE2; ke2++ )
            //                {
            //                    de2 = ke2 + e2;
            //                    if ( de2 < 0 )
            //                    {
            //                        de2 += E2;
            //                    }
            //                    else if ( de2 >= E2 )
            //                    {
            //                        de2 -= E2;
            //                    }

            //                    for ( ke1=-halfKE1; ke1<=halfKE1; ke1++ )
            //                    {
            //                        de1 = ke1 + e1;
            //                        if ( de1 < 0 )
            //                        {
            //                            de1 += E1;
            //                        }
            //                        else if ( de1 >= E1 )
            //                        {
            //                            de1 -= E1;
            //                        }

            //                        for ( kro=-halfKRO; kro<=halfKRO; kro++ )
            //                        {
            //                            dro = kro + ro;
            //                            if ( dro < 0 )
            //                            {
            //                                dro += RO;
            //                            }
            //                            else if ( dro >= RO )
            //                            {
            //                                dro -= RO;
            //                            }

            //                            pZ[ro + e1*RO + e2*RO*E1] += pKer[ kro+halfKRO + (ke1+halfKE1)*(2*halfKRO+1) + (ke2+halfKE2)*(2*halfKRO+1)*(2*halfKE1+1) ] * pX[dro + de1*RO + de2*RO*E1];
            //                        }
            //                    }
            //                }
            //            }
            //        }
            //    }
            //}
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<T>& x, const hoNDArray<T>& y, hoNDArray<T>& z) ... ");
            return false;
        }

        return true;
    }

    // ----------------------------------------------------

    template <typename T> 
    T dotc(const hoNDArray<T>& x, const hoNDArray<T>& y)
    {
        T r = 0;

        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            /*long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pY) reduction(+:r)
            for ( n=0; n<N; n++ )
            {
                r = r + std::conj(pX[n]) *pY[n];
            }*/
            dotc(x, y, r);
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

            /*long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, pY) reduction(+:r)
            for ( n=0; n<N; n++ )
            {
                r = r + pX[n]*pY[n];
            }*/

            Gadgetron::math::dotu(x.get_number_of_elements(), x.begin(), y.begin(), r);
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

            //long long N = (long long)x.get_number_of_elements();
            //long long n;

            //const T* pX = x.begin();
            //const T* pY = y.begin();
            //T* pR = r.begin();

            //#pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR, a)
            //for ( n=0; n<N; n++ )
            //{
            //    pR[n] = a*pX[n] + pY[n];
            //}

            Gadgetron::math::axpy(a, x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
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
            /*long long N = (long long)x.get_number_of_elements();
            long long n;

            T* pX = x.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, a)
            for ( n=0; n<N; n++ )
            {
                pX[n] *= a;
            }*/

            Gadgetron::math::scal(x.get_number_of_elements(), a, x.begin());
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
            /*long long N = (long long)x.get_number_of_elements();
            long long n;

            std::complex<T>* pX = x.begin();

            #pragma omp parallel for default(none) private(n) shared(N, pX, a)
            for ( n=0; n<N; n++ )
            {
                pX[n] *= a;
            }*/

            Gadgetron::math::scal(x.get_number_of_elements(), a, x.begin());
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
            /*long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, a)
            for ( n=0; n<N; n++ )
            {
                x[n] *= a;
            }*/

            Gadgetron::math::scal(N, a, x);
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
            /*long long n;

            #pragma omp parallel for default(none) private(n) shared(N, x, a)
            for ( n=0; n<N; n++ )
            {
                x[n] *= a;
            }*/

            Gadgetron::math::scal(N, a, x);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in scal(T a, std::complex<T>*x, long long N) ... ");
            return false;
        }

        return true;
    }

    //template <typename T> 
    //struct hoNDArrayCompAscending
    //{
    //    bool operator() (T a, T b) { return (a>=b); }
    //};

    //template <typename T> 
    //struct hoNDArrayCompDescending
    //{
    //    bool operator() (T a, T b) { return (a<b); }
    //};

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

        /*if ( isascending )
        {
            hoNDArrayCompAscending<T> obj;
            std::sort(r.begin(), r.end(), obj);
        }
        else
        {
            hoNDArrayCompDescending<T> obj;
            std::sort(r.begin(), r.end(), obj);
        }*/

        Gadgetron::math::sort(x.get_number_of_elements(), x.begin(), r.begin(), isascending);

        return true;
    }
}
