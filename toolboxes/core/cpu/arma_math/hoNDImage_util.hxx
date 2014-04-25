/** \file   hoNDImage_util.hxx
    \brief  operations on the hoNDImage class.
*/

namespace Gadgetron
{
    template<class T, unsigned int D> 
    bool real_imag_to_complex(const hoNDImage<typename realType<T>::Type, D>& real, 
                        const hoNDImage<typename realType<T>::Type, D>& imag, 
                        hoNDImage<T, D>& cplx)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(real.dimensions_equal(imag));

            cplx.createFrom(real);

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

    template<class T, unsigned int D> 
    bool complex_to_real_imag(const hoNDImage<T, D>& cplx, 
                        hoNDImage<typename realType<T>::Type, D>& real, 
                        hoNDImage<typename realType<T>::Type, D>& imag)
    {
        try
        {
            real.createFrom(cplx);
            imag.createFrom(cplx);

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

    template<unsigned int D> 
    bool complex_to_real_imag(const hoNDImage<float, D>& cplx, 
                        hoNDImage<float, D>& real, 
                        hoNDImage<float, D>& imag)
    {
        try
        {
            real.createFrom(cplx);
            imag.createFrom(cplx);

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

    template<unsigned int D> 
    bool complex_to_real_imag(const hoNDImage<double, D>& cplx, 
                        hoNDImage<double, D>& real, 
                        hoNDImage<double, D>& imag)
    {
        try
        {
            real.createFrom(cplx);
            imag.createFrom(cplx);

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

    template<class T, unsigned int D> 
    bool complex_to_real(const hoNDImage<T, D>& cplx, hoNDImage<T, D>& real)
    {
        try
        {
            real.createFrom(cplx);

            const T* pRes = cplx.begin();
            T* pReal = real.begin();

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

    template<class T, unsigned int D> 
    bool complex_to_real(const hoNDImage<T, D>& cplx, 
                        hoNDImage<typename realType<T>::Type, D>& real)
    {
        try
        {
            real.createFrom(cplx);

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

    template<class T, unsigned int D> 
    bool complex_to_imag(const hoNDImage<T, D>& cplx, 
                        hoNDImage<typename realType<T>::Type, D>& imag)
    {
        try
        {
            imag.createFrom(cplx);

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

    template<class T, unsigned int D> 
    bool real_to_complex(const hoNDImage<typename realType<T>::Type, D>& real, hoNDImage<T, D>& cplx)
    {
        try
        {
            cplx.createFrom(real);

            const typename realType<T>::Type* pReal = real.begin();
            T* pRes = cplx.begin();

            size_t N = real.get_number_of_elements();

            long long n;
            #pragma omp parallel for default(none) private(n) shared(N, pRes, pReal)
            for ( n=0; n<(long long)N; n++ )
            {
                pRes[n] = pReal[n];
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors in real_to_complex(...) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool add(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r)
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
            GADGET_ERROR_MSG("Error happened in add(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool subtract(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r)
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
                    pR[n] -= pY[n];
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
            GADGET_ERROR_MSG("Error happened in subtract(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool multiply(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r)
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
                    pR[n] *= pY[n];
                }
            }
            else if ( pR == pY )
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<N; n++ )
                {
                    pR[n] *= pX[n];
                }
            }
            else
            {
                #pragma omp parallel for default(none) private(n) shared(N, pX, pY, pR)
                for ( n=0; n<N; n++ )
                {
                    pR[n] = pX[n] * pY[n];
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in multiply(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool divide(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r)
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
            GADGET_ERROR_MSG("Error happened in divide(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool sqrt(const hoNDImage<T, D>& x, hoNDImage<T, D>& r)
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
            GADGET_ERROR_MSG("Error happened in sqrt(const hoNDImage<T, D>& x, hoNDImage<T, D>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool minAbsolute(const hoNDImage<T, D>& x, T& r, size_t& ind)
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
            GADGET_ERROR_MSG("Error happened in minAbsolute(const hoNDImage<T, D>& x, T& r, size_t& ind) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool maxAbsolute(const hoNDImage<T, D>& x, T& r, size_t& ind)
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
            GADGET_ERROR_MSG("Error happened in maxAbsolute(const hoNDImage<T, D>& x, T& r, size_t& ind) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool multiplyConj(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r)
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
            GADGET_ERROR_MSG("Error happened in multiplyConj(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, hoNDImage<T, D>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool conjugate(const hoNDImage<T, D>& x, hoNDImage<T, D>& r)
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
            GADGET_ERROR_MSG("Error happened in conjugate(const hoNDImage<T, D>& x, hoNDImage<T, D>& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool addEpsilon(hoNDImage<T, D>& x)
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
            GADGET_ERROR_MSG("Error happened in addEpsilon(hoNDImage<T, D>& x) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool norm2(const hoNDImage<T, D>& x, typename realType<T>::Type& r)
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
            GADGET_ERROR_MSG("Error happened in norm2(const hoNDImage<T, D>& x, typename realType<T>::Type& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool norm1(const hoNDImage<T, D>& x, typename realType<T>::Type& r)
    {
        try
        {
            size_t n = x.get_number_of_elements();
            T* pX = x.begin();

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
            GADGET_ERROR_MSG("Error happened in norm1(const hoNDImage<T, D>& x, typename realType<T>::Type& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool dotc(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, T& r)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            long long N = (long long)x.get_number_of_elements();
            long long n;

            const T* pX = x.begin();
            const T* pY = y.begin();
            r = 0;

            T v;

            #pragma omp parallel for default(none) private(n) shared(N, pX, pY) reduction(+:v)
            for ( n=0; n<N; n++ )
            {
                v = v + std::conj(pX[n]) *pY[n];
            }

            r = v;
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Error happened in dotc(const hoNDImage<T, D>& x, const hoNDImage<T, D>& y, T& r) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool absolute(const hoNDImage<T, D>& x, hoNDImage<typename realType<T>::Type, D>& r)
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

    template <typename T, unsigned int D> 
    bool absolute(const hoNDImage< std::complex<T>, D >& x, hoNDImage< std::complex<T>, D >& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        hoNDImage<T, D> rTmp;
        rTmp.create(x.get_dimensions());

        Gadgetron::absolute(x, rTmp);

        r.copyFrom(rTmp);

        return true;
    }

    template <typename T, unsigned int D> 
    bool argument(const hoNDImage<T, D>& x, hoNDImage<typename realType<T>::Type, D>& r)
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

    template <typename T, unsigned int D> 
    bool argument(const hoNDImage< std::complex<T>, D >& x, hoNDImage< std::complex<T>, D >& r)
    {
        if ( r.get_number_of_elements()!=x.get_number_of_elements())
        {
            r.create(x.get_dimensions());
        }

        hoNDImage<T, D> rTmp;
        rTmp.create(x.get_dimensions());

        Gadgetron::argument(x, rTmp);

        r.copyFrom(rTmp);

        return true;
    }

    template <typename T, unsigned int D>
    bool inv(const hoNDImage<T, D>& x, hoNDImage<T, D>& r)
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
            GADGET_ERROR_MSG("Errors happened in inv(const hoNDImage<T, D>& x, hoNDImage<T, D>& r) ... ");
            return false;
        }

        return true;
    }

    template<typename T, typename InterpolatorType, unsigned int D> 
    bool downsampleImage(const hoNDImage<T, D>& in, InterpolatorType& interp, hoNDImage<T, D>& out, float ratio[])
    {
        try
        {
            std::vector<size_t> dim(D);
            in.get_dimensions(dim);

            std::vector<size_t> dim_out(D);

            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                dim_out[ii] = (size_t)(dim[ii]/ratio[ii]);
            }

            return Gadgetron::resampleImage(in, interp, dim_out, out);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in downsampleImage(const hoNDImage<T, D>& in, InterpolatorType& interp, hoNDImage<T, D>& out, float ratio[]) ... ");
            return false;
        }

        return true;
    }

    template<typename T, typename InterpolatorType, unsigned int D> 
    bool upsampleImage(const hoNDImage<T, D>& in, InterpolatorType& interp, hoNDImage<T, D>& out, float ratio[])
    {
        try
        {
            std::vector<size_t> dim(D);
            in.get_dimensions(dim);

            std::vector<size_t> dim_out(D);

            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                dim_out[ii] = (size_t)(dim[ii]*ratio[ii]);
            }

            return Gadgetron::resampleImage(in, interp, dim_out, out);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in upsampleImage(const hoNDImage<T, D>& in, InterpolatorType& interp, hoNDImage<T, D>& out, float ratio[]) ... ");
            return false;
        }

        return true;
    }

    template<typename T, typename InterpolatorType, unsigned int D> 
    bool resampleImage(const hoNDImage<T, D>& in, InterpolatorType& interp, const std::vector<size_t>& dim_out, hoNDImage<T, D>& out)
    {
        try
        {
            typedef typename hoNDImage<T, D>::coord_type coord_type;

            /// get the coordinate parameters
            std::vector<size_t> dim;
            in.get_dimensions(dim);

            std::vector<coord_type> pixelSize;
            in.get_pixel_size(pixelSize);

            std::vector<coord_type> origin;
            in.get_origin(origin);

            typename hoNDImage<T, D>::axis_type axis;
            in.get_axis(axis);

            /// compute new pixel sizes
            std::vector<coord_type> pixelSize_out(D);

            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                if ( dim_out[ii] > 1 )
                {
                    pixelSize_out[ii] = (dim[ii]-1)*pixelSize[ii] / (dim_out[ii]-1);
                }
                else
                {
                    pixelSize_out[ii] = (dim[ii]-1)*pixelSize[ii];
                }
            }

            /// set up the out image
            out.create(dim_out, pixelSize_out, origin, axis);

            /// set up the interpolator
            interp.setArray( const_cast< hoNDImage<T, D>& >(in) );

            /// compute the out image

            size_t N = out.get_number_of_elements();

            if ( D == 2 )
            {
                long long ox = (long long)dim_out[0];
                long long oy = (long long)dim_out[1];

                long long x, y;

                #pragma omp parallel default(none) private(x, y) shared(N, ox, oy, in, out, interp)
                {
                    coord_type px, py, ix_in, iy_in;

                    #pragma omp for 
                    for ( y=0; y<oy; y++ )
                    {
                        for ( x=0; x<ox; x++ )
                        {
                            out.image_to_world( (size_t)x, (size_t)y, px, py);

                            in.world_to_image(px, py, ix_in, iy_in);

                            out( (size_t)(x+y*ox) ) = interp(ix_in, iy_in);
                        }
                    }
                }
            }
            else if ( D == 3 )
            {
                long long ox = (long long)dim_out[0];
                long long oy = (long long)dim_out[1];
                long long oz = (long long)dim_out[2];

                long long x, y, z;

                #pragma omp parallel default(none) private(x, y, z) shared(N, ox, oy, oz, in, out, interp)
                {
                    coord_type ix_in, iy_in, iz_in;
                    coord_type px, py, pz;

                    #pragma omp for 
                    for ( z=0; z<oz; z++ )
                    {
                        for ( y=0; y<oy; y++ )
                        {
                            size_t offset = y*ox + z*ox*oy;

                            for ( x=0; x<ox; x++ )
                            {
                                out.image_to_world( (size_t)x, (size_t)y, (size_t)z, px, py, pz);

                                in.world_to_image(px, py, pz, ix_in, iy_in, iz_in);

                                out( (size_t)(x+offset) ) = interp(ix_in, iy_in, iz_in);
                            }
                        }
                    }
                }
            }
            else if ( D == 4 )
            {
                long long ox = (long long)dim_out[0];
                long long oy = (long long)dim_out[1];
                long long oz = (long long)dim_out[2];
                long long ot = (long long)dim_out[3];

                long long x, y, z, t;

                #pragma omp parallel default(none) private(x, y, z, t) shared(N, ox, oy, oz, ot, in, out, interp)
                {
                    coord_type ix_in, iy_in, iz_in, it_in;
                    coord_type px, py, pz, pt;

                    #pragma omp for 
                    for ( t=0; t<ot; t++ )
                    {
                        for ( z=0; z<oz; z++ )
                        {
                            for ( y=0; y<oy; y++ )
                            {
                                size_t offset = y*ox + z*ox*oy + t*ox*oy*oz;

                                for ( x=0; x<ox; x++ )
                                {
                                    out.image_to_world( (size_t)x, (size_t)y, (size_t)z, (size_t)t, px, py, pz, pt);

                                    in.world_to_image(px, py, pz, pt, ix_in, iy_in, iz_in, it_in);

                                    out( (size_t)(x+offset) ) = interp(ix_in, iy_in, iz_in, it_in);
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                long long n;

                #pragma omp parallel default(none) private(n) shared(N, in, out, interp)
                {
                    std::vector<size_t> ind_o(D);
                    std::vector<coord_type> ind_i(D);

                    std::vector<coord_type> pos(D);

                    #pragma omp for 
                    for ( n=0; n<N; n++ )
                    {
                        out.calculate_index(n, ind_o);
                        out.image_to_world(ind_o, pos);

                        in.world_to_image(pos, ind_i);

                        out(n) = interp(ind_i);
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in resampleImage(const hoNDImage<T, D>& in, InterpolatorType& interp, hoNDImage<T, D>& out, size_t size_out[D]) ... ");
            return false;
        }

        return true;
    }

    template<typename T, typename BoundaryHandlerType, unsigned int D> 
    bool downsampleImageBy2WithAveraging(const hoNDImage<T, D>& in, BoundaryHandlerType& bh, hoNDImage<T, D>& out)
    {
        try
        {
            typedef typename hoNDImage<T, D>::coord_type coord_type;

            bh.setArray( const_cast< hoNDImage<T, D>& >(in) );

            /// get the coordinate parameters
            std::vector<size_t> dim;
            in.get_dimensions(dim);

            std::vector<coord_type> pixelSize;
            in.get_pixel_size(pixelSize);

            std::vector<coord_type> origin;
            in.get_origin(origin);

            typename hoNDImage<T, D>::axis_type axis;
            in.get_axis(axis);

            /// compute out image size and pixel size
            std::vector<size_t> dim_out(D);
            std::vector<coord_type> pixelSize_out(D);

            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                dim_out[ii] = (dim[ii] >> 1);
                pixelSize_out[ii] = 2*pixelSize[ii];
            }

            out.create(dim_out, pixelSize_out, origin, axis);

            if ( D == 2 )
            {
                gt_index_type sx = (gt_index_type)dim_out[0];
                gt_index_type sy = (gt_index_type)dim_out[1];

                T weight = 1.0/5;

                gt_index_type x, y;

                #pragma omp parallel for default(none) private(x, y) shared(sx, sy, bh, out)
                for ( y=0; y<sy; y++ )
                {
                    gt_index_type iy = y<<1;

                    for ( x=0; x<sx; x++ )
                    {
                        gt_index_type ix = x<<1;
                        out( (size_t)(x+y*sx) ) = bh(ix, iy) + ( bh(ix+1, iy) + bh(ix-1, iy) ) + ( bh(ix, iy+1) + bh(ix, iy-1) );
                    }
                }

                Gadgetron::scal(weight, out);
            }
            else if ( D == 3 )
            {
                gt_index_type sx = (gt_index_type)dim_out[0];
                gt_index_type sy = (gt_index_type)dim_out[1];
                gt_index_type sz = (gt_index_type)dim_out[2];

                T weight = 1.0/7;

                gt_index_type x, y, z;

                #pragma omp parallel for default(none) private(x, y, z) shared(sx, sy, sz, bh, out)
                for ( z=0; z<sz; z++ )
                {
                    gt_index_type iz = z<<1;

                    for ( y=0; y<sy; y++ )
                    {
                        gt_index_type iy = y<<1;

                        size_t offset = y*sx + z*sx*sy;

                        for ( x=0; x<sx; x++ )
                        {
                            gt_index_type ix = x<<1;

                            out( (size_t)(x+offset) ) = bh(ix, iy, iz) 
                                        + ( bh(ix+1, iy, iz) + bh(ix-1, iy, iz) ) 
                                        + ( bh(ix, iy+1, iz) + bh(ix, iy-1, iz) )
                                        + ( bh(ix, iy, iz+1) + bh(ix, iy, iz-1) );
                        }
                    }
                }

                Gadgetron::scal(weight, out);
            }
            else if ( D == 4 )
            {
                gt_index_type sx = (gt_index_type)dim_out[0];
                gt_index_type sy = (gt_index_type)dim_out[1];
                gt_index_type sz = (gt_index_type)dim_out[2];
                gt_index_type st = (gt_index_type)dim_out[3];

                T weight = 1.0/9;

                gt_index_type x, y, z, t;

                #pragma omp parallel for default(none) private(x, y, z, t) shared(sx, sy, sz, st, bh, out)
                for ( t=0; t<st; t++ )
                {
                    gt_index_type it = t<<1;

                    for ( z=0; z<sz; z++ )
                    {
                        gt_index_type iz = z<<1;

                        for ( y=0; y<sy; y++ )
                        {
                            gt_index_type iy = y<<1;

                            size_t offset = y*sx + z*sx*sy + t*sx*sy*sz;

                            for ( x=0; x<sx; x++ )
                            {
                                gt_index_type ix = x<<1;

                                out( (size_t)(x+offset) ) = bh(ix, iy, iz, it) 
                                            + ( bh(ix+1, iy, iz, it) + bh(ix-1, iy, iz, it) ) 
                                            + ( bh(ix, iy+1, iz, it) + bh(ix, iy-1, iz, it) )
                                            + ( bh(ix, iy, iz+1, it) + bh(ix, iy, iz-1, it) )
                                            + ( bh(ix, iy, iz, it+1) + bh(ix, iy, iz, it-1) );
                            }
                        }
                    }
                }

                Gadgetron::scal(weight, out);
            }
            else
            {
                T weight = 1.0/(2*D+1);

                gt_index_type N = out.get_number_of_elements();

                gt_index_type n;

                #pragma omp parallel default(none) private(n) shared(N, bh, out, dim_out)
                {
                    std::vector<size_t> ind_out(D);
                    std::vector<gt_index_type> ind_in(D);

                    #pragma omp for 
                    for ( n=0; n<N; n++ )
                    {
                        out.calculate_index(n, ind_out);

                        unsigned int ii;
                        for ( ii=0; ii<D; ii++ )
                        {
                            ind_in[ii] = ind_out[ii]<<1;
                        }

                        T v = bh(ind_in);

                        for ( ii=0; ii<D; ii++ )
                        {
                            ind_in[ii]++;
                            v += bh(ind_in);

                            ind_in[ii]--;
                            ind_in[ii]--;
                            v += bh(ind_in);

                            ind_in[ii]++;
                        }

                        out(n) = v;
                    }
                }

                Gadgetron::scal(weight, out);
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in downsampleImageBy2WithAveraging(const hoNDImage<T, D>& in, hoNDImage<T, D>& out) ... ");
            return false;
        }

        return true;
    }

    template<typename T, typename BoundaryHandlerType, unsigned int D> 
    bool expandImageBy2(const hoNDImage<T, D>& in, BoundaryHandlerType& bh, hoNDImage<T, D>& out)
    {
        try
        {
            typedef typename hoNDImage<T, D>::coord_type coord_type;

            bh.setArray( const_cast< hoNDImage<T, D>& >(in) );

            /// get the coordinate parameters
            std::vector<size_t> dim;
            in.get_dimensions(dim);

            std::vector<coord_type> pixelSize;
            in.get_pixel_size(pixelSize);

            std::vector<coord_type> origin;
            in.get_origin(origin);

            typename hoNDImage<T, D>::axis_type axis;
            in.get_axis(axis);

            /// compute out pixel size
            std::vector<coord_type> pixelSize_out(D);

            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                pixelSize_out[ii] = pixelSize[ii]* (coord_type)0.5;
            }

            out.set_pixel_size(pixelSize_out);
            out.set_origin(origin);
            out.set_axis(axis);

            if ( D == 2 )
            {
                gt_index_type sx = (gt_index_type)dim[0];
                gt_index_type sy = (gt_index_type)dim[1];

                gt_index_type x, y;

                #pragma omp parallel for default(none) private(x, y) shared(sx, sy, bh, out)
                for ( y=0; y<sy; y++ )
                {
                    size_t oy = y<<1;

                    for ( x=0; x<sx; x++ )
                    {
                        size_t ox = x<<1;

                        T p00 = bh(x, y);
                        T p10 = bh(x+1, y);
                        T p01 = bh(x, y+1);
                        T p11 = bh(x+1, y+1);

                        out( ox, oy ) = p00;
                        out( ox+1, oy ) = 0.5*(p00 + p10);
                        out( ox, oy+1 ) = 0.5*(p00 + p01);
                        out( ox+1, oy+1 ) = 0.25*(p00+p10+p01+p11);
                    }
                }

                // if out has odd sizes
                gt_index_type sx_out = (gt_index_type)out.get_size(0);
                gt_index_type sy_out = (gt_index_type)out.get_size(1);

                if ( (2*sx) < sx_out )
                {
                    for ( y=0; y<sy_out; y++ )
                    {
                        size_t offset = y*sx_out + sx_out-1;
                        out(offset) = out(offset-1);
                    }
                }

                if ( (2*sy) < sy_out )
                {
                    memcpy(out.begin()+(sy_out-1)*sx_out, out.begin()+(sy_out-2)*sx_out, sizeof(T)*sx_out);
                }
            }
            else if ( D == 3 )
            {
                gt_index_type sx = (gt_index_type)dim[0];
                gt_index_type sy = (gt_index_type)dim[1];
                gt_index_type sz = (gt_index_type)dim[2];

                gt_index_type x, y, z;

                #pragma omp parallel for default(none) private(x, y, z) shared(sx, sy, sz, bh, out)
                for ( z=0; z<sz; z++ )
                {
                    size_t oz = z<<1;

                    for ( y=0; y<sy; y++ )
                    {
                        size_t oy = y<<1;

                        for ( x=0; x<sx; x++ )
                        {
                            size_t ox = x<<1;

                            T p000 = bh(x, y, z);
                            T p100 = bh(x+1, y, z);
                            T p010 = bh(x, y+1, z);
                            T p110 = bh(x+1, y+1, z);

                            T p001 = bh(x, y, z+1);
                            T p101 = bh(x+1, y, z+1);
                            T p011 = bh(x, y+1, z+1);
                            T p111 = bh(x+1, y+1, z+1);

                            out( ox, oy, oz ) = p000;
                            out( ox+1, oy, oz ) = 0.5*(p000 + p100);
                            out( ox, oy+1, oz ) = 0.5*(p000 + p010);
                            out( ox+1, oy+1, oz ) = 0.25*(p000+p100+p010+p110);

                            out( ox, oy, oz+1 ) = 0.5*(p000 + p001);
                            out( ox+1, oy, oz+1 ) = 0.25*(p000 + p100 + p001 + p101);
                            out( ox, oy+1, oz+1 ) = 0.25*(p000 + p010 + p001 + p011);
                            out( ox+1, oy+1, oz+1 ) = 0.125*(p000+p100+p010+p110+p001+p101+p011+p111);
                        }
                    }
                }

                // if out has odd sizes
                gt_index_type sx_out = (gt_index_type)out.get_size(0);
                gt_index_type sy_out = (gt_index_type)out.get_size(1);
                gt_index_type sz_out = (gt_index_type)out.get_size(2);

                if ( (2*sx) < sx_out )
                {
                    for ( z=0; z<sz_out; z++ )
                    {
                        for ( y=0; y<sy_out; y++ )
                        {
                            size_t offset = y*sx_out + z*sx_out*sy_out;

                            out( size_t(sx_out-1+offset) ) = out( size_t(sx_out-2+offset) );
                        }
                    }
                }

                if ( (2*sy) < sy_out )
                {
                    for ( z=0; z<sz_out; z++ )
                    {
                        size_t offset = z*sx_out*sy_out + (sy_out-1)*sx_out;

                        for ( x=0; x<sx_out; x++ )
                        {
                            out( (size_t)(x+offset) ) = out( (size_t)(x+offset-sx_out) );
                        }
                    }
                }

                if ( (2*sz) < sz_out )
                {
                    memcpy(out.begin()+(sz_out-1)*sx_out*sy_out, out.begin()+(sz_out-2)*sx_out*sy_out, sizeof(T)*sx_out*sy_out);
                }
            }
            else
            {
                hoNDInterpolatorLinear<hoNDImage<T, D> > interp(const_cast< hoNDImage<T, D>& >(in), bh);

                gt_index_type N = (gt_index_type)(out.get_number_of_elements());

                gt_index_type n;

                #pragma omp parallel default(none) private(n) shared(N, bh, in, out, interp)
                {
                    std::vector<size_t> ind_out(D);
                    std::vector<coord_type> ind_in(D);

                    #pragma omp for 
                    for ( n=0; n<N; n++ )
                    {
                        out.calculate_index(n, ind_out);

                        unsigned int ii;
                        for ( ii=0; ii<D; ii++ )
                        {
                            ind_in[ii] = (coord_type)(ind_out[ii]*0.5);
                        }

                        out( (size_t)(n) ) = interp(ind_in);
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in expandImageBy2(const hoNDImage<T, D>& in, BoundaryHandlerType& bh, hoNDImage<T, D>& out) ... ");
            return false;
        }

        return true;
    }

    template<class ArrayType> 
    bool filterMedian(const ArrayType& img, unsigned int w[], ArrayType& img_out)
    {
        try
        {
            typedef typename ArrayType::value_type T;

            size_t D = img.get_number_of_dimensions();

            img_out = img;

            if ( D == 1 )
            {
                long long halfW = w[0]/2;
                long long N = (long long)img.get_number_of_elements();

                long long n, m, t;

                #pragma omp parallel default(none) private(n, m, t) shared(halfW, N, img, img_out)
                {
                    std::vector<T> buf(2*halfW+1);

                    #pragma omp for 
                    for ( n=0; n<N; n++ )
                    {
                        for ( m=-halfW; m<=halfW; m++ )
                        {
                            t = n + m;
                            if ( t<0 ) t = 0;
                            if ( t > N-1 ) t = N-1;
                            buf[m+halfW] = img( (size_t)t );
                        }

                        std::sort(buf.begin(), buf.end());

                        img_out(n) = buf[halfW];
                    }
                }
            }
            else if ( D == 2 )
            {
                long long halfX = w[0]/2;
                long long halfY = w[1]/2;
                long long sx = (long long)img.get_size(0);
                long long sy = (long long)img.get_size(1);

                const T* pImg = img.begin();
                T* pImgOut = img_out.begin();

                long long WX = 2*halfX+1;
                long long WY = 2*halfY+1;

                long long medianInd = WX*WY/2;

                long long x, y, tx, ty, hx, hy;
                #pragma omp parallel default(none) private(x, y, tx, ty, hx, hy) shared(halfX, halfY, sx, sy, WX, WY, pImg, pImgOut, medianInd)
                {
                    std::vector<T> buf(WX*WY);

                    #pragma omp for 
                    for ( y=halfY; y<sy-halfY; y++ )
                    {
                        for ( x=halfX; x<sx-halfX; x++ )
                        {
                            size_t ind(0);
                            for ( hy=-halfY; hy<=halfY; hy++ )
                            {
                                ty = hy + y;

                                for ( hx=-halfX; hx<=halfX; hx++ )
                                {
                                    tx = hx + x;

                                    buf[ind++] = pImg[tx + ty*sx];
                                }
                            }

                            std::sort(buf.begin(), buf.end());

                            pImgOut[x + y*sx] = buf[medianInd];
                        }
                    }
                }

                std::vector<T> buf(WX*WY);

                for ( y=0; y<halfY; y++ )
                {
                    for ( x=0; x<sx; x++ )
                    {
                        size_t ind(0);
                        for ( hy=-halfY; hy<=halfY; hy++ )
                        {
                            ty = hy + y;
                            if ( ty < 0 ) ty = 0;

                            for ( hx=-halfX; hx<=halfX; hx++ )
                            {
                                tx = hx + x;
                                if ( tx < 0 ) tx = 0;
                                if ( tx > sx-1 ) tx = sx-1;

                                buf[ind++] = pImg[tx + ty*sx];
                            }
                        }

                        std::sort(buf.begin(), buf.end());

                        pImgOut[x + y*sx] = buf[medianInd];
                    }
                }

                for ( y=sx-halfY; y<sx; y++ )
                {
                    for ( x=0; x<sx; x++ )
                    {
                        size_t ind(0);
                        for ( hy=-halfY; hy<=halfY; hy++ )
                        {
                            ty = hy + y;
                            if ( ty > sy-1 ) ty = sy-1;

                            for ( hx=-halfX; hx<=halfX; hx++ )
                            {
                                tx = hx + x;
                                if ( tx < 0 ) tx = 0;
                                if ( tx > sx-1 ) tx = sx-1;

                                buf[ind++] = pImg[tx + ty*sx];
                            }
                        }

                        std::sort(buf.begin(), buf.end());

                        pImgOut[x + y*sx] = buf[medianInd];
                    }
                }
            }
            else if ( D == 3 )
            {
                long long halfX = w[0]/2;
                long long halfY = w[1]/2;
                long long halfZ = w[2]/2;
                long long sx = (long long)img.get_size(0);
                long long sy = (long long)img.get_size(1);
                long long sz = (long long)img.get_size(2);

                const T* pImg = img.begin();
                T* pImgOut = img_out.begin();

                long long WX = 2*halfX+1;
                long long WY = 2*halfY+1;
                long long WZ = 2*halfZ+1;

                long long medianInd = WX*WY*WZ/2;

                long long x, y, z, tx, ty, tz, hx, hy, hz;
                #pragma omp parallel default(none) private(x, y, z, tx, ty, tz, hx, hy, hz) shared(halfX, halfY, halfZ, sx, sy, sz, WX, WY, WZ, pImg, pImgOut, medianInd)
                {
                    std::vector<T> buf(WX*WY*WZ);

                    #pragma omp for 
                    for ( z=halfZ; z<sz-halfZ; z++ )
                    {
                        for ( y=halfY; y<sy-halfY; y++ )
                        {
                            for ( x=halfX; x<sx-halfX; x++ )
                            {
                                size_t ind(0);
                                for ( hz=-halfZ; hz<=halfZ; hz++ )
                                {
                                    tz = hz + z;

                                    for ( hy=-halfY; hy<=halfY; hy++ )
                                    {
                                        ty = hy + y;

                                        for ( hx=-halfX; hx<=halfX; hx++ )
                                        {
                                            tx = hx + x;

                                            buf[ind++] = pImg[tx + ty*sx + tz*sx*sy];
                                        }
                                    }
                                }

                                std::sort(buf.begin(), buf.end());

                                pImgOut[x + y*sx + z*sx*sy] = buf[medianInd];
                            }
                        }
                    }
                }

                std::vector<T> buf(WX*WY*WZ);

                for ( z=0; z<halfZ; z++ )
                {
                    for ( y=0; y<sy; y++ )
                    {
                        for ( x=0; x<sx; x++ )
                        {
                            size_t ind(0);
                            for ( hz=-halfZ; hz<=halfZ; hz++ )
                            {
                                tz = hz + z;
                                if ( tz < 0 ) tz = 0;

                                for ( hy=-halfY; hy<=halfY; hy++ )
                                {
                                    ty = hy + y;
                                    if ( ty < 0 ) ty = 0;
                                    if ( ty > sy-1 ) ty = sy-1;

                                    for ( hx=-halfX; hx<=halfX; hx++ )
                                    {
                                        tx = hx + x;
                                        if ( tx < 0 ) tx = 0;
                                        if ( tx > sx-1 ) tx = sx-1;

                                        buf[ind++] = pImg[tx + ty*sx + tz*sx*sy];
                                    }
                                }
                            }

                            std::sort(buf.begin(), buf.end());

                            pImgOut[x + y*sx + z*sx*sy] = buf[medianInd];
                        }
                    }
                }

                for ( z=sz-halfZ; z<sz; z++ )
                {
                    for ( y=0; y<sy; y++ )
                    {
                        for ( x=0; x<sx; x++ )
                        {
                            size_t ind(0);
                            for ( hz=-halfZ; hz<=halfZ; hz++ )
                            {
                                tz = hz + z;
                                if ( tz > sz-1 ) tz = sz-1;

                                for ( hy=-halfY; hy<=halfY; hy++ )
                                {
                                    ty = hy + y;
                                    if ( ty < 0 ) ty = 0;
                                    if ( ty > sy-1 ) ty = sy-1;

                                    for ( hx=-halfX; hx<=halfX; hx++ )
                                    {
                                        tx = hx + x;
                                        if ( tx < 0 ) tx = 0;
                                        if ( tx > sx-1 ) tx = sx-1;

                                        buf[ind++] = pImg[tx + ty*sx + tz*sx*sy];
                                    }
                                }
                            }

                            std::sort(buf.begin(), buf.end());

                            pImgOut[x + y*sx + z*sx*sy] = buf[medianInd];
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in filterMedian(const ArrayType& img, unsigned int w[], ArrayType& img_out) ... ");
            return false;
        }

        return true;
    }
}
