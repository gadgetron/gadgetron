#include "hoNDArray_blas.h"

namespace Gadgetron{

    template<class T> T dot( hoNDArray<T> *x, hoNDArray<T> *y, bool cc )
    {
        if( x == 0x0 || y == 0x0 )
            throw std::runtime_error("Gadgetron::dot(): Invalid input array");

        if( x->get_number_of_elements() != y->get_number_of_elements() )
            throw std::runtime_error("Gadgetron::dot(): Array sizes mismatch");

        arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
        arma::Col<typename stdType<T>::Type> yM = as_arma_col(y);
        typename stdType<T>::Type res = (cc) ? arma::cdot(xM,yM) : arma::dot(xM,yM);
        return *((T*)(&res));
    }

    template<class T> typename realType<T>::Type asum( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::asum(): Invalid input array");

        typedef typename realType<T>::Type realT;
        arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
        return realT(arma::norm(xM,1));
    }

    template<class T> T asum( hoNDArray< std::complex<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::asum(): Invalid input array");

        return arma::norm(arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x))),1);
    }

    template<class T> T asum( hoNDArray< complext<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::asum(): Invalid input array");

        return arma::norm(arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x))),1);
    }

    template<class T> typename realType<T>::Type nrm2( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::nrm2(): Invalid input array");

        typedef typename realType<T>::Type realT;
        arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
        return realT(arma::norm(xM,2));
    }

    template<class T> typename realType<T>::Type nrm1( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            BOOST_THROW_EXCEPTION(runtime_error("Gadgetron::nrm1(): Invalid input array"));

        typedef typename realType<T>::Type realT;
        arma::Col<typename stdType<T>::Type> xM = as_arma_col(x);
        return realT(arma::norm(xM,1));
    }

    template<class T> unsigned long long amin( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amin(): Invalid input array");

        typedef typename realType<T>::Type realT;
        arma::Col<realT> xM = arma::abs(as_arma_col(x));
	arma::uword idx;
        realT min = xM.min(idx);
        return idx;
    }

    template<class T> unsigned long long amin( hoNDArray< std::complex<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amin(): Invalid input array");

        arma::Col<T> xM = arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x)));
	arma::uword idx;
        T min = xM.min(idx);
        return idx;
    }

    template<class T> unsigned long long amin( hoNDArray< complext<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amin(): Invalid input array");

        arma::Col<T> xM = arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x)));
	arma::uword idx;
        T min = xM.min(idx);
        return idx;
    }

    template<class T> unsigned long long amax( hoNDArray<T> *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amax(): Invalid input array");

        typedef typename realType<T>::Type realT;
        arma::Col<realT> xM = arma::abs(as_arma_col(x));
	arma::uword idx;
        realT max = xM.max(idx);
        return idx;
    }

    template<class T> unsigned long long amax( hoNDArray< std::complex<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amax(): Invalid input array");

        arma::Col<T> xM = arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x)));
	arma::uword idx;
        T max = xM.max(idx);
        return idx;
    }

    template<class T> unsigned long long amax( hoNDArray< complext<T> > *x )
    {
        if( x == 0x0 )
            throw std::runtime_error("Gadgetron::amax(): Invalid input array");

        arma::Col<T> xM = arma::abs(real(as_arma_col(x)))+arma::abs(imag(as_arma_col(x)));
	arma::uword idx;
        T max = xM.max(idx);
        return idx;
    }

    template<class T> void axpy( T a, hoNDArray<T> *x, hoNDArray<T> *y )
    {
        if( x == 0x0 || y == 0x0 )
            throw std::runtime_error("Gadgetron::axpy(): Invalid input array");

        if( x->get_number_of_elements() != y->get_number_of_elements() )
            throw std::runtime_error("Gadgetron::axpy(): Array sizes mismatch");

        typedef typename stdType<T>::Type stdT;
        arma::Col<stdT> xM = as_arma_col(x);
        arma::Col<stdT> yM = as_arma_col(y);
        stdT a2 = *((stdT*)(&a));
        yM += (a2*xM);
    }

    #ifdef USE_MKL

    template<> float nrm1( hoNDArray<float> *x )
    {
        if ( x == NULL ) return 0;
        MKL_INT N = x->get_number_of_elements();
        MKL_INT incx = 1;
        return(sasum(&N, x->begin(), &incx));
    }

    template<> double nrm1( hoNDArray<double> *x )
    {
        if ( x == NULL ) return 0;
        MKL_INT N = x->get_number_of_elements();
        MKL_INT incx = 1;
        return(dasum(&N, x->begin(), &incx));
    }

    // BLAS dotc and dotu
    // res = conj(x) dot y
    GT_Complex8 dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y)
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

    GT_Complex16 dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y)
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

    // res = x dot y
    GT_Complex8 dotu(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y)
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

    GT_Complex16 dotu(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y)
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

    // other variants for axpy
    // r = a*x+y
    bool axpy(float a, const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
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

    bool axpy(double a, const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
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

    bool axpy(const GT_Complex8& a, const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
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

    bool axpy(const GT_Complex16& a, const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
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

    // vector-scalar product
    // r = a*x
    bool scal(float a, hoNDArray<float>& x)
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

    bool scal(double a, hoNDArray<double>& x)
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

    bool scal(float a, hoNDArray<GT_Complex8>& x)
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

    bool scal(double a, hoNDArray<GT_Complex16>& x)
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

    bool scal(GT_Complex8 a, hoNDArray<GT_Complex8>& x)
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

    bool scal(GT_Complex16 a, hoNDArray<GT_Complex16>& x)
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

    bool scal(float a, float*x, long long N)
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

    bool scal(double a, double*x, long long N)
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

    bool scal(float a, GT_Complex8*x, long long N)
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

    bool scal(double a, GT_Complex16*x, long long N)
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

    bool scal(GT_Complex8 a, GT_Complex8*x, long long N)
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

    bool scal(GT_Complex16 a, GT_Complex16*x, long long N)
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

    // sort the vector
    // isascending: true for ascending and false for descending
    bool sort(const hoNDArray<float>& x, hoNDArray<float>& r, bool isascending)
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

    bool sort(const hoNDArray<double>& x, hoNDArray<double>& r, bool isascending)
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

    //
    // Instantiation
    //

    template EXPORTCPUCOREMATH float dot<float>( hoNDArray<float>*, hoNDArray<float>*, bool );
    template EXPORTCPUCOREMATH float asum<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH float nrm2<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH unsigned long long amin<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH unsigned long long amax<float>( hoNDArray<float>* );
    template EXPORTCPUCOREMATH void axpy<float>( float, hoNDArray<float>*, hoNDArray<float>* );

    template EXPORTCPUCOREMATH double dot<double>( hoNDArray<double>*, hoNDArray<double>*, bool );
    template EXPORTCPUCOREMATH double asum<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH double nrm2<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH unsigned long long amin<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH unsigned long long amax<double>( hoNDArray<double>* );
    template EXPORTCPUCOREMATH void axpy<double>( double, hoNDArray<double>*, hoNDArray<double>* );

    template EXPORTCPUCOREMATH std::complex<float> dot< std::complex<float> >( hoNDArray< std::complex<float> >*, hoNDArray< std::complex<float> >*, bool );
    template EXPORTCPUCOREMATH float asum<float>( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH float nrm2< std::complex<float> >( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH unsigned long long amin<float>( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH unsigned long long amax<float>( hoNDArray< std::complex<float> >* );
    template EXPORTCPUCOREMATH void axpy< std::complex<float> >( std::complex<float> , hoNDArray< std::complex<float> >*, hoNDArray< std::complex<float> >* );

    template EXPORTCPUCOREMATH std::complex<double> dot< std::complex<double> >( hoNDArray< std::complex<double> >*, hoNDArray< std::complex<double> >*, bool );
    template EXPORTCPUCOREMATH double asum<double>( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH double nrm2< std::complex<double> >( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH unsigned long long amin<double>( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH unsigned long long amax<double>( hoNDArray< std::complex<double> >* );
    template EXPORTCPUCOREMATH void axpy< std::complex<double> >( std::complex<double> , hoNDArray< std::complex<double> >*, hoNDArray< std::complex<double> >* );

    template EXPORTCPUCOREMATH complext<float> dot< complext<float> >( hoNDArray< complext<float> >*, hoNDArray< complext<float> >*, bool );
    template EXPORTCPUCOREMATH float asum<float>( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH float nrm2< complext<float> >( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH unsigned long long amin<float>( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH unsigned long long amax<float>( hoNDArray< complext<float> >* );
    template EXPORTCPUCOREMATH void axpy< complext<float> >( complext<float> , hoNDArray< complext<float> >*, hoNDArray< complext<float> >* );

    template EXPORTCPUCOREMATH complext<double> dot< complext<double> >( hoNDArray< complext<double> >*, hoNDArray< complext<double> >*, bool );
    template EXPORTCPUCOREMATH double asum<double>( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH double nrm2< complext<double> >( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH unsigned long long amin<double>( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH unsigned long long amax<double>( hoNDArray< complext<double> >* );
    template EXPORTCPUCOREMATH void axpy< complext<double> >( complext<double> , hoNDArray< complext<double> >*, hoNDArray< complext<double> >* );
}
