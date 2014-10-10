#include "hoNDArray_math_util.h"

namespace Gadgetron
{

// #ifdef USE_MKL

    // // ----------------------------------------------------------------------------------------
    // // float
    // // ----------------------------------------------------------------------------------------

    // EXPORTCPUCOREMATH bool add(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vsAdd(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool add(size_t N, const float* x, const float* y, float* r)
    // {
        // vsAdd(N, x, y, r);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool subtract(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vsSub(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool multiply(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vsMul(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool divide(const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vsDiv(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool absolute(const hoNDArray<float>& x, hoNDArray<float>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vsAbs(x.get_number_of_elements(), x.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool argument(const hoNDArray<float>& x, hoNDArray<float>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // memset(r.begin(), 0, r.get_number_of_bytes());

        // return true;
    // }

    // EXPORTCPUCOREMATH bool sqrt(const hoNDArray<float>& x, hoNDArray<float>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vsSqrt(x.get_number_of_elements(), x.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<float>& x, float& r, size_t& ind)
    // {
        // try
        // {
            // /*MKL_INT n = x.get_number_of_elements();
            // MKL_INT incx = 1;
            // ind = (size_t)(cblas_isamin(n, x.begin(), incx));
            // r = x.at(ind);*/

           // size_t n = x.get_number_of_elements();
           // GADGET_CHECK_RETURN_FALSE(n>0);
           // if ( n == 1 )
           // {
               // r = x(0);
               // ind = 0;
               // return true;
           // }

           // const float* pX = x.begin();
           // r = pX[0];
           // ind = 0;

           // for ( size_t ii=1; ii<n; ii++ )
           // {
               // if ( GT_ABS(pX[ii]) < r )
               // {
                   // r = pX[ii];
                   // ind = ii;
               // }
           // }
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<float>& x, float& r, size_t& ind)
    // {
        // try
        // {
            // /*MKL_INT n = x.get_number_of_elements();
            // MKL_INT incx = 1;
            // ind = (size_t)(cblas_isamax(n, x.begin(), incx));
            // r = x.at(ind);*/

           // size_t n = x.get_number_of_elements();
           // GADGET_CHECK_RETURN_FALSE(n>0);
           // if ( n == 1 )
           // {
               // r = x(0);
               // ind = 0;
               // return true;
           // }

           // const float* pX = x.begin();
           // r = pX[0];
           // ind = 0;

           // for ( size_t ii=1; ii<n; ii++ )
           // {
               // if ( GT_ABS(pX[ii]) > r )
               // {
                   // r = pX[ii];
                   // ind = ii;
               // }
           // }
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<float>& x)
    // {
        // try
        // {
            // size_t n = x.get_number_of_elements();
            // float* pX = x.begin();

            // long long i;

            // #pragma omp parallel for default(none) private(i) shared(n, pX)
            // for (i=0; i<(long long)n; i++ )
            // {
                // if ( GT_ABS(pX[i]) < FLT_EPSILON )
                // {
                    // pX[i] += GT_SGN(pX[i])*FLT_EPSILON;
                // }
            // }
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool norm2(const hoNDArray<float>& x, float& r)
    // {
        // try
        // {
            // MKL_INT incx = 1;
            // MKL_INT n = x.get_number_of_elements();
            // r = snrm2(&n, x.begin(), &incx);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool norm1(const hoNDArray<float>& x, float& r)
    // {
        // try
        // {
            // MKL_INT incx = 1;
            // MKL_INT n = x.get_number_of_elements();
            // r = sasum(&n, x.begin(), &incx);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool inv(const hoNDArray<float>& x, hoNDArray<float>& r)
    // {
        // try
        // {
            // if ( !r.dimensions_equal(&x) )
            // {
                // r = x;
            // }

            // long long n = x.get_number_of_elements();
            // vsInv(n, x.begin(), r.begin());
            // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<float>& x, hoNDArray<float>& r) ... ");
            // return false;
        // }

        // return true;
    // }

    // // ----------------------------------------------------------------------------------------
    // // double
    // // ----------------------------------------------------------------------------------------

    // EXPORTCPUCOREMATH bool add(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vdAdd(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool add(size_t N, const double* x, const double* y, double* r)
    // {
        // vdAdd(N, x, y, r);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool subtract(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vdSub(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool multiply(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vdMul(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool divide(const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vdDiv(x.get_number_of_elements(), x.begin(), y.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool absolute(const hoNDArray<double>& x, hoNDArray<double>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vdAbs(x.get_number_of_elements(), x.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool argument(const hoNDArray<double>& x, hoNDArray<double>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // memset(r.begin(), 0, r.get_number_of_bytes());

        // return true;
    // }

    // EXPORTCPUCOREMATH bool sqrt(const hoNDArray<double>& x, hoNDArray<double>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vdSqrt(x.get_number_of_elements(), x.begin(), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<double>& x, double& r, size_t& ind)
    // {
        // try
        // {
            // //MKL_INT n = x.get_number_of_elements();
            // //MKL_INT incx = 1;
            // //ind = (size_t)(idamin(&n, x.begin(), &incx));
            // //r = x.at(ind);

           // size_t n = x.get_number_of_elements();
           // GADGET_CHECK_RETURN_FALSE(n>0);
           // if ( n == 1 )
           // {
               // r = x(0);
               // ind = 0;
               // return true;
           // }

           // const double* pX = x.begin();
           // r = pX[0];
           // ind = 0;

           // for ( size_t ii=1; ii<n; ii++ )
           // {
               // if ( GT_ABS(pX[ii]) < r )
               // {
                   // r = pX[ii];
                   // ind = ii;
               // }
           // }
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<double>& x, double& r, size_t& ind)
    // {
        // try
        // {
            // //MKL_INT n = x.get_number_of_elements();
            // //MKL_INT incx = 1;
            // //ind = (size_t)(idamax(&n, x.begin(), &incx));
            // //r = x.at(ind);

           // size_t n = x.get_number_of_elements();
           // GADGET_CHECK_RETURN_FALSE(n>0);
           // if ( n == 1 )
           // {
               // r = x(0);
               // ind = 0;
               // return true;
           // }

           // const double* pX = x.begin();
           // r = pX[0];
           // ind = 0;

           // for ( size_t ii=1; ii<n; ii++ )
           // {
               // if ( GT_ABS(pX[ii]) > r )
               // {
                   // r = pX[ii];
                   // ind = ii;
               // }
           // }
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<double>& x)
    // {
        // try
        // {
            // size_t n = x.get_number_of_elements();
            // double* pX = x.begin();

            // long long i;

            // #pragma omp parallel for default(none) private(i) shared(n, pX)
            // for (i=0; i<(long long)n; i++ )
            // {
                // if ( GT_ABS(pX[i]) < DBL_EPSILON )
                // {
                    // pX[i] += GT_SGN(pX[i])*DBL_EPSILON;
                // }
            // }
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool norm2(const hoNDArray<double>& x, double& r)
    // {
        // try
        // {
            // MKL_INT incx = 1;
            // MKL_INT n = x.get_number_of_elements();
            // r = dnrm2(&n, x.begin(), &incx);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool norm1(const hoNDArray<double>& x, double& r)
    // {
        // try
        // {
            // MKL_INT incx = 1;
            // MKL_INT n = x.get_number_of_elements();
            // r = dasum(&n, x.begin(), &incx);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool inv(const hoNDArray<double>& x, hoNDArray<double>& r)
    // {
        // try
        // {
            // if ( !r.dimensions_equal(&x) )
            // {
                // r = x;
            // }

            // long long n = x.get_number_of_elements();
            // vdInv(n, x.begin(), r.begin());
            // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<double>& x, hoNDArray<double>& r) ... ");
            // return false;
        // }

        // return true;
    // }

    // // ----------------------------------------------------------------------------------------
    // // GT_Complex8
    // // ----------------------------------------------------------------------------------------

    // EXPORTCPUCOREMATH bool add(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vcAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool add(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && z!=NULL);
        // vcAdd(N, reinterpret_cast<const MKL_Complex8*>(x), reinterpret_cast<const MKL_Complex8*>(y), reinterpret_cast<MKL_Complex8*>(r));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool subtract(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vcSub(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool subtract(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && z!=NULL);
        // vcSub(N, reinterpret_cast<const MKL_Complex8*>(x), reinterpret_cast<const MKL_Complex8*>(y), reinterpret_cast<MKL_Complex8*>(r));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool multiply(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vcMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool multiply(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8* r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && z!=NULL);
        // vcMul(N, reinterpret_cast<const MKL_Complex8*>(x), reinterpret_cast<const MKL_Complex8*>(y), reinterpret_cast<MKL_Complex8*>(r));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool divide(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vcDiv(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool sqrt(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r.create(x.get_dimensions());
        // }

        // vcSqrt(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind)
    // {
        // try
        // {
            // MKL_INT n = x.get_number_of_elements();
            // MKL_INT incx = 1;
            // ind = (size_t)(icamin(&n, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx));
            // r = x.at(ind);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<GT_Complex8>& x, GT_Complex8& r, size_t& ind)
    // {
        // try
        // {
            // MKL_INT n = x.get_number_of_elements();
            // MKL_INT incx = 1;
            // ind = (size_t)(icamax(&n, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx));
            // r = x.at(ind);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vcMulByConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<const MKL_Complex8*>(y.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool conjugate(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r.create(x.get_dimensions());
        // }

        // vcConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<GT_Complex8>& x)
    // {
        // try
        // {
            // size_t n = x.get_number_of_elements();
            // GT_Complex8* pX = x.begin();

            // long long i;

            // #pragma omp parallel for default(none) private(i) shared(n, pX)
            // for (i=0; i<(long long)n; i++ )
            // {
                // if ( std::abs(pX[i]) < FLT_EPSILON )
                // {
                    // pX[i] += FLT_EPSILON;
                // }
            // }
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool norm2(const hoNDArray<GT_Complex8>& x, float& r)
    // {
        // try
        // {
            // MKL_INT incx = 1;
            // MKL_INT n = x.get_number_of_elements();
            // r = scnrm2(&n, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool norm1(const hoNDArray<GT_Complex8>& x, float& r)
    // {
        // try
        // {
            // hoNDArray<float> a;
            // GADGET_CHECK_RETURN_FALSE(absolute(x, a));
            // GADGET_CHECK_RETURN_FALSE(norm1(a, r));
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, GT_Complex8& r)
    // {
        // try
        // {
            // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            // MKL_INT N = x.get_number_of_elements();
            // MKL_INT incx(1), incy(1);
            // cdotc(reinterpret_cast<MKL_Complex8*>(&r), &N, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx, 
                    // reinterpret_cast<const MKL_Complex8*>(y.begin()), &incy);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r.create(x.get_dimensions());
        // }

        // vcAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool argument(const hoNDArray<GT_Complex8>& x, hoNDArray<float>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r.create(x.get_dimensions());
        // }

        // vcArg(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool inv(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r)
    // {
        // try
        // {
            // if ( !r.dimensions_equal(&x) )
            // {
                // r = x;
            // }

            // const GT_Complex8* pX = x.begin();
            // GT_Complex8* pR = r.begin();

            // GT_Complex8 v(1.0);
            // long long n = x.get_number_of_elements();
            // long long ii;

            // #pragma omp parallel for default(none) private(ii) shared(n, pX, pR, v)
            // for ( ii=0; ii<n; ii++ )
            // {
                // pR[ii] = v/pX[ii];
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<GT_Complex8>& x, hoNDArray<GT_Complex8>& r) ... ");
            // return false;
        // }

        // return true;
    // }

    // // ----------------------------------------------------------------------------------------
    // // GT_Complex16
    // // ----------------------------------------------------------------------------------------

    // EXPORTCPUCOREMATH bool add(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vzAdd(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool add(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && r!=NULL);
        // vzAdd(N, reinterpret_cast<const MKL_Complex16*>(x), reinterpret_cast<const MKL_Complex16*>(y), reinterpret_cast<MKL_Complex16*>(r));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool subtract(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vzSub(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool subtract(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && r!=NULL);
        // vzSub(N, reinterpret_cast<const MKL_Complex16*>(x), reinterpret_cast<const MKL_Complex16*>(y), reinterpret_cast<MKL_Complex16*>(r));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool multiply(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vzMul(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool multiply(size_t N, const GT_Complex16* x, const GT_Complex16* y, GT_Complex16* r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x!=NULL && y!=NULL && r!=NULL);
        // vzMul(N, reinterpret_cast<const MKL_Complex16*>(x), reinterpret_cast<const MKL_Complex16*>(y), reinterpret_cast<MKL_Complex16*>(r));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);
        // return true;
    // }

    // EXPORTCPUCOREMATH bool divide(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vzDiv(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r.create(x.get_dimensions());
        // }

        // vzAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool absolute(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r.create(x.get_dimensions());
        // }

        // hoNDArray<double> rTmp;
        // rTmp.create(x.get_dimensions());

        // vzAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), rTmp.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // //GADGET_CHECK_RETURN_FALSE(r.copyFrom(rTmp));
        // r.copyFrom(rTmp);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool sqrt(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r.create(x.get_dimensions());
        // }

        // vzSqrt(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool minAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind)
    // {
        // try
        // {
            // MKL_INT n = x.get_number_of_elements();
            // MKL_INT incx = 1;
            // ind = (size_t)(izamin(&n, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx));
            // r = x.at(ind);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool maxAbsolute(const hoNDArray<GT_Complex16>& x, GT_Complex16& r, size_t& ind)
    // {
        // try
        // {
            // MKL_INT n = x.get_number_of_elements();
            // MKL_INT incx = 1;
            // ind = (size_t)(izamax(&n, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx));
            // r = x.at(ind);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool multiplyConj(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    // {
        // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r = x;
        // }

        // vzMulByConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<const MKL_Complex16*>(y.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // bool argument(const hoNDArray<GT_Complex16>& x, hoNDArray<double>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r.create(x.get_dimensions());
        // }

        // vzArg(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), r.begin());
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool conjugate(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    // {
        // if ( r.get_number_of_elements()!=x.get_number_of_elements())
        // {
            // r.create(x.get_dimensions());
        // }

        // vzConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
        // GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

        // return true;
    // }

    // EXPORTCPUCOREMATH bool addEpsilon(hoNDArray<GT_Complex16>& x)
    // {
        // try
        // {
            // size_t n = x.get_number_of_elements();
            // GT_Complex16* pX = x.begin();

            // long long i;

            // #pragma omp parallel for default(none) private(i) shared(n, pX)
            // for (i=0; i<(long long)n; i++ )
            // {
                // if ( std::abs(pX[i]) < DBL_EPSILON )
                // {
                    // pX[i] += DBL_EPSILON;
                // }
            // }
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool norm2(const hoNDArray<GT_Complex16>& x, double& r)
    // {
        // try
        // {
            // MKL_INT incx = 1;
            // MKL_INT n = x.get_number_of_elements();
            // r = dznrm2(&n, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool norm1(const hoNDArray<GT_Complex16>& x, double& r)
    // {
        // try
        // {
            // hoNDArray<double> a;
            // GADGET_CHECK_RETURN_FALSE(absolute(x, a));
            // GADGET_CHECK_RETURN_FALSE(norm1(a, r));
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, GT_Complex16& r)
    // {
        // try
        // {
            // GADGET_DEBUG_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            // MKL_INT N = x.get_number_of_elements();
            // MKL_INT incx(1), incy(1);
            // zdotc(reinterpret_cast<MKL_Complex16*>(&r), &N, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx, 
                    // reinterpret_cast<const MKL_Complex16*>(y.begin()), &incy);
        // }
        // catch(...)
        // {
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool inv(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r)
    // {
        // try
        // {
            // if ( !r.dimensions_equal(&x) )
            // {
                // r = x;
            // }

            // const GT_Complex16* pX = x.begin();
            // GT_Complex16* pR = r.begin();

            // GT_Complex16 v(1.0);
            // long long n = x.get_number_of_elements();
            // long long ii;

            // #pragma omp parallel for default(none) private(ii) shared(n, pX, pR, v)
            // for ( ii=0; ii<n; ii++ )
            // {
                // pR[ii] = v/pX[ii];
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in inv(const hoNDArray<GT_Complex16>& x, hoNDArray<GT_Complex16>& r) ... ");
            // return false;
        // }

        // return true;
    // }

    // // ----------------------------------------------------------------------------------------
    // // other functions
    // // ----------------------------------------------------------------------------------------

    // EXPORTCPUCOREMATH GT_Complex8 dotc(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y)
    // {
        // if ( x.get_number_of_elements() != y.get_number_of_elements() )
        // {
            // GADGET_ERROR_MSG("dotc(x, y), inputs have differnet length ...");
            // return 0.0;
        // }

        // MKL_INT N = x.get_number_of_elements();
        // MKL_INT incx(1), incy(1);
        // GT_Complex8 r;
        // cdotc(reinterpret_cast<MKL_Complex8*>(&r), &N, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx, reinterpret_cast<const MKL_Complex8*>(y.begin()), &incy);
        // return r;
    // }

    // EXPORTCPUCOREMATH GT_Complex16 dotc(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y)
    // {
        // if ( x.get_number_of_elements() != y.get_number_of_elements() )
        // {
            // GADGET_ERROR_MSG("dotc(x, y), inputs have differnet length ...");
            // return 0;
        // }

        // MKL_INT N = x.get_number_of_elements();
        // MKL_INT incx(1), incy(1);
        // GT_Complex16 r;
        // zdotc(reinterpret_cast<MKL_Complex16*>(&r), &N, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx, reinterpret_cast<const MKL_Complex16*>(y.begin()), &incy);
        // return r;
    // }

    // EXPORTCPUCOREMATH GT_Complex8 dotu(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y)
    // {
        // if ( x.get_number_of_elements() != y.get_number_of_elements() )
        // {
            // GADGET_ERROR_MSG("dotu(x, y), inputs have differnet length ...");
            // return 0;
        // }

        // MKL_INT N = x.get_number_of_elements();
        // MKL_INT incx(1), incy(1);
        // GT_Complex8 r;
        // cdotu(reinterpret_cast<MKL_Complex8*>(&r), &N, reinterpret_cast<const MKL_Complex8*>(x.begin()), &incx, reinterpret_cast<const MKL_Complex8*>(y.begin()), &incy);
        // return r;
    // }

    // EXPORTCPUCOREMATH GT_Complex16 dotu(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y)
    // {
        // if ( x.get_number_of_elements() != y.get_number_of_elements() )
        // {
            // GADGET_ERROR_MSG("dotu(x, y), inputs have differnet length ...");
            // return 0;
        // }

        // MKL_INT N = x.get_number_of_elements();
        // MKL_INT incx(1), incy(1);
        // GT_Complex16 r;
        // zdotu(reinterpret_cast<MKL_Complex16*>(&r), &N, reinterpret_cast<const MKL_Complex16*>(x.begin()), &incx, reinterpret_cast<const MKL_Complex16*>(y.begin()), &incy);
        // return r;
    // }

    // EXPORTCPUCOREMATH bool axpy(float a, const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r)
    // {
        // try
        // {
            // GADGET_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            // if ( r.get_number_of_elements() != x.get_number_of_elements() )
            // {
                // r = y;
            // }
            // else
            // {
                // if ( &r != &y )
                // {
                    // memcpy(r.begin(), y.begin(), r.get_number_of_bytes());
                // }
            // }

            // MKL_INT N = (MKL_INT)(x.get_number_of_elements());
            // const MKL_INT incX(1), incY(1);

            // cblas_saxpy (N, a, x.begin(), incX, r.begin(), incY);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in axpy(float a, const hoNDArray<float>& x, const hoNDArray<float>& y, hoNDArray<float>& r) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool axpy(double a, const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r)
    // {
        // try
        // {
            // GADGET_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            // if ( r.get_number_of_elements() != x.get_number_of_elements() )
            // {
                // r = y;
            // }
            // else
            // {
                // if ( &r != &y )
                // {
                    // memcpy(r.begin(), y.begin(), r.get_number_of_bytes());
                // }
            // }

            // MKL_INT N = (MKL_INT)(x.get_number_of_elements());
            // const MKL_INT incX(1), incY(1);

            // cblas_daxpy (N, a, x.begin(), incX, r.begin(), incY);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in axpy(double a, const hoNDArray<double>& x, const hoNDArray<double>& y, hoNDArray<double>& r) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool axpy(const GT_Complex8& a, const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r)
    // {
        // try
        // {
            // GADGET_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            // if ( r.get_number_of_elements() != x.get_number_of_elements() )
            // {
                // r = y;
            // }
            // else
            // {
                // if ( &r != &y )
                // {
                    // memcpy(r.begin(), y.begin(), r.get_number_of_bytes());
                // }
            // }

            // MKL_INT N = (MKL_INT)(x.get_number_of_elements());
            // const MKL_INT incX(1), incY(1);

            // cblas_caxpy (N, &a, x.begin(), incX, r.begin(), incY);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in axpy(const GT_Complex8& a, const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& y, hoNDArray<GT_Complex8>& r) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool axpy(const GT_Complex16& a, const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r)
    // {
        // try
        // {
            // GADGET_CHECK_RETURN_FALSE(x.get_number_of_elements()==y.get_number_of_elements());

            // if ( r.get_number_of_elements() != x.get_number_of_elements() )
            // {
                // r = y;
            // }
            // else
            // {
                // if ( &r != &y )
                // {
                    // memcpy(r.begin(), y.begin(), r.get_number_of_bytes());
                // }
            // }

            // MKL_INT N = (MKL_INT)(x.get_number_of_elements());
            // const MKL_INT incX(1), incY(1);

            // cblas_zaxpy (N, &a, x.begin(), incX, r.begin(), incY);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in axpy(const GT_Complex16& a, const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& y, hoNDArray<GT_Complex16>& r) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(float a, hoNDArray<float>& x)
    // {
        // try
        // {
            // cblas_sscal ((MKL_INT)(x.get_number_of_elements()), a, x.begin(), 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(float a, hoNDArray<float>& x) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(double a, hoNDArray<double>& x)
    // {
        // try
        // {
            // cblas_dscal ((MKL_INT)(x.get_number_of_elements()), a, x.begin(), 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(double a, hoNDArray<double>& x) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(float a, hoNDArray<GT_Complex8>& x)
    // {
        // try
        // {
            // GT_Complex8 alpha = GT_Complex8(a);
            // cblas_cscal (x.get_number_of_elements(), &alpha, x.begin(), 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(float a, hoNDArray<GT_Complex8>& x) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(double a, hoNDArray<GT_Complex16>& x)
    // {
        // try
        // {
            // GT_Complex16 alpha = GT_Complex16(a);
            // cblas_zscal (x.get_number_of_elements(), &alpha, x.begin(), 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(double a, hoNDArray<GT_Complex16>& x) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(GT_Complex8 a, hoNDArray<GT_Complex8>& x)
    // {
        // try
        // {
            // cblas_cscal (x.get_number_of_elements(), &a, x.begin(), 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(GT_Complex8 a, hoNDArray<GT_Complex8>& x) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(GT_Complex16 a, hoNDArray<GT_Complex16>& x)
    // {
        // try
        // {
            // cblas_zscal (x.get_number_of_elements(), &a, x.begin(), 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(GT_Complex16 a, hoNDArray<GT_Complex16>& x) ... ");
            // return false;
        // }

        // return true;
    // }

    // // -----------------------

    // EXPORTCPUCOREMATH bool scal(float a, float*x, long long N)
    // {
        // try
        // {
            // cblas_sscal ((MKL_INT)(N), a, x, 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(float a, float*x, long long N) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(double a, double*x, long long N)
    // {
        // try
        // {
            // cblas_dscal ((MKL_INT)(N), a, x, 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(double a, double*x, long long N) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(float a, GT_Complex8*x, long long N)
    // {
        // try
        // {
            // GT_Complex8 alpha = GT_Complex8(a);
            // cblas_cscal (N, &alpha, x, 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(float a, GT_Complex8*x, long long N) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(double a, GT_Complex16*x, long long N)
    // {
        // try
        // {
            // GT_Complex16 alpha = GT_Complex16(a);
            // cblas_zscal (N, &alpha, x, 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(double a, GT_Complex16*x, long long N) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(GT_Complex8 a, GT_Complex8*x, long long N)
    // {
        // try
        // {
            // cblas_cscal (N, &a, x, 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(GT_Complex8 a, GT_Complex8*x, long long N) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool scal(GT_Complex16 a, GT_Complex16*x, long long N)
    // {
        // try
        // {
            // cblas_zscal (N, &a, x, 1);
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors in scal(GT_Complex16 a, GT_Complex16*x, long long N) ... ");
            // return false;
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool sort(const hoNDArray<float>& x, hoNDArray<float>& r, bool isascending)
    // {
        // if ( &r != &x )
        // {
            // if ( r.get_number_of_elements()!=x.get_number_of_elements())
            // {
                // r = x;
            // }
            // else
            // {
                // memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
            // }
        // }

        // if ( isascending )
        // {
            // GADGET_CHECK_RETURN_FALSE(LAPACKE_slasrt('I', r.get_number_of_elements(), r.begin())==0);
        // }
        // else
        // {
            // GADGET_CHECK_RETURN_FALSE(LAPACKE_slasrt('D', r.get_number_of_elements(), r.begin())==0);
        // }

        // return true;
    // }

    // EXPORTCPUCOREMATH bool sort(const hoNDArray<double>& x, hoNDArray<double>& r, bool isascending)
    // {
        // if ( &r != &x )
        // {
            // if ( r.get_number_of_elements()!=x.get_number_of_elements())
            // {
                // r = x;
            // }
            // else
            // {
                // memcpy(r.begin(), x.begin(), x.get_number_of_bytes());
            // }
        // }

        // if ( isascending )
        // {
            // GADGET_CHECK_RETURN_FALSE(LAPACKE_dlasrt('I', r.get_number_of_elements(), r.begin())==0);
        // }
        // else
        // {
            // GADGET_CHECK_RETURN_FALSE(LAPACKE_dlasrt('D', r.get_number_of_elements(), r.begin())==0);
        // }

        // return true;
    // }

// #endif // USE_MKL

// #ifdef USE_MKL

    // // ----------------------------------------------------------------------------------------
    // // float
    // // ----------------------------------------------------------------------------------------

    // bool conv2(const hoNDArray<float>& x, const hoNDArray<float>& ker, hoNDArray<float>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);

            // size_t num = x.get_number_of_elements()/(RO*E1);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[2];
            // kerShape[0] = kerRO; kerShape[1] = kerE1;

            // MKL_INT xshape[2];
            // xshape[0] = RO; xshape[1] = E1;

            // MKL_INT start[2];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;

            // MKL_INT kerStride[2], xstride[2], zstride[2];
            // kerStride[0] = 1; kerStride[1] = kerRO;
            // xstride[0] = 1; xstride[1] = RO;
            // zstride[0] = 1; zstride[1] = RO;

            // const float* pX = x.begin();
            // const float* pKer = ker.begin();
            // float* pZ = z.begin();

            // if ( num == 1 )
            // {
                // status = vslsConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslsConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslConvDeleteTask(&task);
            // }
            // else
            // {
                // status = vslsConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslsConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslConvDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<float>& x, const hoNDArray<float>& ker, hoNDArray<float>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // bool conv3(const hoNDArray<float>& x, const hoNDArray<float>& ker, hoNDArray<float>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);
            // size_t E2 = x.get_size(2);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);
            // size_t kerE2 = ker.get_size(2);

            // size_t num = x.get_number_of_elements()/(RO*E1*E2);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[3];
            // kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            // MKL_INT xshape[3];
            // xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            // MKL_INT start[3];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;
            // start[2] = kerE2/2;

            // MKL_INT kerStride[3], xstride[3], zstride[3];
            // kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
            // xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
            // zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

            // const float* pX = x.begin();
            // const float* pKer = ker.begin();
            // float* pZ = z.begin();

            // if ( num == 1 )
            // {
                // status = vslsConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslsConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslConvDeleteTask(&task);
            // }
            // else
            // {
                // status = vslsConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslsConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslConvDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<float>& x, const hoNDArray<float>& ker, hoNDArray<float>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // // ----------------------------------------------------------------------------------------
    // // double
    // // ----------------------------------------------------------------------------------------

    // bool conv2(const hoNDArray<double>& x, const hoNDArray<double>& ker, hoNDArray<double>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);

            // size_t num = x.get_number_of_elements()/(RO*E1);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[2];
            // kerShape[0] = kerRO; kerShape[1] = kerE1;

            // MKL_INT xshape[2];
            // xshape[0] = RO; xshape[1] = E1;

            // MKL_INT start[2];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;

            // MKL_INT kerStride[2], xstride[2], zstride[2];
            // kerStride[0] = 1; kerStride[1] = kerRO;
            // xstride[0] = 1; xstride[1] = RO;
            // zstride[0] = 1; zstride[1] = RO;

            // const double* pX = x.begin();
            // const double* pKer = ker.begin();
            // double* pZ = z.begin();

            // if ( num == 1 )
            // {
                // status = vsldConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vsldConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslConvDeleteTask(&task);
            // }
            // else
            // {
                // status = vsldConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vsldConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslConvDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<double>& x, const hoNDArray<double>& ker, hoNDArray<double>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // bool conv3(const hoNDArray<double>& x, const hoNDArray<double>& ker, hoNDArray<double>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);
            // size_t E2 = x.get_size(2);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);
            // size_t kerE2 = ker.get_size(2);

            // size_t num = x.get_number_of_elements()/(RO*E1*E2);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[3];
            // kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            // MKL_INT xshape[3];
            // xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            // MKL_INT start[3];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;
            // start[2] = kerE2/2;

            // MKL_INT kerStride[3], xstride[3], zstride[3];
            // kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
            // xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
            // zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

            // const double* pX = x.begin();
            // const double* pKer = ker.begin();
            // double* pZ = z.begin();

            // if ( num == 1 )
            // {
                // status = vsldConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vsldConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslConvDeleteTask(&task);
            // }
            // else
            // {
                // status = vsldConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vsldConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslConvDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<double>& x, const hoNDArray<double>& ker, hoNDArray<double>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // // ----------------------------------------------------------------------------------------
    // // GT_Complex8
    // // ----------------------------------------------------------------------------------------

    // bool conv2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);

            // size_t num = x.get_number_of_elements()/(RO*E1);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[2];
            // kerShape[0] = kerRO; kerShape[1] = kerE1;

            // MKL_INT xshape[2];
            // xshape[0] = RO; xshape[1] = E1;

            // MKL_INT start[2];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;

            // MKL_INT kerStride[2], xstride[2], zstride[2];
            // kerStride[0] = 1; kerStride[1] = kerRO;
            // xstride[0] = 1; xstride[1] = RO;
            // zstride[0] = 1; zstride[1] = RO;

            // const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
            // const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
            // MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

            // if ( num == 1 )
            // {
                // status = vslcConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetInternalPrecision(task, VSL_CONV_PRECISION_SINGLE);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslcConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslConvDeleteTask(&task);
            // }
            // else
            // {
                // status = vslcConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetInternalPrecision(task, VSL_CONV_PRECISION_SINGLE);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslcConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslConvDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // bool conv3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);
            // size_t E2 = x.get_size(2);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);
            // size_t kerE2 = ker.get_size(2);

            // size_t num = x.get_number_of_elements()/(RO*E1*E2);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[3];
            // kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            // MKL_INT xshape[3];
            // xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            // MKL_INT start[3];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;
            // start[2] = kerE2/2;

            // MKL_INT kerStride[3], xstride[3], zstride[3];
            // kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
            // xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
            // zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

            // const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
            // const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
            // MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

            // if ( num == 1 )
            // {
                // status = vslcConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetInternalPrecision(task, VSL_CONV_PRECISION_SINGLE);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslcConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslConvDeleteTask(&task);
            // }
            // else
            // {
                // status = vslcConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetInternalPrecision(task, VSL_CONV_PRECISION_SINGLE);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslcConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslConvDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // bool corr2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);

            // size_t num = x.get_number_of_elements()/(RO*E1);

            // int status;
            // VSLCorrTaskPtr task;

            // MKL_INT kerShape[2];
            // kerShape[0] = kerRO; kerShape[1] = kerE1;

            // MKL_INT xshape[2];
            // xshape[0] = RO; xshape[1] = E1;

            // MKL_INT start[2];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;

            // MKL_INT decimation[2];
            // decimation[0] = 1;
            // decimation[1] = 1;

            // MKL_INT kerStride[2], xstride[2], zstride[2];
            // kerStride[0] = 1; kerStride[1] = kerRO;
            // xstride[0] = 1; xstride[1] = RO;
            // zstride[0] = 1; zstride[1] = RO;

            // const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
            // const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
            // MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

            // if ( num == 1 )
            // {
                // status = vslcCorrNewTask(&task, VSL_CORR_MODE_AUTO, 2, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetDecimation(task, decimation);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetInternalPrecision(task, VSL_CORR_PRECISION_SINGLE);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslcCorrExec(task, pKer, NULL, pX, NULL, pZ, NULL);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslCorrDeleteTask(&task);
            // }
            // else
            // {
                // status = vslcCorrNewTaskX(&task, VSL_CORR_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetDecimation(task, decimation);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetInternalPrecision(task, VSL_CORR_PRECISION_SINGLE);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslcCorrExecX(task, pX+n*RO*E1, NULL, pZ+n*RO*E1, NULL);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslCorrDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in corr2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // bool corr3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);
            // size_t E2 = x.get_size(2);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);
            // size_t kerE2 = ker.get_size(2);

            // size_t num = x.get_number_of_elements()/(RO*E1*E2);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[3];
            // kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            // MKL_INT xshape[3];
            // xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            // MKL_INT start[3];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;
            // start[2] = kerE2/2;

            // MKL_INT kerStride[3], xstride[3], zstride[3];
            // kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
            // xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
            // zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

            // const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
            // const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
            // MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

            // if ( num == 1 )
            // {
                // status = vslcCorrNewTask(&task, VSL_CORR_MODE_AUTO, 3, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetInternalPrecision(task, VSL_CORR_PRECISION_SINGLE);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslcCorrExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslCorrDeleteTask(&task);
            // }
            // else
            // {
                // status = vslcCorrNewTaskX(&task, VSL_CORR_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetInternalPrecision(task, VSL_CORR_PRECISION_SINGLE);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslcCorrExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslCorrDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in corr3(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // // ----------------------------------------------------------------------------------------
    // // GT_Complex16
    // // ----------------------------------------------------------------------------------------

    // bool conv2(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);

            // size_t num = x.get_number_of_elements()/(RO*E1);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[2];
            // kerShape[0] = kerRO; kerShape[1] = kerE1;

            // MKL_INT xshape[2];
            // xshape[0] = RO; xshape[1] = E1;

            // MKL_INT start[2];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;

            // MKL_INT kerStride[2], xstride[2], zstride[2];
            // kerStride[0] = 1; kerStride[1] = kerRO;
            // xstride[0] = 1; xstride[1] = RO;
            // zstride[0] = 1; zstride[1] = RO;

            // const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
            // const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
            // MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

            // if ( num == 1 )
            // {
                // status = vslzConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslzConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslConvDeleteTask(&task);
            // }
            // else
            // {
                // status = vslzConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslzConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslConvDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in conv2(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // bool conv3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);
            // size_t E2 = x.get_size(2);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);
            // size_t kerE2 = ker.get_size(2);

            // size_t num = x.get_number_of_elements()/(RO*E1*E2);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[3];
            // kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            // MKL_INT xshape[3];
            // xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            // MKL_INT start[3];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;
            // start[2] = kerE2/2;

            // MKL_INT kerStride[3], xstride[3], zstride[3];
            // kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
            // xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
            // zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

            // const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
            // const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
            // MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

            // if ( num == 1 )
            // {
                // status = vslzConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslzConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslConvDeleteTask(&task);
            // }
            // else
            // {
                // status = vslzConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslConvSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslzConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslConvDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in conv3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // bool corr2(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);

            // size_t num = x.get_number_of_elements()/(RO*E1);

            // int status;
            // VSLCorrTaskPtr task;

            // MKL_INT kerShape[2];
            // kerShape[0] = kerRO; kerShape[1] = kerE1;

            // MKL_INT xshape[2];
            // xshape[0] = RO; xshape[1] = E1;

            // MKL_INT start[2];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;

            // MKL_INT kerStride[2], xstride[2], zstride[2];
            // kerStride[0] = 1; kerStride[1] = kerRO;
            // xstride[0] = 1; xstride[1] = RO;
            // zstride[0] = 1; zstride[1] = RO;

            // const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
            // const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
            // MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

            // if ( num == 1 )
            // {
                // status = vslzCorrNewTask(&task, VSL_CORR_MODE_AUTO, 2, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslzCorrExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslCorrDeleteTask(&task);
            // }
            // else
            // {
                // status = vslzCorrNewTaskX(&task, VSL_CORR_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslzCorrExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslCorrDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in corr2(const hoNDArray<GT_Complex8>& x, const hoNDArray<GT_Complex8>& ker, hoNDArray<GT_Complex8>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // bool corr3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z)
    // {
        // try
        // {
            // if ( !z.dimensions_equal(&x) )
            // {
                // z = x;
            // }

            // size_t RO = x.get_size(0);
            // size_t E1 = x.get_size(1);
            // size_t E2 = x.get_size(2);

            // size_t kerRO = ker.get_size(0);
            // size_t kerE1 = ker.get_size(1);
            // size_t kerE2 = ker.get_size(2);

            // size_t num = x.get_number_of_elements()/(RO*E1*E2);

            // int status;
            // VSLConvTaskPtr task;

            // MKL_INT kerShape[3];
            // kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

            // MKL_INT xshape[3];
            // xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

            // MKL_INT start[3];
            // start[0] = kerRO/2;
            // start[1] = kerE1/2;
            // start[2] = kerE2/2;

            // MKL_INT kerStride[3], xstride[3], zstride[3];
            // kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerE1;
            // xstride[0] = 1; xstride[1] = RO; xstride[2] = E1;
            // zstride[0] = 1; zstride[1] = RO; zstride[2] = E1;

            // const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
            // const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
            // MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

            // if ( num == 1 )
            // {
                // status = vslzCorrNewTask(&task, VSL_CORR_MODE_AUTO, 3, kerShape, xshape, xshape);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslzCorrExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                 // vslCorrDeleteTask(&task);
            // }
            // else
            // {
                // status = vslzCorrNewTaskX(&task, VSL_CORR_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // status = vslCorrSetStart(task, start);
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // long long n;

                // #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
                // for ( n=0; n<(long long)num; n++ )
                // {
                    // status = vslzCorrExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
                // }
                // GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                // vslCorrDeleteTask(&task);
            // }
        // }
        // catch(...)
        // {
            // GADGET_ERROR_MSG("Errors happened in corr3(const hoNDArray<GT_Complex16>& x, const hoNDArray<GT_Complex16>& ker, hoNDArray<GT_Complex16>& z) ... ");
            // return false;
        // }

        // return true;
    // }

    // #endif // USE_MKL

    // //
    // // Instantiation
    // //

    // // -----------------------------------------------------------

    // #ifdef USE_MKL

    // template <unsigned int D> 
    // inline bool scal(float a, hoNDImage<float, D>& x)
    // {
        // long long N = (long long)(x.get_number_of_elements());
        // return scal(a, x.begin(), N);
    // }

    // template <unsigned int D> 
    // inline bool scal(double a, hoNDImage<double, D>& x)
    // {
        // long long N = (long long)(x.get_number_of_elements());
        // return scal(a, x.begin(), N);
    // }

    // template <unsigned int D> 
    // inline bool scal(float a, hoNDImage<GT_Complex8, D>& x)
    // {
        // long long N = (long long)(x.get_number_of_elements());
        // return scal(a, x.begin(), N);
    // }

    // template <unsigned int D> 
    // inline bool scal(double a, hoNDImage<GT_Complex16, D>& x)
    // {
        // long long N = (long long)(x.get_number_of_elements());
        // return scal(a, x.begin(), N);
    // }

    // template <unsigned int D> 
    // inline bool scal(GT_Complex8 a, hoNDImage<GT_Complex8, D>& x)
    // {
        // long long N = (long long)(x.get_number_of_elements());
        // return scal(a, x.begin(), N);
    // }

    // template <unsigned int D> 
    // inline bool scal(GT_Complex16 a, hoNDImage<GT_Complex16, D>& x)
    // {
        // long long N = (long long)(x.get_number_of_elements());
        // return scal(a, x.begin(), N);
    // }

    // #endif // USE_MKL

    //
    // Instantiation
    //

    // template EXPORTCPUCOREMATH bool scal(float a, hoNDImage<float, 1>& x);
    // template EXPORTCPUCOREMATH bool scal(double a, hoNDImage<double, 1>& x);
    // template EXPORTCPUCOREMATH bool scal(float a, hoNDImage<GT_Complex8, 1>& x);
    // template EXPORTCPUCOREMATH bool scal(double a, hoNDImage<GT_Complex16, 1>& x);
    // template EXPORTCPUCOREMATH bool scal(GT_Complex8 a, hoNDImage<GT_Complex8, 1>& x);
    // template EXPORTCPUCOREMATH bool scal(GT_Complex16 a, hoNDImage<GT_Complex16, 1>& x);

    // template EXPORTCPUCOREMATH bool scal(float a, hoNDImage<float, 2>& x);
    // template EXPORTCPUCOREMATH bool scal(double a, hoNDImage<double, 2>& x);
    // template EXPORTCPUCOREMATH bool scal(float a, hoNDImage<GT_Complex8, 2>& x);
    // template EXPORTCPUCOREMATH bool scal(double a, hoNDImage<GT_Complex16, 2>& x);
    // template EXPORTCPUCOREMATH bool scal(GT_Complex8 a, hoNDImage<GT_Complex8, 2>& x);
    // template EXPORTCPUCOREMATH bool scal(GT_Complex16 a, hoNDImage<GT_Complex16, 2>& x);

    // template EXPORTCPUCOREMATH bool scal(float a, hoNDImage<float, 3>& x);
    // template EXPORTCPUCOREMATH bool scal(double a, hoNDImage<double, 3>& x);
    // template EXPORTCPUCOREMATH bool scal(float a, hoNDImage<GT_Complex8, 3>& x);
    // template EXPORTCPUCOREMATH bool scal(double a, hoNDImage<GT_Complex16, 3>& x);
    // template EXPORTCPUCOREMATH bool scal(GT_Complex8 a, hoNDImage<GT_Complex8, 3>& x);
    // template EXPORTCPUCOREMATH bool scal(GT_Complex16 a, hoNDImage<GT_Complex16, 3>& x);
}
