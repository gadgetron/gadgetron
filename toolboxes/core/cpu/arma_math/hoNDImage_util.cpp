/** \file   hoNDImage_util.hxx
    \brief  operations on the hoNDImage class.
*/

#include "hoNDImage_util.h"
#include "hoNDBoundaryHandler.h"

namespace Gadgetron
{

template<class T, unsigned int D> 
bool gradient(const hoNDImage<T, D>& x, hoNDImage<T, D> gx[])
{
    try
    {
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            if ( !gx[ii].dimensions_equal(x) )
            {
                gx[ii] = x;
            }
        }

        if ( D == 1 )
        {
            long long sx = (long long)x.get_size(0);
            const T* pX = x.begin();
            T* pGx = gx[0].begin();

            long long x;

            #pragma omp parallel for default(none) private(x) shared(sx, pX, pGx)
            for ( x=1; x<sx-1; x++ )
            {
                pGx[x] = pX[x+1] - pX[x-1];
            }

            pGx[0] = pX[1] - pX[0];
            pGx[sx-1] = pX[sx-1] - pX[sx-2];
        }
        else if ( D == 2 )
        {
            long long sx = (long long)x.get_size(0);
            long long sy = (long long)x.get_size(1);

            const T* pX = x.begin();
            T* pGx = gx[0].begin();
            T* pGy = gx[1].begin();

            long long x, y;

            #pragma omp parallel for default(none) private(x, y) shared(sx, sy, pX, pGx, pGy)
            for ( y=1; y<sy-1; y++ )
            {
                for ( x=1; x<sx-1; x++ )
                {
                    size_t offset = x + y*sx;

                    pGx[offset] = pX[offset+1] - pX[offset-1];
                    pGy[offset] = pX[offset+sx] - pX[offset-sx];
                }
            }

            #pragma omp parallel for default(none) private(x) shared(sx, sy, pX, pGx, pGy)
            for ( x=1; x<sx-1; x++ )
            {
                pGx[x] = pX[x+1] - pX[x-1];

                size_t offset = x + (sy-1)*sx;
                pGx[offset] = pX[offset+1] - pX[offset-1];

                pGy[x] = pX[x+sx] - pX[x];
                pGy[x + (sy-1)*sx] = pX[x + (sy-1)*sx] - pX[x + (sy-2)*sx];
            }

            #pragma omp parallel for default(none) private(y) shared(sx, sy, pX, pGx, pGy)
            for ( y=1; y<sy-1; y++ )
            {
                size_t offset = y*sx;
                pGy[offset] = pX[offset+sx] - pX[offset-sx];

                pGx[offset] = pX[offset+1] - pX[offset];

                offset = sx-1 + y*sx;
                pGy[offset] = pX[offset+sx] - pX[offset-sx];

                pGx[offset] = pX[offset] - pX[offset-1];
            }

            pGx[0] = pX[1]-pX[0];
            pGx[sx-1] = pX[sx-1]-pX[sx-2];
            pGx[(sy-1)*sx] = pX[(sy-1)*sx+1]-pX[(sy-1)*sx];
            pGx[sx*sy-1] = pX[sx*sy-1]-pX[sx*sy-2];

            pGy[0] = pX[sx]-pX[0];
            pGy[sx-1] = pX[2*sx-1]-pX[sx-1];
            pGy[(sy-1)*sx] = pX[(sy-1)*sx] - pX[(sy-2)*sx];
            pGy[sx*sy-1] = pX[sx*sy-1] - pX[sx*sy-1-sx];
        }
        else if ( D == 3 )
        {
            long long sx = (long long)x.get_size(0);
            long long sy = (long long)x.get_size(1);
            long long sz = (long long)x.get_size(2);

            const T* pX = x.begin();
            T* pGx = gx[0].begin();
            T* pGy = gx[1].begin();
            T* pGz = gx[2].begin();

            long long x, y, z;

            #pragma omp parallel default(none) private(x, y, z) shared(sx, sy, sz, pX, pGx, pGy, pGz)
            {
                long long z_positive, z_negative, y_positive, y_negative;
                size_t offset, offset_z_positive, offset_z_negative, offset_y_positive, offset_y_negative;

                #pragma omp for 
                for ( z=0; z<sz; z++ )
                {
                    z_positive = z+1;
                    z_positive = (z_positive==sz) ? sz-1 : z_positive;

                    z_negative = z-1;
                    z_negative = (z_negative==-1) ? 0 : z_negative;

                    for ( y=0; y<sy; y++ )
                    {

                        y_positive = y+1;
                        y_positive = (y_positive==sy) ? sy-1 : y_positive;

                        y_negative = y-1;
                        y_negative = (y_negative==-1) ? 0 : y_negative;

                        offset = y*sx + z*sx*sy;

                        offset_z_positive = y*sx + z_positive*sx*sy;
                        offset_z_negative = y*sx + z_negative*sx*sy;

                        offset_y_positive = y_positive*sx + z*sx*sy;
                        offset_y_negative = y_negative*sx + z*sx*sy;

                        for ( x=1; x<sx-1; x++ )
                        {
                            pGx[offset+x] = pX[offset+x+1] - pX[offset+x-1];
                            pGy[offset+x] = pX[offset_y_positive+x] - pX[offset_y_negative+x];
                            pGz[offset+x] = pX[offset_z_positive+x] - pX[offset_z_negative+x];
                        }

                        // x = 0
                        pGx[offset] = pX[offset+1] - pX[offset];
                        pGy[offset] = pX[offset_y_positive] - pX[offset_y_negative];
                        pGz[offset] = pX[offset_z_positive] - pX[offset_z_negative];

                        // x = sx-1
                        pGx[offset+sx-1] = pX[offset+sx-1] - pX[offset+sx-2];
                        pGy[offset+sx-1] = pX[offset_y_positive+sx-1] - pX[offset_y_negative+sx-1];
                        pGz[offset+sx-1] = pX[offset_z_positive+sx-1] - pX[offset_z_negative+sx-1];
                    }
                }
            }
        }
        else
        {
            size_t N = x.get_number_of_elements();

            long long n;

            std::vector<size_t> dim(D);
            x.get_dimensions(dim);

            #pragma omp parallel default(none) private(n) shared(N, dim, x, gx)
            {
                size_t ind[D];
                size_t ind_positive[D];
                size_t ind_negative[D];
                bool inside = true;
                unsigned int ii;

                #pragma omp for 
                for ( n=0; n<(long long)N; n++ )
                {
                    x.calculate_index(n, ind);

                    inside = true;
                    for ( ii=0; ii<D; ii++ )
                    {
                        if ( ind[ii]==0 || ind[ii]==dim[ii]-1 )
                        {
                            inside = false;
                            break;
                        }
                    }

                    if ( inside )
                    {
                        for ( ii=0; ii<D; ii++ )
                        {
                            memcpy(ind_positive, ind, sizeof(size_t)*D);
                            memcpy(ind_negative, ind, sizeof(size_t)*D);

                            ind_positive[ii] = ind[ii] + 1;
                            ind_negative[ii] = ind[ii] - 1;

                            gx[ii](n) = x(ind_positive) - x(ind_negative);
                        }
                    }
                    else
                    {
                        for ( ii=0; ii<D; ii++ )
                        {
                            memcpy(ind_positive, ind, sizeof(size_t)*D);
                            memcpy(ind_negative, ind, sizeof(size_t)*D);

                            ind_positive[ii] = ind[ii] + 1;
                            ind_positive[ii] = (ind_positive[ii]==dim[ii]) ? dim[ii]-1 : ind_positive[ii];

                            ind_negative[ii] = ind[ii] - 1;
                            ind_negative[ii] = (ind_negative[ii]==-1) ? 0 : ind_negative[ii];

                            gx[ii](n) = x(ind_positive) - x(ind_negative);
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in gradient(const hoNDImage<T, D>& x, hoNDImage<T, D> gx[D]) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gaussianKernel(T sigma, double kerWidthInUnitOfSigma, double deltaKer, hoNDArray<T>& ker)
{
    try
    {
        long long N  =  (long long)(2*std::ceil(kerWidthInUnitOfSigma*sigma/deltaKer) + 1);

        ker.create(N);

        T kerSum = 0;

        T D = (deltaKer*deltaKer)/(2*sigma*sigma);

        long long ii;
        for ( ii=-N/2; ii<=N/2; ii++ )
        {
            ker(ii+N/2) = exp( -(ii*ii*D) );
            kerSum += ker(ii+N/2);
        }

        T GNorm = 1/std::sqrt(2*3.141592653579*sigma*sigma);
        GNorm /= kerSum;

        Gadgetron::scal(GNorm, ker);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in gaussianKernel(T sigma, double kerWidthInUnitOfSigma, double deltaKer, hoNDArray<T>& ker) ... ");
        return false;
    }
    return true;
}

template <typename T> 
struct gaussianPara
{
    T alpha, ea, e2a, k; 

    gaussianPara(T sigma)
    {
        if ( sigma > 0 )
        {
            alpha = T(1.4105)/sigma;
        }
        else
        {
            alpha = T(1.4105)/T(1e-5);
        }
        ea = static_cast<T>(exp( (double)(-alpha) ));
        e2a = static_cast<T>(exp( (double)(-2*alpha) ));
        k = static_cast<T>( (1.0 - 2.0*ea + e2a) / (1.0 + 2.0*alpha*ea - e2a) );
    }
};

#ifdef USE_MKL

// ----------------------------------------------------------------------------------------
// float
// ----------------------------------------------------------------------------------------

template <unsigned int D> bool add(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& r)
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

template <unsigned int D> bool subtract(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& r)
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

template <unsigned int D> bool multiply(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& r)
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

template <unsigned int D> bool divide(const hoNDImage<float, D>& x, const hoNDImage<float, D>& y, hoNDImage<float, D>& r)
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

template <unsigned int D> bool absolute(const hoNDImage<float, D>& x, hoNDImage<float, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r = x;
    }

    vsAbs(x.get_number_of_elements(), x.begin(), r.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool argument(const hoNDImage<float, D>& x, hoNDImage<float, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r = x;
    }

    memset(r.begin(), 0, r.get_number_of_bytes());

    return true;
}

template <unsigned int D> bool sqrt(const hoNDImage<float, D>& x, hoNDImage<float, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r = x;
    }

    vsSqrt(x.get_number_of_elements(), x.begin(), r.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool minAbsolute(const hoNDImage<float, D>& x, float& r, size_t& ind)
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

template <unsigned int D> bool maxAbsolute(const hoNDImage<float, D>& x, float& r, size_t& ind)
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

template <unsigned int D> bool addEpsilon(hoNDImage<float, D>& x)
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

template <unsigned int D> bool norm2(const hoNDImage<float, D>& x, float& r)
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

template <unsigned int D> bool norm1(const hoNDImage<float, D>& x, float& r)
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

template <unsigned int D> bool conv2(const hoNDImage<float, D>& x, const hoNDImage<float, D>& ker, hoNDImage<float, D>& z)
{
    try
    {
        if ( !z.dimensions_equal(&x) )
        {
            z = x;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);

        size_t num = x.get_number_of_elements()/(RO*E1);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[2];
        kerShape[0] = kerRO; kerShape[1] = kerE1;

        MKL_INT xshape[2];
        xshape[0] = RO; xshape[1] = E1;

        MKL_INT start[2];
        start[0] = kerRO/2;
        start[1] = kerE1/2;

        MKL_INT kerStride[2], xstride[2], zstride[2];
        kerStride[0] = 1; kerStride[1] = kerRO;
        xstride[0] = 1; xstride[1] = RO;
        zstride[0] = 1; zstride[1] = RO;

        const float* pX = x.begin();
        const float* pKer = ker.begin();
        float* pZ = z.begin();

        if ( num == 1 )
        {
            status = vslsConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslsConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
        }
        else
        {
            status = vslsConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslsConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslConvDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in conv2(const hoNDImage<float, D>& x, const hoNDImage<float, D>& ker, hoNDImage<float, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool conv3(const hoNDImage<float, D>& x, const hoNDImage<float, D>& ker, hoNDImage<float, D>& z)
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

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);
        size_t kerE2 = ker.get_size(2);

        size_t num = x.get_number_of_elements()/(RO*E1*E2);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[3];
        kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

        MKL_INT xshape[3];
        xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

        MKL_INT start[3];
        start[0] = kerRO/2;
        start[1] = kerE1/2;
        start[2] = kerE2/2;

        MKL_INT kerStride[3], xstride[3], zstride[3];
        kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
        xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
        zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

        const float* pX = x.begin();
        const float* pKer = ker.begin();
        float* pZ = z.begin();

        if ( num == 1 )
        {
            status = vslsConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslsConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
        }
        else
        {
            status = vslsConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslsConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslConvDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in conv3(const hoNDImage<float, D>& x, const hoNDImage<float, D>& ker, hoNDImage<float, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool inv(const hoNDImage<float, D>& x, hoNDImage<float, D>& r)
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
        GADGET_ERROR_MSG("Errors happened in inv(const hoNDImage<float, D>& x, hoNDImage<float, D>& r) ... ");
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------------------
// double
// ----------------------------------------------------------------------------------------

template <unsigned int D> bool add(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& r)
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

template <unsigned int D> bool subtract(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& r)
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

template <unsigned int D> bool multiply(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& r)
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

template <unsigned int D> bool divide(const hoNDImage<double, D>& x, const hoNDImage<double, D>& y, hoNDImage<double, D>& r)
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

template <unsigned int D> bool absolute(const hoNDImage<double, D>& x, hoNDImage<double, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r = x;
    }

    vdAbs(x.get_number_of_elements(), x.begin(), r.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool argument(const hoNDImage<double, D>& x, hoNDImage<double, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r = x;
    }

    memset(r.begin(), 0, r.get_number_of_bytes());

    return true;
}

template <unsigned int D> bool sqrt(const hoNDImage<double, D>& x, hoNDImage<double, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r = x;
    }

    vdSqrt(x.get_number_of_elements(), x.begin(), r.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool minAbsolute(const hoNDImage<double, D>& x, double& r, size_t& ind)
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

template <unsigned int D> bool maxAbsolute(const hoNDImage<double, D>& x, double& r, size_t& ind)
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

template <unsigned int D> bool addEpsilon(hoNDImage<double, D>& x)
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

template <unsigned int D> bool norm2(const hoNDImage<double, D>& x, double& r)
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

template <unsigned int D> bool norm1(const hoNDImage<double, D>& x, double& r)
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

template <unsigned int D> bool conv2(const hoNDImage<double, D>& x, const hoNDImage<double, D>& ker, hoNDImage<double, D>& z)
{
    try
    {
        if ( !z.dimensions_equal(&x) )
        {
            z = x;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);

        size_t num = x.get_number_of_elements()/(RO*E1);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[2];
        kerShape[0] = kerRO; kerShape[1] = kerE1;

        MKL_INT xshape[2];
        xshape[0] = RO; xshape[1] = E1;

        MKL_INT start[2];
        start[0] = kerRO/2;
        start[1] = kerE1/2;

        MKL_INT kerStride[2], xstride[2], zstride[2];
        kerStride[0] = 1; kerStride[1] = kerRO;
        xstride[0] = 1; xstride[1] = RO;
        zstride[0] = 1; zstride[1] = RO;

        const double* pX = x.begin();
        const double* pKer = ker.begin();
        double* pZ = z.begin();

        if ( num == 1 )
        {
            status = vsldConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vsldConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
        }
        else
        {
            status = vsldConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vsldConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslConvDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in conv2(const hoNDImage<double, D>& x, const hoNDImage<double, D>& ker, hoNDImage<double, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool conv3(const hoNDImage<double, D>& x, const hoNDImage<double, D>& ker, hoNDImage<double, D>& z)
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

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);
        size_t kerE2 = ker.get_size(2);

        size_t num = x.get_number_of_elements()/(RO*E1*E2);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[3];
        kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

        MKL_INT xshape[3];
        xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

        MKL_INT start[3];
        start[0] = kerRO/2;
        start[1] = kerE1/2;
        start[2] = kerE2/2;

        MKL_INT kerStride[3], xstride[3], zstride[3];
        kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
        xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
        zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

        const double* pX = x.begin();
        const double* pKer = ker.begin();
        double* pZ = z.begin();

        if ( num == 1 )
        {
            status = vsldConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vsldConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
        }
        else
        {
            status = vsldConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vsldConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslConvDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in conv3(const hoNDImage<double, D>& x, const hoNDImage<double, D>& ker, hoNDImage<double, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool inv(const hoNDImage<double, D>& x, hoNDImage<double, D>& r)
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
        GADGET_ERROR_MSG("Errors happened in inv(const hoNDImage<double, D>& x, hoNDImage<double, D>& r) ... ");
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------------------
// GT_Complex8
// ----------------------------------------------------------------------------------------

template <unsigned int D> bool add(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r)
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

template <unsigned int D> bool subtract(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r)
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

template <unsigned int D> bool multiply(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r)
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

template <unsigned int D> bool divide(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r)
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

template <unsigned int D> bool absolute(const hoNDImage<GT_Complex8, D>& x, hoNDImage<float, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    vcAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), r.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool absolute(const hoNDImage<GT_Complex8, D>& x, hoNDImage<GT_Complex8, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    hoNDImage<float, D> rTmp;
    rTmp.create(x.get_dimensions());

    vcAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), rTmp.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    //GADGET_CHECK_RETURN_FALSE(r.copyFrom(rTmp));
    r.copyFrom(rTmp);

    return true;
}

template <unsigned int D> bool sqrt(const hoNDImage<GT_Complex8, D>& x, hoNDImage<GT_Complex8, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    vcSqrt(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool minAbsolute(const hoNDImage<GT_Complex8, D>& x, GT_Complex8& r, size_t& ind)
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

template <unsigned int D> bool maxAbsolute(const hoNDImage<GT_Complex8, D>& x, GT_Complex8& r, size_t& ind)
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

template <unsigned int D> bool multiplyConj(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, hoNDImage<GT_Complex8, D>& r)
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

template <unsigned int D> bool argument(const hoNDImage<GT_Complex8, D>& x, hoNDImage<float, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    vcArg(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), r.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool conjugate(const hoNDImage<GT_Complex8, D>& x, hoNDImage<GT_Complex8, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    vcConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex8*>(x.begin()), reinterpret_cast<MKL_Complex8*>(r.begin()));
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool addEpsilon(hoNDImage<GT_Complex8, D>& x)
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

template <unsigned int D> bool norm2(const hoNDImage<GT_Complex8, D>& x, float& r)
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

template <unsigned int D> bool norm1(const hoNDImage<GT_Complex8, D>& x, float& r)
{
    try
    {
        hoNDImage<float, D> a;
        GADGET_CHECK_RETURN_FALSE(absolute(x, a));
        GADGET_CHECK_RETURN_FALSE(norm1(a, r));
    }
    catch(...)
    {
        return false;
    }

    return true;
}

template <unsigned int D> bool dotc(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& y, GT_Complex8& r)
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

template <unsigned int D> bool conv2(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& ker, hoNDImage<GT_Complex8, D>& z)
{
    try
    {
        if ( !z.dimensions_equal(&x) )
        {
            z = x;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);

        size_t num = x.get_number_of_elements()/(RO*E1);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[2];
        kerShape[0] = kerRO; kerShape[1] = kerE1;

        MKL_INT xshape[2];
        xshape[0] = RO; xshape[1] = E1;

        MKL_INT start[2];
        start[0] = kerRO/2;
        start[1] = kerE1/2;

        MKL_INT kerStride[2], xstride[2], zstride[2];
        kerStride[0] = 1; kerStride[1] = kerRO;
        xstride[0] = 1; xstride[1] = RO;
        zstride[0] = 1; zstride[1] = RO;

        const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
        const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
        MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

        if ( num == 1 )
        {
            status = vslcConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetInternalPrecision(task, VSL_CONV_PRECISION_SINGLE);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslcConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
        }
        else
        {
            status = vslcConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetInternalPrecision(task, VSL_CONV_PRECISION_SINGLE);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslcConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslConvDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in conv2(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& ker, hoNDImage<GT_Complex8, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool conv3(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& ker, hoNDImage<GT_Complex8, D>& z)
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

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);
        size_t kerE2 = ker.get_size(2);

        size_t num = x.get_number_of_elements()/(RO*E1*E2);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[3];
        kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

        MKL_INT xshape[3];
        xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

        MKL_INT start[3];
        start[0] = kerRO/2;
        start[1] = kerE1/2;
        start[2] = kerE2/2;

        MKL_INT kerStride[3], xstride[3], zstride[3];
        kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
        xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
        zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

        const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
        const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
        MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

        if ( num == 1 )
        {
            status = vslcConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetInternalPrecision(task, VSL_CONV_PRECISION_SINGLE);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslcConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
        }
        else
        {
            status = vslcConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetInternalPrecision(task, VSL_CONV_PRECISION_SINGLE);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslcConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslConvDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in conv3(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& ker, hoNDImage<GT_Complex8, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool corr2(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& ker, hoNDImage<GT_Complex8, D>& z)
{
    try
    {
        if ( !z.dimensions_equal(&x) )
        {
            z = x;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);

        size_t num = x.get_number_of_elements()/(RO*E1);

        int status;
        VSLCorrTaskPtr task;

        MKL_INT kerShape[2];
        kerShape[0] = kerRO; kerShape[1] = kerE1;

        MKL_INT xshape[2];
        xshape[0] = RO; xshape[1] = E1;

        MKL_INT start[2];
        start[0] = kerRO/2;
        start[1] = kerE1/2;

        MKL_INT decimation[2];
        decimation[0] = 1;
        decimation[1] = 1;

        MKL_INT kerStride[2], xstride[2], zstride[2];
        kerStride[0] = 1; kerStride[1] = kerRO;
        xstride[0] = 1; xstride[1] = RO;
        zstride[0] = 1; zstride[1] = RO;

        const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
        const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
        MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

        if ( num == 1 )
        {
            status = vslcCorrNewTask(&task, VSL_CORR_MODE_AUTO, 2, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetDecimation(task, decimation);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetInternalPrecision(task, VSL_CORR_PRECISION_SINGLE);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslcCorrExec(task, pKer, NULL, pX, NULL, pZ, NULL);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslCorrDeleteTask(&task);
        }
        else
        {
            status = vslcCorrNewTaskX(&task, VSL_CORR_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetDecimation(task, decimation);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetInternalPrecision(task, VSL_CORR_PRECISION_SINGLE);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslcCorrExecX(task, pX+n*RO*E1, NULL, pZ+n*RO*E1, NULL);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslCorrDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in corr2(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& ker, hoNDImage<GT_Complex8, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool corr3(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& ker, hoNDImage<GT_Complex8, D>& z)
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

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);
        size_t kerE2 = ker.get_size(2);

        size_t num = x.get_number_of_elements()/(RO*E1*E2);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[3];
        kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

        MKL_INT xshape[3];
        xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

        MKL_INT start[3];
        start[0] = kerRO/2;
        start[1] = kerE1/2;
        start[2] = kerE2/2;

        MKL_INT kerStride[3], xstride[3], zstride[3];
        kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
        xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
        zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

        const MKL_Complex8* pX = reinterpret_cast<const MKL_Complex8*>(x.begin());
        const MKL_Complex8* pKer = reinterpret_cast<const MKL_Complex8*>(ker.begin());
        MKL_Complex8* pZ = reinterpret_cast<MKL_Complex8*>(z.begin());

        if ( num == 1 )
        {
            status = vslcCorrNewTask(&task, VSL_CORR_MODE_AUTO, 3, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetInternalPrecision(task, VSL_CORR_PRECISION_SINGLE);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslcCorrExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslCorrDeleteTask(&task);
        }
        else
        {
            status = vslcCorrNewTaskX(&task, VSL_CORR_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetInternalPrecision(task, VSL_CORR_PRECISION_SINGLE);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslcCorrExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslCorrDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in corr3(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& ker, hoNDImage<GT_Complex8, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool inv(const hoNDImage<GT_Complex8, D>& x, hoNDImage<GT_Complex8, D>& r)
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
        GADGET_ERROR_MSG("Errors happened in inv(const hoNDImage<GT_Complex8, D>& x, hoNDImage<GT_Complex8, D>& r) ... ");
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------------------
// GT_Complex16
// ----------------------------------------------------------------------------------------

template <unsigned int D> bool add(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r)
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

template <unsigned int D> bool subtract(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r)
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

template <unsigned int D> bool multiply(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r)
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

template <unsigned int D> bool divide(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r)
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

template <unsigned int D> bool absolute(const hoNDImage<GT_Complex16, D>& x, hoNDImage<double, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    vzAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), r.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool absolute(const hoNDImage<GT_Complex16, D>& x, hoNDImage<GT_Complex16, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    hoNDImage<double, D> rTmp;
    rTmp.create(x.get_dimensions());

    vzAbs(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), rTmp.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    //GADGET_CHECK_RETURN_FALSE(r.copyFrom(rTmp));
r.copyFrom(rTmp);

    return true;
}

template <unsigned int D> bool sqrt(const hoNDImage<GT_Complex16, D>& x, hoNDImage<GT_Complex16, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    vzSqrt(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool minAbsolute(const hoNDImage<GT_Complex16, D>& x, GT_Complex16& r, size_t& ind)
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

template <unsigned int D> bool maxAbsolute(const hoNDImage<GT_Complex16, D>& x, GT_Complex16& r, size_t& ind)
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

template <unsigned int D> bool multiplyConj(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, hoNDImage<GT_Complex16, D>& r)
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

template <unsigned int D> bool argument(const hoNDImage<GT_Complex16, D>& x, hoNDImage<double, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    vzArg(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), r.begin());
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool conjugate(const hoNDImage<GT_Complex16, D>& x, hoNDImage<GT_Complex16, D>& r)
{
    if ( r.get_number_of_elements()!=x.get_number_of_elements())
    {
        r.create(x.get_dimensions());
    }

    vzConj(x.get_number_of_elements(), reinterpret_cast<const MKL_Complex16*>(x.begin()), reinterpret_cast<MKL_Complex16*>(r.begin()));
    GADGET_CHECK_RETURN_FALSE(vmlGetErrStatus()==0);

    return true;
}

template <unsigned int D> bool addEpsilon(hoNDImage<GT_Complex16, D>& x)
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

template <unsigned int D> bool norm2(const hoNDImage<GT_Complex16, D>& x, double& r)
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

template <unsigned int D> bool norm1(const hoNDImage<GT_Complex16, D>& x, double& r)
{
    try
    {
        hoNDImage<double, D> a;
        GADGET_CHECK_RETURN_FALSE(absolute(x, a));
        GADGET_CHECK_RETURN_FALSE(norm1(a, r));
    }
    catch(...)
    {
        return false;
    }

    return true;
}

template <unsigned int D> bool dotc(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& y, GT_Complex16& r)
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

template <unsigned int D> bool conv2(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& ker, hoNDImage<GT_Complex16, D>& z)
{
    try
    {
        if ( !z.dimensions_equal(&x) )
        {
            z = x;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);

        size_t num = x.get_number_of_elements()/(RO*E1);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[2];
        kerShape[0] = kerRO; kerShape[1] = kerE1;

        MKL_INT xshape[2];
        xshape[0] = RO; xshape[1] = E1;

        MKL_INT start[2];
        start[0] = kerRO/2;
        start[1] = kerE1/2;

        MKL_INT kerStride[2], xstride[2], zstride[2];
        kerStride[0] = 1; kerStride[1] = kerRO;
        xstride[0] = 1; xstride[1] = RO;
        zstride[0] = 1; zstride[1] = RO;

        const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
        const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
        MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

        if ( num == 1 )
        {
            status = vslzConvNewTask(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslzConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
        }
        else
        {
            status = vslzConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslzConvExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslConvDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in conv2(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& ker, hoNDImage<GT_Complex16, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool conv3(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& ker, hoNDImage<GT_Complex16, D>& z)
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

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);
        size_t kerE2 = ker.get_size(2);

        size_t num = x.get_number_of_elements()/(RO*E1*E2);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[3];
        kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

        MKL_INT xshape[3];
        xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

        MKL_INT start[3];
        start[0] = kerRO/2;
        start[1] = kerE1/2;
        start[2] = kerE2/2;

        MKL_INT kerStride[3], xstride[3], zstride[3];
        kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerRO*kerE1;
        xstride[0] = 1; xstride[1] = RO; xstride[2] = RO*E1;
        zstride[0] = 1; zstride[1] = RO; zstride[2] = RO*E1;

        const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
        const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
        MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

        if ( num == 1 )
        {
            status = vslzConvNewTask(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslzConvExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslConvDeleteTask(&task);
        }
        else
        {
            status = vslzConvNewTaskX(&task, VSL_CONV_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslConvSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslzConvExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslConvDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in conv3(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& ker, hoNDImage<GT_Complex16, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool corr2(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& ker, hoNDImage<GT_Complex16, D>& z)
{
    try
    {
        if ( !z.dimensions_equal(&x) )
        {
            z = x;
        }

        size_t RO = x.get_size(0);
        size_t E1 = x.get_size(1);

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);

        size_t num = x.get_number_of_elements()/(RO*E1);

        int status;
        VSLCorrTaskPtr task;

        MKL_INT kerShape[2];
        kerShape[0] = kerRO; kerShape[1] = kerE1;

        MKL_INT xshape[2];
        xshape[0] = RO; xshape[1] = E1;

        MKL_INT start[2];
        start[0] = kerRO/2;
        start[1] = kerE1/2;

        MKL_INT kerStride[2], xstride[2], zstride[2];
        kerStride[0] = 1; kerStride[1] = kerRO;
        xstride[0] = 1; xstride[1] = RO;
        zstride[0] = 1; zstride[1] = RO;

        const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
        const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
        MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

        if ( num == 1 )
        {
            status = vslzCorrNewTask(&task, VSL_CORR_MODE_AUTO, 2, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslzCorrExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslCorrDeleteTask(&task);
        }
        else
        {
            status = vslzCorrNewTaskX(&task, VSL_CORR_MODE_AUTO, 2, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslzCorrExecX(task, pX+n*RO*E1, xstride, pZ+n*RO*E1, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslCorrDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in corr2(const hoNDImage<GT_Complex8, D>& x, const hoNDImage<GT_Complex8, D>& ker, hoNDImage<GT_Complex8, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool corr3(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& ker, hoNDImage<GT_Complex16, D>& z)
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

        size_t kerRO = ker.get_size(0);
        size_t kerE1 = ker.get_size(1);
        size_t kerE2 = ker.get_size(2);

        size_t num = x.get_number_of_elements()/(RO*E1*E2);

        int status;
        VSLConvTaskPtr task;

        MKL_INT kerShape[3];
        kerShape[0] = kerRO; kerShape[1] = kerE1; kerShape[2] = kerE2;

        MKL_INT xshape[3];
        xshape[0] = RO; xshape[1] = E1; xshape[2] = E2;

        MKL_INT start[3];
        start[0] = kerRO/2;
        start[1] = kerE1/2;
        start[2] = kerE2/2;

        MKL_INT kerStride[3], xstride[3], zstride[3];
        kerStride[0] = 1; kerStride[1] = kerRO; kerStride[2] = kerE1;
        xstride[0] = 1; xstride[1] = RO; xstride[2] = E1;
        zstride[0] = 1; zstride[1] = RO; zstride[2] = E1;

        const MKL_Complex16* pX = reinterpret_cast<const MKL_Complex16*>(x.begin());
        const MKL_Complex16* pKer = reinterpret_cast<const MKL_Complex16*>(ker.begin());
        MKL_Complex16* pZ = reinterpret_cast<MKL_Complex16*>(z.begin());

        if ( num == 1 )
        {
            status = vslzCorrNewTask(&task, VSL_CORR_MODE_AUTO, 3, kerShape, xshape, xshape);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslzCorrExec(task, pKer, kerStride, pX, xstride, pZ, zstride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

                vslCorrDeleteTask(&task);
        }
        else
        {
            status = vslzCorrNewTaskX(&task, VSL_CORR_MODE_AUTO, 3, kerShape, xshape, xshape, pKer, kerStride);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            status = vslCorrSetStart(task, start);
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            long long n;

            #pragma omp parallel for default(none) private(n) shared(num, task, pX, RO, E1, E2, status, xstride, pZ, zstride)
            for ( n=0; n<(long long)num; n++ )
            {
                status = vslzCorrExecX(task, pX+n*RO*E1*E2, xstride, pZ+n*RO*E1*E2, zstride);
            }
            GADGET_CHECK_RETURN_FALSE(status==VSL_STATUS_OK);

            vslCorrDeleteTask(&task);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors happened in corr3(const hoNDImage<GT_Complex16, D>& x, const hoNDImage<GT_Complex16, D>& ker, hoNDImage<GT_Complex16, D>& z) ... ");
        return false;
    }

    return true;
}

template <unsigned int D> bool inv(const hoNDImage<GT_Complex16, D>& x, hoNDImage<GT_Complex16, D>& r)
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
        GADGET_ERROR_MSG("Errors happened in inv(const hoNDImage<GT_Complex16, D>& x, hoNDImage<GT_Complex16, D>& r) ... ");
        return false;
    }

    return true;
}

#endif // USE_MKL

template EXPORTCPUCOREMATH bool gaussianKernel(float sigma, double kerWidthInUnitOfSigma, double deltaKer, hoNDArray<float>& ker);
template EXPORTCPUCOREMATH bool gaussianKernel(double sigma, double kerWidthInUnitOfSigma, double deltaKer, hoNDArray<double>& ker);

#define DimImage 1
#include "hoNDImage_util_instantiate.hxx"
#undef DimImage

#define DimImage 2
#include "hoNDImage_util_instantiate.hxx"
#undef DimImage

#define DimImage 3
#include "hoNDImage_util_instantiate.hxx"
#undef DimImage

#define DimImage 4
#include "hoNDImage_util_instantiate.hxx"
#undef DimImage

#define DimImage 5
#include "hoNDImage_util_instantiate.hxx"
#undef DimImage

#define DimImage 6
#include "hoNDImage_util_instantiate.hxx"
#undef DimImage

#define DimImage 7
#include "hoNDImage_util_instantiate.hxx"
#undef DimImage

#define DimImage 8
#include "hoNDImage_util_instantiate.hxx"
#undef DimImage

#define DimImage 9
#include "hoNDImage_util_instantiate.hxx"
#undef DimImage

}
