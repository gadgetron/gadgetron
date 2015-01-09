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

            // #pragma omp parallel for default(none) private(x) shared(sx, pX, pGx)
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

            // #pragma omp parallel for default(none) private(x, y) shared(sx, sy, pX, pGx, pGy)
            for ( y=1; y<sy-1; y++ )
            {
                for ( x=1; x<sx-1; x++ )
                {
                    size_t offset = x + y*sx;

                    pGx[offset] = pX[offset+1] - pX[offset-1];
                    pGy[offset] = pX[offset+sx] - pX[offset-sx];
                }
            }

            // #pragma omp parallel for default(none) private(x) shared(sx, sy, pX, pGx, pGy)
            for ( x=1; x<sx-1; x++ )
            {
                pGx[x] = pX[x+1] - pX[x-1];

                size_t offset = x + (sy-1)*sx;
                pGx[offset] = pX[offset+1] - pX[offset-1];

                pGy[x] = pX[x+sx] - pX[x];
                pGy[x + (sy-1)*sx] = pX[x + (sy-1)*sx] - pX[x + (sy-2)*sx];
            }

            // #pragma omp parallel for default(none) private(y) shared(sx, sy, pX, pGx, pGy)
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
        GERROR_STREAM("Errors happened in gradient(const hoNDImage<T, D>& x, hoNDImage<T, D> gx[D]) ... ");
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

        T D = (T)( (deltaKer*deltaKer)/(2*sigma*sigma) );

        long long ii;
        for ( ii=-N/2; ii<=N/2; ii++ )
        {
            ker(ii+N/2) = exp( -(ii*ii*D) );
            kerSum += ker(ii+N/2);
        }

        T GNorm = (T)(1/std::sqrt(2*3.141592653579*sigma*sigma));
        GNorm /= kerSum;

        Gadgetron::scal(GNorm, ker);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in gaussianKernel(T sigma, double kerWidthInUnitOfSigma, double deltaKer, hoNDArray<T>& ker) ... ");
        return false;
    }
    return true;
}

// As well-know in the computer vision, the gaussian filter is implemented as the DERICHE filter
// therefore, the computation cost is independent from the sigma
// [1] Deriche, R., 1992, Recursively implementing the Gaussian and its derivatives: Proceedings of the 2nd International Conference on Image Processing, Singapore, p. 263–267.
// [2] http://en.wikipedia.org/wiki/Deriche_edge_detector gives details about this filter
// this implementation is based on this webpage

template <class T, class T2>
inline void DericheSmoothing(T* pData, size_t N, T* mem, T2 sigma, size_t offset=0)
{
    typedef typename realType<T>::Type real_type;

    if ( sigma < 1e-6 ) sigma = (T2)(1e-6);

    // following the note of http://en.wikipedia.org/wiki/Deriche_edge_detector

    real_type alpha = (real_type)(1.4105/sigma); // this value 1.4105 is from equation 37 of ref [1]
    real_type e_alpha = (real_type)( exp( (double)(-alpha) ) );
    real_type e_alpha_sqr = e_alpha*e_alpha;
    real_type k = ( (1-e_alpha)*(1-e_alpha) ) / ( 1 + 2*alpha*e_alpha - e_alpha_sqr );

    real_type a1 = k;
    real_type a2 = k * e_alpha * (alpha-1);
    real_type a3 = k * e_alpha * (alpha+1);
    real_type a4 = -k * e_alpha_sqr;

    real_type b1 = 2 * e_alpha;
    real_type b2 = -e_alpha_sqr;

    // compute the left to right filtering and the right to left filtering
    // for the speed, just use the zero boundary condition
    // TODO: try out other boundary conditions
    T* forward = mem;
    T* reverse = mem + N;

    if ( offset == 0 )
    {
        forward[0] = a1 * pData[0];
        reverse[N-1] = 0;

        size_t ii;

        if ( N > 1 )
        {
            forward[1] = a1 * pData[1] + a2*pData[0] + b1 * forward[0];
            reverse[N-2] = a3 * pData[N-1] + b1 * reverse[N-1];

            for ( ii=2; ii<N; ii++ )
            {
                forward[ii] = (a1*pData[ii] + a2*pData[ii-1]) + (b1*forward[ii-1] + b2*forward[ii-2]);
                reverse[N-1-ii] = (a3*pData[N-ii] + a4*pData[N-ii+1]) + (b1*reverse[N-ii] + b2*reverse[N-ii+1]);
            }
        }

        // Gadgetron::math::add(N, forward, reverse, pData);

        for ( ii=0; ii<N; ii++ )
        {
            pData[ii] = forward[ii] + reverse[ii];
        }
    }
    else
    {
        forward[0] = a1 * pData[0];
        reverse[N-1] = 0;

        if ( N > 1 )
        {
            forward[1] = a1 * pData[offset] + a2*pData[0] + b1 * forward[0];
            reverse[N-2] = a3 * pData[(N-1)*offset] + b1 * reverse[N-1];

            size_t ii;
            for ( ii=2; ii<N; ii++ )
            {
                forward[ii] = (a1*pData[ii*offset] + a2*pData[(ii-1)*offset]) + (b1*forward[ii-1] + b2*forward[ii-2]);
                reverse[N-1-ii] = (a3*pData[(N-ii)*offset] + a4*pData[(N-ii+1)*offset]) + (b1*reverse[N-ii] + b2*reverse[N-ii+1]);
            }

            for ( ii=0; ii<N; ii++ )
            {
                pData[ii*offset] = forward[ii] + reverse[ii];
            }
        }
    }
}

template<class ArrayType, class T2> 
bool filterGaussian(ArrayType& img, T2 sigma[], typename ArrayType::value_type* mem)
{
    try
    {
        typedef typename ArrayType::value_type T;

        size_t D = img.get_number_of_dimensions();

        if ( D == 1 )
        {
            if ( sigma[0] > 0 )
            {
                size_t sx = img.get_size(0);

                bool allocate = false;
                if ( mem == NULL )
                {
                    mem = new T[2*sx];
                    allocate = true;
                }

                Gadgetron::DericheSmoothing(img.begin(), sx, mem, sigma[0]);

                if ( allocate ) delete [] mem;
            }
        }
        else if ( D == 2 )
        {
            long long sx = (long long)img.get_size(0);
            long long sy = (long long)img.get_size(1);

            T* pData = img.begin();

            long long x, y;

            if ( mem != NULL )
            {
                if ( sigma[0] > 0 )
                {
                    // filter along x
                    {
                        for ( y=0; y<sy; y++ )
                        {
                            Gadgetron::DericheSmoothing(pData+y*sx, sx, mem, sigma[0]);
                        }
                    }
                }

                if ( sigma[1] > 0 )
                {
                    // filter along y
                    {
                        for ( x=0; x<sx; x++ )
                        {
                            Gadgetron::DericheSmoothing(pData+x, sy, mem, sigma[1], sx);
                        }
                    }
                }
            }
            else
            {
                if ( sigma[0] > 0 )
                {
                    // filter along x
                    // #pragma omp parallel default(none) private(y) shared(sx, sy, pData, sigma)
                    {
                        T* mem = new T[2*sx];

                        // #pragma omp for 
                        for ( y=0; y<sy; y++ )
                        {
                            Gadgetron::DericheSmoothing(pData+y*sx, sx, mem, sigma[0]);
                        }

                        delete [] mem;
                    }
                }

                if ( sigma[1] > 0 )
                {
                    // filter along y
                    //#pragma omp parallel default(none) private(x) shared(sx, sy, pData, sigma)
                    {
                        T* mem = new T[2*sy];

                        // #pragma omp for 
                        for ( x=0; x<sx; x++ )
                        {
                            Gadgetron::DericheSmoothing(pData+x, sy, mem, sigma[1], sx);
                        }

                        delete [] mem;
                    }
                }
            }
        }
        else if ( D == 3 )
        {
            long long sx = (long long)img.get_size(0);
            long long sy = (long long)img.get_size(1);
            long long sz = (long long)img.get_size(2);

            T* pData = img.begin();

            long long x, y, z;

            if ( sigma[0] > 0 )
            {
                // filter along x
                #pragma omp parallel default(none) private(y, z) shared(sx, sy, sz, pData, sigma)
                {
                    T* mem = new T[2*sx];

                    #pragma omp for 
                    for ( z=0; z<sz; z++ )
                    {
                        for ( y=0; y<sy; y++ )
                        {
                            Gadgetron::DericheSmoothing(pData+y*sx+z*sx*sy, sx, mem, sigma[0]);
                        }
                    }

                    delete [] mem;
                }
            }

            if ( sigma[1] > 0 )
            {
                // filter along y
                #pragma omp parallel default(none) private(x, y, z) shared(sx, sy, sz, pData, sigma)
                {
                    T* buf = new T[3*sy];
                    T* mem = buf + sy;

                    #pragma omp for 
                    for ( z=0; z<sz; z++ )
                    {
                        for ( x=0; x<sx; x++ )
                        {
                            size_t offset = x + z*sx*sy;

                            for ( y=0; y<sy; y++ )
                            {
                                buf[y] = pData[offset + y*sx];
                            }

                            Gadgetron::DericheSmoothing(buf, sy, mem, sigma[1]);

                            for ( y=0; y<sy; y++ )
                            {
                                pData[offset + y*sx] = buf[y];
                            }
                        }
                    }

                    delete [] buf;
                }
            }

            if ( sigma[2] > 0 )
            {
                // filter along z
                #pragma omp parallel default(none) private(x, y, z) shared(sx, sy, sz, pData, sigma)
                {
                    T* buf = new T[3*sz];
                    T* mem = buf + sz;

                    #pragma omp for 
                    for ( y=0; y<sy; y++ )
                    {
                        for ( x=0; x<sx; x++ )
                        {
                            size_t offset = x + y*sx;

                            for ( z=0; z<sz; z++ )
                            {
                                buf[z] = pData[offset + z*sx*sy];
                            }

                            Gadgetron::DericheSmoothing(buf, sz, mem, sigma[2]);

                            for ( z=0; z<sz; z++ )
                            {
                                pData[offset + z*sx*sy] = buf[z];
                            }
                        }
                    }

                    delete [] buf;
                }
            }
        }
        else if ( D == 4 )
        {
            long long sx = (long long)img.get_size(0);
            long long sy = (long long)img.get_size(1);
            long long sz = (long long)img.get_size(2);
            long long st = (long long)img.get_size(3);

            T* pData = img.begin();

            long long x, y, z, t;

            if ( sigma[0] > 0 )
            {
                // filter along x
                #pragma omp parallel default(none) private(y, z, t) shared(sx, sy, sz, st, pData, sigma)
                {
                    T* mem = new T[2*sx];

                    #pragma omp for 
                    for ( t=0; t<st; t++ )
                    {
                        for ( z=0; z<sz; z++ )
                        {
                            for ( y=0; y<sy; y++ )
                            {
                                Gadgetron::DericheSmoothing(pData+y*sx+z*sx*sy+t*sx*sy*sz, sx, mem, sigma[0]);
                            }
                        }
                    }

                    delete [] mem;
                }
            }

            if ( sigma[1] > 0 )
            {
                // filter along y
                #pragma omp parallel default(none) private(x, y, z, t) shared(sx, sy, sz, st, pData, sigma)
                {
                    T* buf = new T[3*sy];
                    T* mem = buf + sy;

                    #pragma omp for 
                    for ( t=0; t<st; t++ )
                    {
                        for ( z=0; z<sz; z++ )
                        {
                            for ( x=0; x<sx; x++ )
                            {
                                size_t offset = x + z*sx*sy + t*sx*sy*sz;

                                for ( y=0; y<sy; y++ )
                                {
                                    buf[y] = pData[offset + y*sx];
                                }

                                Gadgetron::DericheSmoothing(buf, sy, mem, sigma[1]);

                                for ( y=0; y<sy; y++ )
                                {
                                    pData[offset + y*sx] = buf[y];
                                }
                            }
                        }
                    }

                    delete [] buf;
                }
            }

            if ( sigma[2] > 0 )
            {
                // filter along z
                #pragma omp parallel default(none) private(x, y, z, t) shared(sx, sy, sz, st, pData, sigma)
                {
                    T* buf = new T[3*sz];
                    T* mem = buf + sz;

                    #pragma omp for 
                    for ( t=0; t<st; t++ )
                    {
                        for ( y=0; y<sy; y++ )
                        {
                            for ( x=0; x<sx; x++ )
                            {
                                size_t offset = x + y*sx + t*sx*sy*sz;

                                for ( z=0; z<sz; z++ )
                                {
                                    buf[z] = pData[offset + z*sx*sy];
                                }

                                Gadgetron::DericheSmoothing(buf, sz, mem, sigma[2]);

                                for ( z=0; z<sz; z++ )
                                {
                                    pData[offset + z*sx*sy] = buf[z];
                                }
                            }
                        }
                    }

                    delete [] buf;
                }
            }

            if ( sigma[3] > 0 )
            {
                // filter along t
                #pragma omp parallel default(none) private(x, y, z, t) shared(sx, sy, sz, st, pData, sigma)
                {
                    T* buf = new T[3*st];
                    T* mem = buf + st;

                    #pragma omp for 
                    for ( z=0; z<sz; z++ )
                    {
                        for ( y=0; y<sy; y++ )
                        {
                            for ( x=0; x<sx; x++ )
                            {
                                size_t offset = x + y*sx + z*sx*sy;

                                for ( t=0; t<st; t++ )
                                {
                                    buf[t] = pData[offset + t*sx*sy*sz];
                                }

                                Gadgetron::DericheSmoothing(buf, st, mem, sigma[3]);

                                for ( t=0; t<st; t++ )
                                {
                                    pData[offset + t*sx*sy*sz] = buf[t];
                                }
                            }
                        }
                    }

                    delete [] buf;
                }
            }
        }
        else
        {
            std::vector<long long> dim(D);

            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                dim[ii] = (long long)img.get_size(ii);
            }

            T* pData = img.begin();

            long long N = (long long)img.get_number_of_elements();

            std::vector<size_t> offsetFactor(D);
            img.get_offset_factor(offsetFactor);

            // filter along every dimension
            for ( ii=0; ii<D; ii++ )
            {
                if ( sigma[ii] > 0 )
                {
                    long long num = N/dim[ii];

                    long long n;

                    if ( ii == 0 )
                    {
                        #pragma omp parallel default(none) private(n) shared(num, dim, pData, sigma)
                        {
                            T* mem = new T[ 2*dim[0] ];

                            #pragma omp for 
                            for ( n=0; n<num; n++ )
                            {
                                Gadgetron::DericheSmoothing(pData+n*dim[0], dim[0], mem, sigma[0]);
                            }

                            delete [] mem;
                        }
                    }
                    else
                    {
                        std::vector<size_t> dimCurr(D-1);

                        unsigned int jj;
                        for ( jj=0; jj<D; jj++ )
                        {
                            if ( jj < ii )
                            {
                                dimCurr[jj] = dim[jj];
                            }

                            if ( jj > ii )
                            {
                                dimCurr[jj-1] = dim[jj];
                            }
                        }

                        std::vector<size_t> offsetFactorCurr(D-1);
                        NDArray<T>::calculate_offset_factors(dimCurr, offsetFactorCurr);

                        #pragma omp parallel default(none) private(n) shared(D, num, dim, img, pData, sigma, ii, offsetFactor, offsetFactorCurr)
                        {
                            T* buf = new T[ 3*dim[ii] ];
                            T* mem = buf + dim[ii];

                            std::vector<size_t> ind(D);
                            std::vector<size_t> indCurr(D-1);

                            std::vector<size_t> offset(dim[ii]);

                            #pragma omp for 
                            for ( n=0; n<num; n++ )
                            {
                                NDArray<T>::calculate_index(n, offsetFactorCurr, indCurr);

                                unsigned int jj;
                                for ( jj=0; jj<D; jj++ )
                                {
                                    if ( jj < ii )
                                    {
                                        ind[jj] = indCurr[jj];
                                    }

                                    if ( jj > ii )
                                    {
                                        ind[jj] = indCurr[jj-1];
                                    }
                                }

                                ind[ii] = 0;
                                offset[0] = img.calculate_offset(ind);
                                buf[0] = pData[ offset[0] ];

                                long long d;
                                for ( d=1; d<dim[ii]; d++ )
                                {
                                    offset[d] = offset[d-1] + offsetFactor[ii];
                                    buf[d] = pData[ offset[d] ];
                                }

                                Gadgetron::DericheSmoothing(buf, dim[ii], mem, sigma[ii]);

                                for ( d=0; d<dim[ii]; d++ )
                                {
                                    pData[ offset[d] ] = buf[d];
                                }
                            }

                            delete [] buf;
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in filterGaussian(const hoNDImage<T, D>& x, T sigma[], typename ArrayType::value_type* mem) ... ");
        return false;
    }

    return true;
}

template EXPORTCPUCOREMATH bool filterGaussian(hoNDArray<float>& img, float sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDArray<float>& img, double sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDArray<double>& img, double sigma[], double* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDArray<double>& img, float sigma[], double* mem);

template EXPORTCPUCOREMATH bool filterGaussian(hoNDArray< std::complex<float> >& img, float sigma[],  std::complex<float> * mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDArray< std::complex<double> >& img, double sigma[],  std::complex<double> * mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDArray< std::complex<float> >& img, double sigma[],  std::complex<float> * mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoNDArray< std::complex<double> >& img, float sigma[],  std::complex<double> * mem);

template EXPORTCPUCOREMATH bool filterGaussian(ho2DArray<float>& img, float sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(ho2DArray<float>& img, double sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(ho2DArray<double>& img, float sigma[], double* mem);
template EXPORTCPUCOREMATH bool filterGaussian(ho2DArray<double>& img, double sigma[], double* mem);

template EXPORTCPUCOREMATH bool filterGaussian(hoMatrix<float>& img, float sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoMatrix<float>& img, double sigma[], float* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoMatrix<double>& img, double sigma[], double* mem);
template EXPORTCPUCOREMATH bool filterGaussian(hoMatrix<double>& img, float sigma[], double* mem);

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
