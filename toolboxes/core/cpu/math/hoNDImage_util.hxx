/** \file   hoNDImage_util.hxx
    \brief  operations on the hoNDImage class.
*/

namespace Gadgetron
{
    template <typename ImageType> 
    bool corrCoef(const ImageType& a, const ImageType& b, typename ImageType::value_type& r)
    {
        typedef typename ImageType::value_type T;
        unsigned int D = ImageType::NDIM;

        try
        {
            GADGET_CHECK_RETURN_FALSE(a.dimensions_equal(&b));

            r = -1;

            T ma, mb;
            ma = Gadgetron::mean( const_cast< ImageType* >(&a) );
            mb = Gadgetron::mean( const_cast< ImageType* >(&b) );

            size_t N = a.get_number_of_elements();

            const T* pA = a.begin();
            const T* pB = b.begin();

            size_t n;

            double x(0), y(0), z(0);
            for ( n=0; n<N; n++ )
            {
                x += (pA[n]-ma)*(pA[n]-ma);
                y += (pB[n]-mb)*(pB[n]-mb);
                z += (pA[n]-ma)*(pB[n]-mb);
            }

            double p = std::sqrt(x*y);
            if ( p > 0 )
            {
                r = (T)(z/p);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in corrCoef(const ImageType& a, const ImageType& b, T& r) ... ");
            return false;
        }

        return true;
    }

    template<typename ImageType, typename InterpolatorType> 
    bool downsampleImage(const ImageType& in, InterpolatorType& interp, ImageType& out, float ratio[])
    {
        unsigned int D = ImageType::NDIM;

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
            GERROR_STREAM("Errors happened in downsampleImage(const ImageType& in, InterpolatorType& interp, ImageType& out, float ratio[]) ... ");
            return false;
        }

        return true;
    }

    template<typename ImageType, typename InterpolatorType> 
    bool upsampleImage(const ImageType& in, InterpolatorType& interp, ImageType& out, float ratio[])
    {
        unsigned int D = ImageType::NDIM;

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
            GERROR_STREAM("Errors happened in upsampleImage(const ImageType& in, InterpolatorType& interp, ImageType& out, float ratio[]) ... ");
            return false;
        }

        return true;
    }

    template<typename ImageType, typename InterpolatorType> 
    bool resampleImage(const ImageType& in, InterpolatorType& interp, const std::vector<size_t>& dim_out, ImageType& out)
    {
        unsigned int D = ImageType::NDIM;

        try
        {
            typedef typename ImageType::coord_type coord_type;

            /// get the coordinate parameters
            std::vector<size_t> dim;
            in.get_dimensions(dim);

            std::vector<coord_type> pixelSize;
            in.get_pixel_size(pixelSize);

            std::vector<coord_type> origin;
            in.get_origin(origin);

            typename ImageType::axis_type axis;
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
            interp.setArray( const_cast< ImageType& >(in) );

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

                #pragma omp parallel private(n) shared(N, in, out, interp)
                {
                    std::vector<size_t> ind_o(ImageType::NDIM);
                    std::vector<coord_type> ind_i(ImageType::NDIM);

                    std::vector<coord_type> pos(ImageType::NDIM);

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
            GERROR_STREAM("Errors happened in resampleImage(const ImageType& in, InterpolatorType& interp, ImageType& out, size_t size_out[D]) ... ");
            return false;
        }

        return true;
    }

    template<typename ImageType, typename BoundaryHandlerType> 
    bool downsampleImageBy2WithAveraging(const ImageType& in, BoundaryHandlerType& bh, ImageType& out)
    {
        typedef typename ImageType::value_type T;
        unsigned int D = ImageType::NDIM;

        try
        {
            typedef typename ImageType::coord_type coord_type;

            bh.setArray( const_cast< ImageType& >(in) );

            /// get the coordinate parameters
            std::vector<size_t> dim;
            in.get_dimensions(dim);

            std::vector<coord_type> pixelSize;
            in.get_pixel_size(pixelSize);

            std::vector<coord_type> origin;
            in.get_origin(origin);

            typename ImageType::axis_type axis;
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
                size_t sx = dim_out[0];
                size_t sy = dim_out[1];

                T weight = 1.0/5;

                long long x, y;

                #pragma omp parallel for default(none) private(x, y) shared(sx, sy, bh, out)
                for ( y=0; y<(long long)sy; y++ )
                {
                    long long iy = y<<1;

                    for ( x=0; x<(long long)sx; x++ )
                    {
                        long long ix = x<<1;
                        out( (size_t)(x+y*sx) ) = bh(ix, iy) + ( bh(ix+1, iy) + bh(ix-1, iy) ) + ( bh(ix, iy+1) + bh(ix, iy-1) );
                    }
                }

                Gadgetron::scal(weight, out);
            }
            else if ( D == 3 )
            {
                size_t sx = dim_out[0];
                size_t sy = dim_out[1];
                size_t sz = dim_out[2];

                T weight = 1.0/7;

                long long x, y, z;

                #pragma omp parallel for default(none) private(x, y, z) shared(sx, sy, sz, bh, out)
                for ( z=0; z<sz; z++ )
                {
                    long long iz = z<<1;

                    for ( y=0; y<sy; y++ )
                    {
                        long long iy = y<<1;

                        size_t offset = y*sx + z*sx*sy;

                        for ( x=0; x<sx; x++ )
                        {
                            long long ix = x<<1;

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
                size_t sx = dim_out[0];
                size_t sy = dim_out[1];
                size_t sz = dim_out[2];
                size_t st = dim_out[3];

                T weight = 1.0/9;

                long long x, y, z, t;

                #pragma omp parallel for default(none) private(x, y, z, t) shared(sx, sy, sz, st, bh, out)
                for ( t=0; t<st; t++ )
                {
                    long long it = t<<1;

                    for ( z=0; z<sz; z++ )
                    {
                        long long iz = z<<1;

                        for ( y=0; y<sy; y++ )
                        {
                            long long iy = y<<1;

                            size_t offset = y*sx + z*sx*sy + t*sx*sy*sz;

                            for ( x=0; x<sx; x++ )
                            {
                                long long ix = x<<1;

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

                size_t N = out.get_number_of_elements();

                long long n;

                #pragma omp parallel private(n) shared(N, bh, out, dim_out)
                {
                    std::vector<size_t> ind_out(ImageType::NDIM);
                    std::vector<long long> ind_in(ImageType::NDIM);

                    #pragma omp for 
                    for ( n=0; n<N; n++ )
                    {
                        out.calculate_index(n, ind_out);

                        unsigned int ii;
                        for ( ii=0; ii<ImageType::NDIM; ii++ )
                        {
                            ind_in[ii] = ind_out[ii]<<1;
                        }

                        T v = bh(ind_in);

                        for ( ii=0; ii<ImageType::NDIM; ii++ )
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
            GERROR_STREAM("Errors happened in downsampleImageBy2WithAveraging(const ImageType& in, ImageType& out) ... ");
            return false;
        }

        return true;
    }

    template<typename ImageType, typename BoundaryHandlerType> 
    bool expandImageBy2(const ImageType& in, BoundaryHandlerType& bh, ImageType& out)
    {
        typedef typename ImageType::value_type T;
        unsigned int D = ImageType::NDIM;

        try
        {
            typedef typename ImageType::coord_type coord_type;

            bh.setArray( const_cast< ImageType& >(in) );

            /// get the coordinate parameters
            std::vector<size_t> dim;
            in.get_dimensions(dim);

            std::vector<coord_type> pixelSize;
            in.get_pixel_size(pixelSize);

            std::vector<coord_type> origin;
            in.get_origin(origin);

            typename ImageType::axis_type axis;
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
                size_t sx = dim[0];
                size_t sy = dim[1];

                long long x, y;

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
                size_t sx_out = out.get_size(0);
                size_t sy_out = out.get_size(1);

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
                size_t sx = dim[0];
                size_t sy = dim[1];
                size_t sz = dim[2];

                long long x, y, z;

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
                size_t sx_out = out.get_size(0);
                size_t sy_out = out.get_size(1);
                size_t sz_out = out.get_size(2);

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
                hoNDInterpolatorLinear<ImageType > interp(const_cast< ImageType& >(in), bh);

                size_t N = out.get_number_of_elements();

                long long n;

                #pragma omp parallel private(n) shared(N, bh, in, out, interp)
                {
                    std::vector<size_t> ind_out(ImageType::NDIM);
                    std::vector<coord_type> ind_in(ImageType::NDIM);

                    #pragma omp for 
                    for ( n=0; n<N; n++ )
                    {
                        out.calculate_index(n, ind_out);

                        unsigned int ii;
                        for ( ii=0; ii<ImageType::NDIM; ii++ )
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
            GERROR_STREAM("Errors happened in expandImageBy2(const ImageType& in, BoundaryHandlerType& bh, ImageType& out) ... ");
            return false;
        }

        return true;
    }

    template<class ArrayType> 
    bool filterMedian(const ArrayType& img, size_t w[], ArrayType& img_out)
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

                #pragma omp parallel private(n, m, t) shared(halfW, N, img, img_out)
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
                #pragma omp parallel private(x, y, tx, ty, hx, hy) shared(halfX, halfY, sx, sy, WX, WY, pImg, pImgOut, medianInd)
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

                for ( y=sy-halfY; y<sy; y++ )
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
                #pragma omp parallel private(x, y, z, tx, ty, tz, hx, hy, hz) shared(halfX, halfY, halfZ, sx, sy, sz, WX, WY, WZ, pImg, pImgOut, medianInd)
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
            GERROR_STREAM("Errors happened in filterMedian(const ArrayType& img, size_t w[], ArrayType& img_out) ... ");
            return false;
        }

        return true;
    }

    template<class ArrayType> 
    bool filter1D(const ArrayType& img, const hoNDArray<typename realType<typename ArrayType::value_type>::Type>& ker, GT_BOUNDARY_CONDITION bh, ArrayType& img_out)
    {
        try
        {
            typedef typename ArrayType::value_type T;
            typedef typename realType<T>::Type real_value_type;

            typedef hoNDImage<T, 1> Array1DType;

            long long RO = (long long)img.get_size(0);
            long long num = (long long)(img.get_number_of_elements()/RO);

            img_out = img;

            long long kerLen = (long long)ker.get_size(0);
            long long kerHalfLen = kerLen/2;

            const real_value_type* pKer = ker.begin();

            long long ii;

            #pragma omp parallel default(none) private(ii) shared(bh, num, RO, img, img_out, kerLen, kerHalfLen, pKer)
            {
                hoNDBoundaryHandler<Array1DType>* pBH = NULL;

                hoNDBoundaryHandlerFixedValue<Array1DType> bhFixedValue;
                hoNDBoundaryHandlerBorderValue<Array1DType> bhBorderValue;
                hoNDBoundaryHandlerPeriodic<Array1DType> bhPeriodic;
                hoNDBoundaryHandlerMirror<Array1DType> bhMirror;

                pBH = &bhBorderValue;

                if ( bh == GT_BOUNDARY_CONDITION_FIXEDVALUE )
                {
                    pBH = &bhFixedValue;
                }
                else if ( bh == GT_BOUNDARY_CONDITION_BORDERVALUE )
                {
                    pBH = &bhBorderValue;
                }
                else if ( bh == GT_BOUNDARY_CONDITION_PERIODIC )
                {
                    pBH = &bhPeriodic;
                }
                else if ( bh == GT_BOUNDARY_CONDITION_MIRROR )
                {
                    pBH = &bhMirror;
                }

                #pragma omp for 
                for ( ii=0; ii<num; ii++ )
                {
                    Array1DType img1D(RO, const_cast<T*>(img.begin()+ii*RO));
                    pBH->setArray(img1D);

                    Array1DType img_out1D(RO, img_out.begin()+ii*RO);

                    long long k, j;
                    for ( k=0; k<RO; k++ )
                    {
                        T v = 0;
                        for ( j=0; j<kerLen; j++ )
                        {
                            v += (*pBH)(k+j-kerHalfLen) * pKer[kerLen-j-1];
                        }

                        img_out1D(k) = v;
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in filter1D(const ArrayType& img, const hoNDArray<T>& ker, GT_BOUNDARY_CONDITION bh, ArrayType& img_out) ... ");
            return false;
        }

        return true;
    }

    template<typename ImageType> 
    bool gradient(const ImageType& x, ImageType gx[])
    {
        typedef typename ImageType::value_type T;
        unsigned int D = ImageType::NDIM;

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
                    size_t ind[ImageType::NDIM];
                    size_t ind_positive[ImageType::NDIM];
                    size_t ind_negative[ImageType::NDIM];
                    bool inside = true;
                    unsigned int ii;

#pragma omp for 
                    for ( n=0; n<(long long)N; n++ )
                    {
                        x.calculate_index(n, ind);

                        inside = true;
                        for ( ii=0; ii<ImageType::NDIM; ii++ )
                        {
                            if ( ind[ii]==0 || ind[ii]==dim[ii]-1 )
                            {
                                inside = false;
                                break;
                            }
                        }

                        if ( inside )
                        {
                            for ( ii=0; ii<ImageType::NDIM; ii++ )
                            {
                                memcpy(ind_positive, ind, sizeof(size_t)*ImageType::NDIM);
                                memcpy(ind_negative, ind, sizeof(size_t)*ImageType::NDIM);

                                ind_positive[ii] = ind[ii] + 1;
                                ind_negative[ii] = ind[ii] - 1;

                                gx[ii](n) = x(ind_positive) - x(ind_negative);
                            }
                        }
                        else
                        {
                            for ( ii=0; ii<ImageType::NDIM; ii++ )
                            {
                                memcpy(ind_positive, ind, sizeof(size_t)*ImageType::NDIM);
                                memcpy(ind_negative, ind, sizeof(size_t)*ImageType::NDIM);

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
            GERROR_STREAM("Errors happened in gradient(const ImageType& x, ImageType gx[D]) ... ");
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
                ker(ii+N/2) = std::exp( -(ii*ii*D) );
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
    // [1] Deriche, R., 1992, Recursively implementing the Gaussian and its derivatives: Proceedings of the 2nd International Conference on Image Processing, Singapore, p. 263�267.
    // [2] http://en.wikipedia.org/wiki/Deriche_edge_detector gives details about this filter
    // this implementation is based on this webpage

    template <class T, class T2>
    inline void DericheSmoothing(T* pData, size_t N, T* mem, T2 sigma, size_t offset=0)
    {
        typedef typename realType<T>::Type real_type;

        if ( sigma < 1e-6 ) sigma = (T2)(1e-6);

        // following the note of http://en.wikipedia.org/wiki/Deriche_edge_detector

        real_type alpha = (real_type)(1.4105/sigma); // this value 1.4105 is from equation 37 of ref [1]
        real_type e_alpha = (real_type)( std::exp( (double)(-alpha) ) );
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


                for (auto ii=0; ii<D; ii++ )
                {
                    dim[ii] = (long long)img.get_size(ii);
                }

                T* pData = img.begin();

                long long N = (long long)img.get_number_of_elements();

                std::vector<size_t> offsetFactor(D);
                img.get_offset_factor(offsetFactor);

                // filter along every dimension
                for (auto ii=0; ii<D; ii++ )
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

#pragma omp parallel private(n) shared(D, num, dim, img, pData, sigma, ii, offsetFactor, offsetFactorCurr)
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
            GERROR_STREAM("Errors happened in filterGaussian(const ImageType& x, T sigma[], typename ArrayType::value_type* mem) ... ");
            return false;
        }

        return true;
    }
}
