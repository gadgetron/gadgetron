/** \file       hoNDInterpolatorLinear.h
    \brief      N-dimensional linear interpolator

                Designed to work with hoNDArray and hoNDImage

    \author     Hui Xue
*/

#ifdef _WIN32
    #include "malloc.h"
#else
    #include "alloca.h"
#endif // _WIN32

namespace Gadgetron
{
    /// hoNDInterpolatorLinear

    template <typename ArrayType> 
    typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( const coord_type* pos )
    {
        unsigned int D = array_->get_number_of_dimensions();

        gt_index_type* anchor = reinterpret_cast<gt_index_type*>(alloca(D * sizeof(gt_index_type)));
        coord_type* weights = reinterpret_cast<coord_type*>(alloca(D * sizeof(coord_type)));
        coord_type* weightsMinusOne = reinterpret_cast<coord_type*>(alloca(D * sizeof(coord_type)));

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            anchor[ii] = static_cast<gt_index_type>(std::floor(pos[ii]));
            weights[ii] = pos[ii] - anchor[ii];
            weightsMinusOne[ii] = coord_type(1.0) - weights[ii];
        }

        T res(0);

        coord_type weightAll(1.0);

        bool inRange = true;
        for ( ii=0; ii<D; ii++ )
        {
            if ( anchor[ii]<0 || anchor[ii]>=array_->get_size(ii)-1 )
            {
                inRange = false;
                break;
            }
        }

        if( inRange )
        {
            std::vector<size_t> ind(D);

            unsigned int n;
            for ( n=0; n<number_of_points_; n++ )
            {
                unsigned int lastDigit = n;
                weightAll = coord_type(1.0);

                for ( ii=0; ii<D; ii++ )
                {
                    if ( lastDigit & 1 )
                    {
                        ind[ii] = anchor[ii]+1;
                        weightAll *= weights[ii];
                    }
                    else
                    {
                        ind[ii] = anchor[ii];
                        weightAll *= weightsMinusOne[ii];
                    }

                    // shift one digit
                    lastDigit >>= 1;
                }

                res += weightAll * (*array_)(ind);
            }
        }
        else
        {
            std::vector<gt_index_type> ind(D);

            unsigned int n;
            for ( n=0; n<number_of_points_; n++ )
            {
                unsigned int lastDigit = n;
                weightAll = coord_type(1.0);

                for ( ii=0; ii<D; ii++ )
                {
                    if ( lastDigit & 1 )
                    {
                        ind[ii] = anchor[ii]+1;
                        weightAll *= weights[ii];
                    }
                    else
                    {
                        ind[ii] = anchor[ii];
                        weightAll *= weightsMinusOne[ii];
                    }

                    // shift one digit
                    lastDigit >>= 1;
                }

                res += weightAll * (*bh_)(ind);
            }
        }

        return res;
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( const std::vector<coord_type>& pos )
    {
        return this->operator()(&pos[0]);
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( coord_type x )
    {
        gt_index_type ix = static_cast<gt_index_type>(std::floor(x));
        coord_type dx = x - ix;

        if ( ix>=0 && ix<(gt_index_type)sx_-1 )
        {
            return ( (*array_)( size_t(ix) )*(1-dx) + (*array_)( size_t(ix)+1 )*dx );
        }
        else
        {
            return ( (*bh_)(ix)*(1-dx) + (*bh_)(ix+1)*dx );
        }
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( coord_type x, coord_type y )
    {
        gt_index_type ix = static_cast<gt_index_type>(std::floor(x));
        coord_type dx = x - ix;
        coord_type dx_prime = coord_type(1.0)-dx;

        gt_index_type iy = static_cast<gt_index_type>(std::floor(y));
        coord_type dy = y - iy;
        coord_type dy_prime = coord_type(1.0)-dy;

        if ( ix>=0 && ix<sx_-1 && iy>=0 && iy<sy_-1 )
        {
            size_t offset = ix + iy*sx_;
            T* data = array_->begin();

            //return (    ( data_[offset]   *   dx_prime     *dy_prime 
            //        +   data_[offset+1]   *   dx           *dy_prime)
            //        +   (data_[offset+sx_]   *   dx_prime     *dy
            //        +   data_[offset+sx_+1]   *   dx           *dy) );

            /*return (    ((*array_)(size_t(ix), size_t(iy)       )   *   dx_prime     *dy_prime
                    +   (*array_)(size_t(ix)+1, size_t(iy)      )   *   dx           *dy_prime)
                    +   ((*array_)(size_t(ix), size_t(iy)+1     )   *   dx_prime     *dy
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1    )   *   dx           *dy) );*/

            return (    (data[offset]       *   dx_prime     *dy_prime
                    +   data[offset+1]      *   dx           *dy_prime)
                    +   (data[offset+sx_]   *   dx_prime     *dy
                    +   data[offset+sx_+1]  *   dx           *dy) );
        }
        else
        {
            return (    ((*bh_)(ix, iy       )   *   dx_prime    *dy_prime 
                    +   (*bh_)(ix+1, iy      )   *   dx          *dy_prime)
                    +   ((*bh_)(ix, iy+1     )   *   dx_prime    *dy
                    +   (*bh_)(ix+1, iy+1    )   *   dx          *dy) );
        }
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( coord_type x, coord_type y, coord_type z )
    {
        gt_index_type ix = static_cast<gt_index_type>(std::floor(x));
        coord_type dx = x - ix;
        coord_type dx_prime = coord_type(1.0)-dx;

        gt_index_type iy = static_cast<gt_index_type>(std::floor(y));
        coord_type dy = y - iy;
        coord_type dy_prime = coord_type(1.0)-dy;

        gt_index_type iz = static_cast<gt_index_type>(std::floor(z));
        coord_type dz = z - iz;
        coord_type dz_prime = coord_type(1.0)-dz;

        if ( ix>=0 && ix<sx_-1 
            && iy>=0 && iy<sy_-1 
            && iz>=0 && iz<sz_-1 )
        {
            /*return (    ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)   )   *   dx_prime     *dy_prime   *dz_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)    )   *   dx           *dy_prime   *dz_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)   )   *   dx_prime     *dy         *dz_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)    )   *   dx           *dy         *dz_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1 )   *   dx_prime     *dy_prime   *dz 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1  )   *   dx           *dy_prime   *dz) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1 )   *   dx_prime     *dy         *dz 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1  )   *   dx           *dy         *dz) );*/

            size_t offset = ix + iy*sx_ + iz*sx_*sy_;

            return (    (data_[offset]              *   dx_prime     *dy_prime   *dz_prime 
                    +   data_[offset+1]             *   dx           *dy_prime   *dz_prime) 
                    +   (data_[offset+sx_]          *   dx_prime     *dy         *dz_prime 
                    +   data_[offset+sx_+1]         *   dx           *dy         *dz_prime) 
                    +   (data_[offset+sx_*sy_]      *   dx_prime     *dy_prime   *dz 
                    +   data_[offset+sx_*sy_+1]     *   dx           *dy_prime   *dz) 
                    +   (data_[offset+sx_*sy_+sx_]  *   dx_prime     *dy         *dz 
                    +   data_[offset+sx_*sy_+sx_+1] *   dx           *dy         *dz) );
        }
        else
        {
            return (    ((*bh_)(ix,   iy,     iz   )   *   dx_prime     *dy_prime   *dz_prime 
                    +   (*bh_)(ix+1, iy,     iz    )   *   dx           *dy_prime   *dz_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz   )   *   dx_prime     *dy         *dz_prime 
                    +   (*bh_)(ix+1, iy+1,   iz    )   *   dx           *dy         *dz_prime) 
                    +   ((*bh_)(ix,   iy,     iz+1 )   *   dx_prime     *dy_prime   *dz 
                    +   (*bh_)(ix+1, iy,     iz+1  )   *   dx           *dy_prime   *dz) 
                    +   ((*bh_)(ix,   iy+1,   iz+1 )   *   dx_prime     *dy         *dz 
                    +   (*bh_)(ix+1, iy+1,   iz+1  )   *   dx           *dy         *dz) );
        }
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s )
    {
        gt_index_type ix = static_cast<gt_index_type>(std::floor(x));
        coord_type dx = x - ix;
        coord_type dx_prime = coord_type(1.0)-dx;

        gt_index_type iy = static_cast<gt_index_type>(std::floor(y));
        coord_type dy = y - iy;
        coord_type dy_prime = coord_type(1.0)-dy;

        gt_index_type iz = static_cast<gt_index_type>(std::floor(z));
        coord_type dz = z - iz;
        coord_type dz_prime = coord_type(1.0)-dz;

        gt_index_type is = static_cast<gt_index_type>(std::floor(s));
        coord_type ds = s - is;
        coord_type ds_prime = coord_type(1.0)-ds;

        if ( ix>=0 && ix<(gt_index_type)array_->get_size(0)-1 
            && iy>=0 && iy<(gt_index_type)array_->get_size(1)-1 
            && iz>=0 && iz<(gt_index_type)array_->get_size(2)-1 
            && is>=0 && is<(gt_index_type)array_->get_size(3)-1 )
        {
            return (    ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is) )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is) )   *   dx           *dy_prime   *dz_prime   *ds_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is) )   *   dx_prime     *dy         *dz_prime   *ds_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is) )   *   dx           *dy         *dz_prime   *ds_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is) )   *   dx_prime     *dy_prime   *dz         *ds_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is) )   *   dx           *dy_prime   *dz         *ds_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is) )   *   dx_prime     *dy         *dz         *ds_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is) )   *   dx           *dy         *dz         *ds_prime)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is)+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is)+1 )   *   dx           *dy_prime   *dz_prime   *ds) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is)+1 )   *   dx_prime     *dy         *dz_prime   *ds 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is)+1 )   *   dx           *dy         *dz_prime   *ds) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is)+1 )   *   dx_prime     *dy_prime   *dz         *ds 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is)+1 )   *   dx           *dy_prime   *dz         *ds) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is)+1 )   *   dx_prime     *dy         *dz         *ds 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is)+1 )   *   dx           *dy         *dz         *ds) );
        }
        else
        {
            return (    ((*bh_)(ix,   iy,     iz,    is )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime 
                    +   (*bh_)(ix+1, iy,     iz,     is )   *   dx           *dy_prime   *dz_prime   *ds_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is )   *   dx_prime     *dy         *dz_prime   *ds_prime 
                    +   (*bh_)(ix+1, iy+1,   iz,     is )   *   dx           *dy         *dz_prime   *ds_prime) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is )   *   dx_prime     *dy_prime   *dz         *ds_prime 
                    +   (*bh_)(ix+1, iy,     iz+1,   is )   *   dx           *dy_prime   *dz         *ds_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is )   *   dx_prime     *dy         *dz         *ds_prime 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is )   *   dx           *dy         *dz         *ds_prime)
                    +   ((*bh_)(ix,   iy,     iz,    is+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds 
                    +   (*bh_)(ix+1, iy,     iz,     is+1 )   *   dx           *dy_prime   *dz_prime   *ds) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is+1 )   *   dx_prime     *dy         *dz_prime   *ds 
                    +   (*bh_)(ix+1, iy+1,   iz,     is+1 )   *   dx           *dy         *dz_prime   *ds) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is+1 )   *   dx_prime     *dy_prime   *dz         *ds 
                    +   (*bh_)(ix+1, iy,     iz+1,   is+1 )   *   dx           *dy_prime   *dz         *ds) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is+1 )   *   dx_prime     *dy         *dz         *ds 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is+1 )   *   dx           *dy         *dz         *ds) );
        }
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p )
    {
        gt_index_type ix = static_cast<gt_index_type>(std::floor(x));
        coord_type dx = x - ix;
        coord_type dx_prime = coord_type(1.0)-dx;

        gt_index_type iy = static_cast<gt_index_type>(std::floor(y));
        coord_type dy = y - iy;
        coord_type dy_prime = coord_type(1.0)-dy;

        gt_index_type iz = static_cast<gt_index_type>(std::floor(z));
        coord_type dz = z - iz;
        coord_type dz_prime = coord_type(1.0)-dz;

        gt_index_type is = static_cast<gt_index_type>(std::floor(s));
        coord_type ds = s - is;
        coord_type ds_prime = coord_type(1.0)-ds;

        gt_index_type ip = static_cast<gt_index_type>(std::floor(p));
        coord_type dp = p - ip;
        coord_type dp_prime = coord_type(1.0)-dp;

        if ( ix>=0 && ix<(gt_index_type)array_->get_size(0)-1 
            && iy>=0 && iy<(gt_index_type)array_->get_size(1)-1 
            && iz>=0 && iz<(gt_index_type)array_->get_size(2)-1 
            && is>=0 && is<(gt_index_type)array_->get_size(3)-1 
            && ip>=0 && ip<(gt_index_type)array_->get_size(4)-1 )
        {
            return (    ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is),     size_t(ip) )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp_prime
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is),     size_t(ip) )   *   dx           *dy_prime   *dz_prime   *ds_prime   *dp_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is),     size_t(ip) )   *   dx_prime     *dy         *dz_prime   *ds_prime   *dp_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is),     size_t(ip) )   *   dx           *dy         *dz_prime   *ds_prime   *dp_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is),     size_t(ip) )   *   dx_prime     *dy_prime   *dz         *ds_prime   *dp_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is),     size_t(ip) )   *   dx           *dy_prime   *dz         *ds_prime   *dp_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is),     size_t(ip) )   *   dx_prime     *dy         *dz         *ds_prime   *dp_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is),     size_t(ip) )   *   dx           *dy         *dz         *ds_prime   *dp_prime)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is)+1,   size_t(ip) )   *   dx_prime     *dy_prime   *dz_prime   *ds   *dp_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is)+1,   size_t(ip) )   *   dx           *dy_prime   *dz_prime   *ds   *dp_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is)+1,   size_t(ip) )   *   dx_prime     *dy         *dz_prime   *ds   *dp_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is)+1,   size_t(ip) )   *   dx           *dy         *dz_prime   *ds   *dp_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is)+1,   size_t(ip) )   *   dx_prime     *dy_prime   *dz         *ds   *dp_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is)+1,   size_t(ip) )   *   dx           *dy_prime   *dz         *ds   *dp_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is)+1,   size_t(ip) )   *   dx_prime     *dy         *dz         *ds   *dp_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is)+1,   size_t(ip) )   *   dx           *dy         *dz         *ds   *dp_prime)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is),     size_t(ip)+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is),     size_t(ip)+1 )   *   dx           *dy_prime   *dz_prime   *ds_prime   *dp) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is),     size_t(ip)+1 )   *   dx_prime     *dy         *dz_prime   *ds_prime   *dp 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is),     size_t(ip)+1 )   *   dx           *dy         *dz_prime   *ds_prime   *dp) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is),     size_t(ip)+1 )   *   dx_prime     *dy_prime   *dz         *ds_prime   *dp 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is),     size_t(ip)+1 )   *   dx           *dy_prime   *dz         *ds_prime   *dp) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is),     size_t(ip)+1 )   *   dx_prime     *dy         *dz         *ds_prime   *dp 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is),     size_t(ip)+1 )   *   dx           *dy         *dz         *ds_prime   *dp)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is)+1,     size_t(ip)+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds   *dp 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is)+1,     size_t(ip)+1 )   *   dx           *dy_prime   *dz_prime   *ds   *dp) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is)+1,     size_t(ip)+1 )   *   dx_prime     *dy         *dz_prime   *ds   *dp 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is)+1,     size_t(ip)+1 )   *   dx           *dy         *dz_prime   *ds   *dp) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is)+1,     size_t(ip)+1 )   *   dx_prime     *dy_prime   *dz         *ds   *dp 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is)+1,     size_t(ip)+1 )   *   dx           *dy_prime   *dz         *ds   *dp) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is)+1,     size_t(ip)+1 )   *   dx_prime     *dy         *dz         *ds   *dp 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is)+1,     size_t(ip)+1 )   *   dx           *dy         *dz         *ds   *dp) );
        }
        else
        {
            return (    ((*bh_)(ix,   iy,     iz,    is,     ip )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp_prime
                    +   (*bh_)(ix+1, iy,     iz,     is,     ip )   *   dx           *dy_prime   *dz_prime   *ds_prime   *dp_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is,     ip )   *   dx_prime     *dy         *dz_prime   *ds_prime   *dp_prime 
                    +   (*bh_)(ix+1, iy+1,   iz,     is,     ip )   *   dx           *dy         *dz_prime   *ds_prime   *dp_prime) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is,     ip )   *   dx_prime     *dy_prime   *dz         *ds_prime   *dp_prime 
                    +   (*bh_)(ix+1, iy,     iz+1,   is,     ip )   *   dx           *dy_prime   *dz         *ds_prime   *dp_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is,     ip )   *   dx_prime     *dy         *dz         *ds_prime   *dp_prime 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is,     ip )   *   dx           *dy         *dz         *ds_prime   *dp_prime)
                    +   ((*bh_)(ix,   iy,     iz,    is+1,   ip )   *   dx_prime     *dy_prime   *dz_prime   *ds   *dp_prime 
                    +   (*bh_)(ix+1, iy,     iz,     is+1,   ip )   *   dx           *dy_prime   *dz_prime   *ds   *dp_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is+1,   ip )   *   dx_prime     *dy         *dz_prime   *ds   *dp_prime 
                    +   (*bh_)(ix+1, iy+1,   iz,     is+1,   ip )   *   dx           *dy         *dz_prime   *ds   *dp_prime) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is+1,   ip )   *   dx_prime     *dy_prime   *dz         *ds   *dp_prime 
                    +   (*bh_)(ix+1, iy,     iz+1,   is+1,   ip )   *   dx           *dy_prime   *dz         *ds   *dp_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is+1,   ip )   *   dx_prime     *dy         *dz         *ds   *dp_prime 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is+1,   ip )   *   dx           *dy         *dz         *ds   *dp_prime)
                    +   ((*bh_)(ix,   iy,     iz,    is,     ip+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp
                    +   (*bh_)(ix+1, iy,     iz,     is,     ip+1 )   *   dx           *dy_prime   *dz_prime   *ds_prime   *dp) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is,     ip+1 )   *   dx_prime     *dy         *dz_prime   *ds_prime   *dp 
                    +   (*bh_)(ix+1, iy+1,   iz,     is,     ip+1 )   *   dx           *dy         *dz_prime   *ds_prime   *dp) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is,     ip+1 )   *   dx_prime     *dy_prime   *dz         *ds_prime   *dp 
                    +   (*bh_)(ix+1, iy,     iz+1,   is,     ip+1 )   *   dx           *dy_prime   *dz         *ds_prime   *dp) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is,     ip+1 )   *   dx_prime     *dy         *dz         *ds_prime   *dp 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is,     ip+1 )   *   dx           *dy         *dz         *ds_prime   *dp)
                    +   ((*bh_)(ix,   iy,     iz,    is+1,   ip+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds   *dp 
                    +   (*bh_)(ix+1, iy,     iz,     is+1,   ip+1 )   *   dx           *dy_prime   *dz_prime   *ds   *dp) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is+1,   ip+1 )   *   dx_prime     *dy         *dz_prime   *ds   *dp 
                    +   (*bh_)(ix+1, iy+1,   iz,     is+1,   ip+1 )   *   dx           *dy         *dz_prime   *ds   *dp) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is+1,   ip+1 )   *   dx_prime     *dy_prime   *dz         *ds   *dp 
                    +   (*bh_)(ix+1, iy,     iz+1,   is+1,   ip+1 )   *   dx           *dy_prime   *dz         *ds   *dp) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is+1,   ip+1 )   *   dx_prime     *dy         *dz         *ds   *dp 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is+1,   ip+1 )   *   dx           *dy         *dz         *ds   *dp) );
        }
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r )
    {
        gt_index_type ix = static_cast<gt_index_type>(std::floor(x));
        coord_type dx = x - ix;
        coord_type dx_prime = coord_type(1.0)-dx;

        gt_index_type iy = static_cast<gt_index_type>(std::floor(y));
        coord_type dy = y - iy;
        coord_type dy_prime = coord_type(1.0)-dy;

        gt_index_type iz = static_cast<gt_index_type>(std::floor(z));
        coord_type dz = z - iz;
        coord_type dz_prime = coord_type(1.0)-dz;

        gt_index_type is = static_cast<gt_index_type>(std::floor(s));
        coord_type ds = s - is;
        coord_type ds_prime = coord_type(1.0)-ds;

        gt_index_type ip = static_cast<gt_index_type>(std::floor(p));
        coord_type dp = p - ip;
        coord_type dp_prime = coord_type(1.0)-dp;

        gt_index_type ir = static_cast<gt_index_type>(std::floor(r));
        coord_type dr = r - ir;
        coord_type dr_prime = coord_type(1.0)-dr;

        if ( ix>=0 && ix<(gt_index_type)array_->get_size(0)-1 
            && iy>=0 && iy<(gt_index_type)array_->get_size(1)-1 
            && iz>=0 && iz<(gt_index_type)array_->get_size(2)-1 
            && is>=0 && is<(gt_index_type)array_->get_size(3)-1 
            && ip>=0 && ip<(gt_index_type)array_->get_size(4)-1 
            && ir>=0 && ir<(gt_index_type)array_->get_size(5)-1 )
        {
            return (    ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is),     size_t(ip),     size_t(ir) )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp_prime   *dr_prime
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is),     size_t(ip),     size_t(ir) )   *   dx           *dy_prime   *dz_prime   *ds_prime   *dp_prime   *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is),     size_t(ip),     size_t(ir) )   *   dx_prime     *dy         *dz_prime   *ds_prime   *dp_prime   *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is),     size_t(ip),     size_t(ir) )   *   dx           *dy         *dz_prime   *ds_prime   *dp_prime   *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is),     size_t(ip),     size_t(ir) )   *   dx_prime     *dy_prime   *dz         *ds_prime   *dp_prime   *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is),     size_t(ip),     size_t(ir) )   *   dx           *dy_prime   *dz         *ds_prime   *dp_prime   *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is),     size_t(ip),     size_t(ir) )   *   dx_prime     *dy         *dz         *ds_prime   *dp_prime   *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is),     size_t(ip),     size_t(ir) )   *   dx           *dy         *dz         *ds_prime   *dp_prime   *dr_prime)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is)+1,   size_t(ip),     size_t(ir) )   *   dx_prime     *dy_prime   *dz_prime   *ds         *dp_prime   *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is)+1,   size_t(ip),     size_t(ir) )   *   dx           *dy_prime   *dz_prime   *ds         *dp_prime   *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is)+1,   size_t(ip),     size_t(ir) )   *   dx_prime     *dy         *dz_prime   *ds         *dp_prime   *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is)+1,   size_t(ip),     size_t(ir) )   *   dx           *dy         *dz_prime   *ds         *dp_prime   *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is)+1,   size_t(ip),     size_t(ir) )   *   dx_prime     *dy_prime   *dz         *ds         *dp_prime   *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is)+1,   size_t(ip),     size_t(ir) )   *   dx           *dy_prime   *dz         *ds         *dp_prime   *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is)+1,   size_t(ip),     size_t(ir) )   *   dx_prime     *dy         *dz         *ds         *dp_prime   *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is)+1,   size_t(ip),     size_t(ir) )   *   dx           *dy         *dz         *ds         *dp_prime   *dr_prime)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is),     size_t(ip)+1,     size_t(ir) ) *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp         *dr_prime
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is),     size_t(ip)+1,     size_t(ir) ) *   dx           *dy_prime   *dz_prime   *ds_prime   *dp         *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is),     size_t(ip)+1,     size_t(ir) ) *   dx_prime     *dy         *dz_prime   *ds_prime   *dp         *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is),     size_t(ip)+1,     size_t(ir) ) *   dx           *dy         *dz_prime   *ds_prime   *dp         *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is),     size_t(ip)+1,     size_t(ir) ) *   dx_prime     *dy_prime   *dz         *ds_prime   *dp         *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is),     size_t(ip)+1,     size_t(ir) ) *   dx           *dy_prime   *dz         *ds_prime   *dp         *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is),     size_t(ip)+1,     size_t(ir) ) *   dx_prime     *dy         *dz         *ds_prime   *dp         *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is),     size_t(ip)+1,     size_t(ir) ) *   dx           *dy         *dz         *ds_prime   *dp         *dr_prime)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is)+1,   size_t(ip)+1,     size_t(ir) ) *   dx_prime     *dy_prime   *dz_prime   *ds         *dp         *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is)+1,   size_t(ip)+1,     size_t(ir) ) *   dx           *dy_prime   *dz_prime   *ds         *dp         *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is)+1,   size_t(ip)+1,     size_t(ir) ) *   dx_prime     *dy         *dz_prime   *ds         *dp         *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is)+1,   size_t(ip)+1,     size_t(ir) ) *   dx           *dy         *dz_prime   *ds         *dp         *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is)+1,   size_t(ip)+1,     size_t(ir) ) *   dx_prime     *dy_prime   *dz         *ds         *dp         *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is)+1,   size_t(ip)+1,     size_t(ir) ) *   dx           *dy_prime   *dz         *ds         *dp         *dr_prime) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is)+1,   size_t(ip)+1,     size_t(ir) ) *   dx_prime     *dy         *dz         *ds         *dp         *dr_prime 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is)+1,   size_t(ip)+1,     size_t(ir) ) *   dx           *dy         *dz         *ds         *dp         *dr_prime)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is),     size_t(ip),     size_t(ir)+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp_prime   *dr
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is),     size_t(ip),     size_t(ir)+1 )   *   dx           *dy_prime   *dz_prime   *ds_prime   *dp_prime   *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is),     size_t(ip),     size_t(ir)+1 )   *   dx_prime     *dy         *dz_prime   *ds_prime   *dp_prime   *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is),     size_t(ip),     size_t(ir)+1 )   *   dx           *dy         *dz_prime   *ds_prime   *dp_prime   *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is),     size_t(ip),     size_t(ir)+1 )   *   dx_prime     *dy_prime   *dz         *ds_prime   *dp_prime   *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is),     size_t(ip),     size_t(ir)+1 )   *   dx           *dy_prime   *dz         *ds_prime   *dp_prime   *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is),     size_t(ip),     size_t(ir)+1 )   *   dx_prime     *dy         *dz         *ds_prime   *dp_prime   *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is),     size_t(ip),     size_t(ir)+1 )   *   dx           *dy         *dz         *ds_prime   *dp_prime   *dr)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is)+1,   size_t(ip),     size_t(ir)+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds         *dp_prime   *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is)+1,   size_t(ip),     size_t(ir)+1 )   *   dx           *dy_prime   *dz_prime   *ds         *dp_prime   *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is)+1,   size_t(ip),     size_t(ir)+1 )   *   dx_prime     *dy         *dz_prime   *ds         *dp_prime   *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is)+1,   size_t(ip),     size_t(ir)+1 )   *   dx           *dy         *dz_prime   *ds         *dp_prime   *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is)+1,   size_t(ip),     size_t(ir)+1 )   *   dx_prime     *dy_prime   *dz         *ds         *dp_prime   *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is)+1,   size_t(ip),     size_t(ir)+1 )   *   dx           *dy_prime   *dz         *ds         *dp_prime   *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is)+1,   size_t(ip),     size_t(ir)+1 )   *   dx_prime     *dy         *dz         *ds         *dp_prime   *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is)+1,   size_t(ip),     size_t(ir)+1 )   *   dx           *dy         *dz         *ds         *dp_prime   *dr)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is),     size_t(ip)+1,     size_t(ir)+1 ) *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp         *dr
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is),     size_t(ip)+1,     size_t(ir)+1 ) *   dx           *dy_prime   *dz_prime   *ds_prime   *dp         *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is),     size_t(ip)+1,     size_t(ir)+1 ) *   dx_prime     *dy         *dz_prime   *ds_prime   *dp         *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is),     size_t(ip)+1,     size_t(ir)+1 ) *   dx           *dy         *dz_prime   *ds_prime   *dp         *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is),     size_t(ip)+1,     size_t(ir)+1 ) *   dx_prime     *dy_prime   *dz         *ds_prime   *dp         *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is),     size_t(ip)+1,     size_t(ir)+1 ) *   dx           *dy_prime   *dz         *ds_prime   *dp         *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is),     size_t(ip)+1,     size_t(ir)+1 ) *   dx_prime     *dy         *dz         *ds_prime   *dp         *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is),     size_t(ip)+1,     size_t(ir)+1 ) *   dx           *dy         *dz         *ds_prime   *dp         *dr)
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz),    size_t(is)+1,   size_t(ip)+1,     size_t(ir)+1 ) *   dx_prime     *dy_prime   *dz_prime   *ds         *dp         *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz),     size_t(is)+1,   size_t(ip)+1,     size_t(ir)+1 ) *   dx           *dy_prime   *dz_prime   *ds         *dp         *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz),    size_t(is)+1,   size_t(ip)+1,     size_t(ir)+1 ) *   dx_prime     *dy         *dz_prime   *ds         *dp         *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz),     size_t(is)+1,   size_t(ip)+1,     size_t(ir)+1 ) *   dx           *dy         *dz_prime   *ds         *dp         *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy),     size_t(iz)+1,  size_t(is)+1,   size_t(ip)+1,     size_t(ir)+1 ) *   dx_prime     *dy_prime   *dz         *ds         *dp         *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy),     size_t(iz)+1,   size_t(is)+1,   size_t(ip)+1,     size_t(ir)+1 ) *   dx           *dy_prime   *dz         *ds         *dp         *dr) 
                    +   ((*array_)(size_t(ix),   size_t(iy)+1,   size_t(iz)+1,  size_t(is)+1,   size_t(ip)+1,     size_t(ir)+1 ) *   dx_prime     *dy         *dz         *ds         *dp         *dr 
                    +   (*array_)(size_t(ix)+1, size_t(iy)+1,   size_t(iz)+1,   size_t(is)+1,   size_t(ip)+1,     size_t(ir)+1 ) *   dx           *dy         *dz         *ds         *dp         *dr) );
        }
        else
        {
            return (    ((*bh_)(ix,   iy,     iz,    is,     ip,     ir )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp_prime   *dr_prime
                    +   (*bh_)(ix+1, iy,     iz,     is,     ip,     ir )   *   dx           *dy_prime   *dz_prime   *ds_prime   *dp_prime   *dr_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is,     ip,     ir )   *   dx_prime     *dy         *dz_prime   *ds_prime   *dp_prime   *dr_prime 
                    +   (*bh_)(ix+1, iy+1,   iz,     is,     ip,     ir )   *   dx           *dy         *dz_prime   *ds_prime   *dp_prime   *dr_prime) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is,     ip,     ir )   *   dx_prime     *dy_prime   *dz         *ds_prime   *dp_prime   *dr_prime 
                    +   (*bh_)(ix+1, iy,     iz+1,   is,     ip,     ir )   *   dx           *dy_prime   *dz         *ds_prime   *dp_prime   *dr_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is,     ip,     ir )   *   dx_prime     *dy         *dz         *ds_prime   *dp_prime   *dr_prime 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is,     ip,     ir )   *   dx           *dy         *dz         *ds_prime   *dp_prime   *dr_prime)
                    +   ((*bh_)(ix,   iy,     iz,    is+1,   ip,     ir )   *   dx_prime     *dy_prime   *dz_prime   *ds         *dp_prime   *dr_prime 
                    +   (*bh_)(ix+1, iy,     iz,     is+1,   ip,     ir )   *   dx           *dy_prime   *dz_prime   *ds         *dp_prime   *dr_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is+1,   ip,     ir )   *   dx_prime     *dy         *dz_prime   *ds         *dp_prime   *dr_prime 
                    +   (*bh_)(ix+1, iy+1,   iz,     is+1,   ip,     ir )   *   dx           *dy         *dz_prime   *ds         *dp_prime   *dr_prime) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is+1,   ip,     ir )   *   dx_prime     *dy_prime   *dz         *ds         *dp_prime   *dr_prime 
                    +   (*bh_)(ix+1, iy,     iz+1,   is+1,   ip,     ir )   *   dx           *dy_prime   *dz         *ds         *dp_prime   *dr_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is+1,   ip,     ir )   *   dx_prime     *dy         *dz         *ds         *dp_prime   *dr_prime 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is+1,   ip,     ir )   *   dx           *dy         *dz         *ds         *dp_prime   *dr_prime)
                    +   ((*bh_)(ix,   iy,     iz,    is,     ip+1,     ir ) *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp         *dr_prime
                    +   (*bh_)(ix+1, iy,     iz,     is,     ip+1,     ir ) *   dx           *dy_prime   *dz_prime   *ds_prime   *dp         *dr_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is,     ip+1,     ir ) *   dx_prime     *dy         *dz_prime   *ds_prime   *dp         *dr_prime 
                    +   (*bh_)(ix+1, iy+1,   iz,     is,     ip+1,     ir ) *   dx           *dy         *dz_prime   *ds_prime   *dp         *dr_prime) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is,     ip+1,     ir ) *   dx_prime     *dy_prime   *dz         *ds_prime   *dp         *dr_prime 
                    +   (*bh_)(ix+1, iy,     iz+1,   is,     ip+1,     ir ) *   dx           *dy_prime   *dz         *ds_prime   *dp         *dr_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is,     ip+1,     ir ) *   dx_prime     *dy         *dz         *ds_prime   *dp         *dr_prime 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is,     ip+1,     ir ) *   dx           *dy         *dz         *ds_prime   *dp         *dr_prime)
                    +   ((*bh_)(ix,   iy,     iz,    is+1,   ip+1,     ir ) *   dx_prime     *dy_prime   *dz_prime   *ds         *dp         *dr_prime 
                    +   (*bh_)(ix+1, iy,     iz,     is+1,   ip+1,     ir ) *   dx           *dy_prime   *dz_prime   *ds         *dp         *dr_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is+1,   ip+1,     ir ) *   dx_prime     *dy         *dz_prime   *ds         *dp         *dr_prime 
                    +   (*bh_)(ix+1, iy+1,   iz,     is+1,   ip+1,     ir ) *   dx           *dy         *dz_prime   *ds         *dp         *dr_prime) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is+1,   ip+1,     ir ) *   dx_prime     *dy_prime   *dz         *ds         *dp         *dr_prime 
                    +   (*bh_)(ix+1, iy,     iz+1,   is+1,   ip+1,     ir ) *   dx           *dy_prime   *dz         *ds         *dp         *dr_prime) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is+1,   ip+1,     ir ) *   dx_prime     *dy         *dz         *ds         *dp         *dr_prime 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is+1,   ip+1,     ir ) *   dx           *dy         *dz         *ds         *dp         *dr_prime)
                    +   ((*bh_)(ix,   iy,     iz,    is,     ip,     ir+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp_prime   *dr
                    +   (*bh_)(ix+1, iy,     iz,     is,     ip,     ir+1 )   *   dx           *dy_prime   *dz_prime   *ds_prime   *dp_prime   *dr) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is,     ip,     ir+1 )   *   dx_prime     *dy         *dz_prime   *ds_prime   *dp_prime   *dr 
                    +   (*bh_)(ix+1, iy+1,   iz,     is,     ip,     ir+1 )   *   dx           *dy         *dz_prime   *ds_prime   *dp_prime   *dr) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is,     ip,     ir+1 )   *   dx_prime     *dy_prime   *dz         *ds_prime   *dp_prime   *dr 
                    +   (*bh_)(ix+1, iy,     iz+1,   is,     ip,     ir+1 )   *   dx           *dy_prime   *dz         *ds_prime   *dp_prime   *dr) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is,     ip,     ir+1 )   *   dx_prime     *dy         *dz         *ds_prime   *dp_prime   *dr 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is,     ip,     ir+1 )   *   dx           *dy         *dz         *ds_prime   *dp_prime   *dr)
                    +   ((*bh_)(ix,   iy,     iz,    is+1,   ip,     ir+1 )   *   dx_prime     *dy_prime   *dz_prime   *ds         *dp_prime   *dr 
                    +   (*bh_)(ix+1, iy,     iz,     is+1,   ip,     ir+1 )   *   dx           *dy_prime   *dz_prime   *ds         *dp_prime   *dr) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is+1,   ip,     ir+1 )   *   dx_prime     *dy         *dz_prime   *ds         *dp_prime   *dr 
                    +   (*bh_)(ix+1, iy+1,   iz,     is+1,   ip,     ir+1 )   *   dx           *dy         *dz_prime   *ds         *dp_prime   *dr) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is+1,   ip,     ir+1 )   *   dx_prime     *dy_prime   *dz         *ds         *dp_prime   *dr 
                    +   (*bh_)(ix+1, iy,     iz+1,   is+1,   ip,     ir+1 )   *   dx           *dy_prime   *dz         *ds         *dp_prime   *dr) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is+1,   ip,     ir+1 )   *   dx_prime     *dy         *dz         *ds         *dp_prime   *dr 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is+1,   ip,     ir+1 )   *   dx           *dy         *dz         *ds         *dp_prime   *dr)
                    +   ((*bh_)(ix,   iy,     iz,    is,     ip+1,     ir+1 ) *   dx_prime     *dy_prime   *dz_prime   *ds_prime   *dp         *dr
                    +   (*bh_)(ix+1, iy,     iz,     is,     ip+1,     ir+1 ) *   dx           *dy_prime   *dz_prime   *ds_prime   *dp         *dr) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is,     ip+1,     ir+1 ) *   dx_prime     *dy         *dz_prime   *ds_prime   *dp         *dr 
                    +   (*bh_)(ix+1, iy+1,   iz,     is,     ip+1,     ir+1 ) *   dx           *dy         *dz_prime   *ds_prime   *dp         *dr) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is,     ip+1,     ir+1 ) *   dx_prime     *dy_prime   *dz         *ds_prime   *dp         *dr 
                    +   (*bh_)(ix+1, iy,     iz+1,   is,     ip+1,     ir+1 ) *   dx           *dy_prime   *dz         *ds_prime   *dp         *dr) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is,     ip+1,     ir+1 ) *   dx_prime     *dy         *dz         *ds_prime   *dp         *dr 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is,     ip+1,     ir+1 ) *   dx           *dy         *dz         *ds_prime   *dp         *dr)
                    +   ((*bh_)(ix,   iy,     iz,    is+1,   ip+1,     ir+1 ) *   dx_prime     *dy_prime   *dz_prime   *ds         *dp         *dr 
                    +   (*bh_)(ix+1, iy,     iz,     is+1,   ip+1,     ir+1 ) *   dx           *dy_prime   *dz_prime   *ds         *dp         *dr) 
                    +   ((*bh_)(ix,   iy+1,   iz,    is+1,   ip+1,     ir+1 ) *   dx_prime     *dy         *dz_prime   *ds         *dp         *dr 
                    +   (*bh_)(ix+1, iy+1,   iz,     is+1,   ip+1,     ir+1 ) *   dx           *dy         *dz_prime   *ds         *dp         *dr) 
                    +   ((*bh_)(ix,   iy,     iz+1,  is+1,   ip+1,     ir+1 ) *   dx_prime     *dy_prime   *dz         *ds         *dp         *dr 
                    +   (*bh_)(ix+1, iy,     iz+1,   is+1,   ip+1,     ir+1 ) *   dx           *dy_prime   *dz         *ds         *dp         *dr) 
                    +   ((*bh_)(ix,   iy+1,   iz+1,  is+1,   ip+1,     ir+1 ) *   dx_prime     *dy         *dz         *ds         *dp         *dr 
                    +   (*bh_)(ix+1, iy+1,   iz+1,   is+1,   ip+1,     ir+1 ) *   dx           *dy         *dz         *ds         *dp         *dr) );
        }
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a )
    {
        gt_index_type anchor[7];
        coord_type d[7];
        coord_type d_prime[7];

        anchor[0] = static_cast<gt_index_type>(std::floor(x));
        anchor[1] = static_cast<gt_index_type>(std::floor(y));
        anchor[2] = static_cast<gt_index_type>(std::floor(z));
        anchor[3] = static_cast<gt_index_type>(std::floor(s));
        anchor[4] = static_cast<gt_index_type>(std::floor(p));
        anchor[5] = static_cast<gt_index_type>(std::floor(r));
        anchor[6] = static_cast<gt_index_type>(std::floor(a));

        d[0] = x - anchor[0];
        d[1] = y - anchor[1];
        d[2] = z - anchor[2];
        d[3] = s - anchor[3];
        d[4] = p - anchor[4];
        d[5] = r - anchor[5];
        d[6] = a - anchor[6];

        unsigned int ii;
        for ( ii=0; ii<7; ii++ )
        {
            d_prime[ii] = coord_type(1.0)-d[ii];
        }

        T res(0);

        coord_type weightAll(1.0);

        unsigned int n;

        if ( anchor[0]>=0 && anchor[0]<(gt_index_type)array_->get_size(0)-1 
            && anchor[1]>=0 && anchor[1]<(gt_index_type)array_->get_size(1)-1 
            && anchor[2]>=0 && anchor[2]<(gt_index_type)array_->get_size(2)-1 
            && anchor[3]>=0 && anchor[3]<(gt_index_type)array_->get_size(3)-1 
            && anchor[4]>=0 && anchor[4]<(gt_index_type)array_->get_size(4)-1 
            && anchor[5]>=0 && anchor[5]<(gt_index_type)array_->get_size(5)-1
            && anchor[6]>=0 && anchor[6]<(gt_index_type)array_->get_size(6)-1 )
        {
            std::vector<size_t> ind(7);

            for ( n=0; n<number_of_points_; n++ )
            {
                unsigned int lastDigit = n;
                weightAll = coord_type(1.0);

                for ( ii=0; ii<7; ii++ )
                {
                    if ( lastDigit & 1 )
                    {
                        ind[ii] = anchor[ii]+1;
                        weightAll *= d[ii];
                    }
                    else
                    {
                        ind[ii] = anchor[ii];
                        weightAll *= d_prime[ii];
                    }

                    // shift one digit
                    lastDigit >>= 1;
                }

                res += weightAll * (*array_)(ind);
            }
        }
        else
        {
            std::vector<gt_index_type> ind(7);

            for ( n=0; n<number_of_points_; n++ )
            {
                unsigned int lastDigit = n;
                weightAll = coord_type(1.0);

                for ( ii=0; ii<7; ii++ )
                {
                    if ( lastDigit & 1 )
                    {
                        ind[ii] = anchor[ii]+1;
                        weightAll *= d[ii];
                    }
                    else
                    {
                        ind[ii] = anchor[ii];
                        weightAll *= d_prime[ii];
                    }

                    // shift one digit
                    lastDigit >>= 1;
                }

                res += weightAll * (*bh_)(ind);
            }
        }

        return res;
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q )
    {
        gt_index_type anchor[8];
        coord_type d[8];
        coord_type d_prime[8];

        anchor[0] = static_cast<gt_index_type>(std::floor(x));
        anchor[1] = static_cast<gt_index_type>(std::floor(y));
        anchor[2] = static_cast<gt_index_type>(std::floor(z));
        anchor[3] = static_cast<gt_index_type>(std::floor(s));
        anchor[4] = static_cast<gt_index_type>(std::floor(p));
        anchor[5] = static_cast<gt_index_type>(std::floor(r));
        anchor[6] = static_cast<gt_index_type>(std::floor(a));
        anchor[7] = static_cast<gt_index_type>(std::floor(q));

        d[0] = x - anchor[0];
        d[1] = y - anchor[1];
        d[2] = z - anchor[2];
        d[3] = s - anchor[3];
        d[4] = p - anchor[4];
        d[5] = r - anchor[5];
        d[6] = a - anchor[6];
        d[7] = q - anchor[7];

        unsigned int ii;
        for ( ii=0; ii<8; ii++ )
        {
            d_prime[ii] = coord_type(1.0)-d[ii];
        }

        T res(0);

        coord_type weightAll(1.0);

        unsigned int n;

        if ( anchor[0]>=0 && anchor[0]<(gt_index_type)array_->get_size(0)-1 
            && anchor[1]>=0 && anchor[1]<(gt_index_type)array_->get_size(1)-1 
            && anchor[2]>=0 && anchor[2]<(gt_index_type)array_->get_size(2)-1 
            && anchor[3]>=0 && anchor[3]<(gt_index_type)array_->get_size(3)-1 
            && anchor[4]>=0 && anchor[4]<(gt_index_type)array_->get_size(4)-1 
            && anchor[5]>=0 && anchor[5]<(gt_index_type)array_->get_size(5)-1
            && anchor[6]>=0 && anchor[6]<(gt_index_type)array_->get_size(6)-1
            && anchor[7]>=0 && anchor[7]<(gt_index_type)array_->get_size(7)-1 )
        {
            std::vector<size_t> ind(8);

            for ( n=0; n<number_of_points_; n++ )
            {
                unsigned int lastDigit = n;
                weightAll = coord_type(1.0);

                for ( ii=0; ii<8; ii++ )
                {
                    if ( lastDigit & 1 )
                    {
                        ind[ii] = anchor[ii]+1;
                        weightAll *= d[ii];
                    }
                    else
                    {
                        ind[ii] = anchor[ii];
                        weightAll *= d_prime[ii];
                    }

                    // shift one digit
                    lastDigit >>= 1;
                }

                res += weightAll * (*array_)(ind);
            }
        }
        else
        {
            std::vector<gt_index_type> ind(8);

            for ( n=0; n<number_of_points_; n++ )
            {
                unsigned int lastDigit = n;
                weightAll = coord_type(1.0);

                for ( ii=0; ii<8; ii++ )
                {
                    if ( lastDigit & 1 )
                    {
                        ind[ii] = anchor[ii]+1;
                        weightAll *= d[ii];
                    }
                    else
                    {
                        ind[ii] = anchor[ii];
                        weightAll *= d_prime[ii];
                    }

                    // shift one digit
                    lastDigit >>= 1;
                }

                res += weightAll * (*bh_)(ind);
            }
        }

        return res;
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorLinear<ArrayType>::T hoNDInterpolatorLinear<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type u )
    {
        gt_index_type anchor[9];
        coord_type d[9];
        coord_type d_prime[9];

        anchor[0] = static_cast<gt_index_type>(std::floor(x));
        anchor[1] = static_cast<gt_index_type>(std::floor(y));
        anchor[2] = static_cast<gt_index_type>(std::floor(z));
        anchor[3] = static_cast<gt_index_type>(std::floor(s));
        anchor[4] = static_cast<gt_index_type>(std::floor(p));
        anchor[5] = static_cast<gt_index_type>(std::floor(r));
        anchor[6] = static_cast<gt_index_type>(std::floor(a));
        anchor[7] = static_cast<gt_index_type>(std::floor(q));
        anchor[8] = static_cast<gt_index_type>(std::floor(u));

        d[0] = x - anchor[0];
        d[1] = y - anchor[1];
        d[2] = z - anchor[2];
        d[3] = s - anchor[3];
        d[4] = p - anchor[4];
        d[5] = r - anchor[5];
        d[6] = a - anchor[6];
        d[7] = q - anchor[7];
        d[8] = u - anchor[8];

        unsigned int ii;
        for ( ii=0; ii<9; ii++ )
        {
            d_prime[ii] = coord_type(1.0)-d[ii];
        }

        T res(0);

        coord_type weightAll(1.0);

        unsigned int n;

        if ( anchor[0]>=0 && anchor[0]<(gt_index_type)array_->get_size(0)-1 
            && anchor[1]>=0 && anchor[1]<(gt_index_type)array_->get_size(1)-1 
            && anchor[2]>=0 && anchor[2]<(gt_index_type)array_->get_size(2)-1 
            && anchor[3]>=0 && anchor[3]<(gt_index_type)array_->get_size(3)-1 
            && anchor[4]>=0 && anchor[4]<(gt_index_type)array_->get_size(4)-1 
            && anchor[5]>=0 && anchor[5]<(gt_index_type)array_->get_size(5)-1
            && anchor[6]>=0 && anchor[6]<(gt_index_type)array_->get_size(6)-1
            && anchor[7]>=0 && anchor[7]<(gt_index_type)array_->get_size(7)-1
            && anchor[8]>=0 && anchor[8]<(gt_index_type)array_->get_size(8)-1 )
        {
            std::vector<size_t> ind(9);

            for ( n=0; n<number_of_points_; n++ )
            {
                unsigned int lastDigit = n;
                weightAll = coord_type(1.0);

                for ( ii=0; ii<9; ii++ )
                {
                    if ( lastDigit & 1 )
                    {
                        ind[ii] = anchor[ii]+1;
                        weightAll *= d[ii];
                    }
                    else
                    {
                        ind[ii] = anchor[ii];
                        weightAll *= d_prime[ii];
                    }

                    // shift one digit
                    lastDigit >>= 1;
                }

                res += weightAll * (*array_)(ind);
            }
        }
        else
        {
            std::vector<gt_index_type> ind(9);

            for ( n=0; n<number_of_points_; n++ )
            {
                unsigned int lastDigit = n;
                weightAll = coord_type(1.0);

                for ( ii=0; ii<9; ii++ )
                {
                    if ( lastDigit & 1 )
                    {
                        ind[ii] = anchor[ii]+1;
                        weightAll *= d[ii];
                    }
                    else
                    {
                        ind[ii] = anchor[ii];
                        weightAll *= d_prime[ii];
                    }

                    // shift one digit
                    lastDigit >>= 1;
                }

                res += weightAll * (*bh_)(ind);
            }
        }

        return res;
    }
}
