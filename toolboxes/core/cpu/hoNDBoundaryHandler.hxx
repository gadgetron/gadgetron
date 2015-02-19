/** \file       hoNDBoundaryHandler.hxx
    \brief      N-dimensional boundary condition handler

                Designed to work with hoNDArray and hoNDImage

    \author     Hui Xue
*/

namespace Gadgetron
{
    /// hoNDBoundaryHandlerFixedValue

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( const std::vector<long long>& ind )
    {
        size_t D = ind.size();

        bool inside = true;

        size_t d;
        for ( d=0; d<D; d++)
        {
            if ( (ind[d]<0) || (ind[d]>=array_->get_size(d)) )
            {
                inside = false;
                break;
            }
        }

        if (inside)
        {
            size_t offset = ind[0];
            for ( d=1; d<D; d++ )
            {
                offset += ind[d]*array_->get_offset_factor(d);
            }

            return (*array_)(offset);
        }
        else
        {
            return value_;
        }
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( long long x )
    {
        return ((x >= 0) && array_->point_in_range( size_t(x) )) ? (*array_)( size_t(x) ) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( long long x, long long y )
    {
        return ((x >=0 ) && (y >= 0) && array_->point_in_range(size_t(x), size_t(y) )) ? (*array_)(size_t(x), size_t(y)) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( long long x, long long y, long long z )
    {
        return ((x >= 0) && (y >= 0) && (z >= 0) && array_->point_in_range(size_t(x), size_t(y), size_t(z))) ? (*array_)(size_t(x), size_t(y), size_t(z)) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( long long x, long long y, long long z, long long s )
    {
        return ((x >= 0) && (y >= 0) && (z >= 0) && (s >= 0) && array_->point_in_range(size_t(x), size_t(y), size_t(z), size_t(s))) ? (*array_)(size_t(x), size_t(y), size_t(z), size_t(s)) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p )
    {
        return ((x >= 0) && (y >= 0) && (z >= 0) && (s >= 0) && (p >= 0) && array_->point_in_range(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p))) ? (*array_)(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p)) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r )
    {
        return ((x >= 0) && (y >= 0) && (z >= 0) && (s >= 0) && (p >= 0) && (r >= 0) && array_->point_in_range(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r))) ? (*array_)(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r)) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a )
    {
        return ((x >= 0) && (y >= 0) && (z >= 0) && (s >= 0) && (p >= 0) && (r >= 0) && (a >= 0) && array_->point_in_range(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r), size_t(a))) ? (*array_)(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r), size_t(a)) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q )
    {
        return ((x >= 0) && (y >= 0) && (z >= 0) && (s >= 0) && (p >= 0) && (r >= 0) && (a >= 0) && (q >= 0) && array_->point_in_range(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r), size_t(a), size_t(q))) ? (*array_)(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r), size_t(a), size_t(q)) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q, long long u )
    {
        return ((x >= 0) && (y >= 0) && (z >= 0) && (s >= 0) && (p >= 0) && (r >= 0) && (a >= 0) && (q >= 0) && (u >= 0) && array_->point_in_range(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r), size_t(a), size_t(q), size_t(u))) ? (*array_)(size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r), size_t(a), size_t(q), size_t(u)) : value_;
    }

    /// hoNDBoundaryHandlerBorderValue

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( const std::vector<long long>& ind )
    {
        std::vector<size_t> indInside(array_->get_number_of_dimensions());
        unsigned int D = array_->get_number_of_dimensions();
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            if ( ind[ii] < 0 )
            {
                indInside[ii] = 0;
            }
            else if ( ind[ii] >= array_->get_size(ii) )
            {
                indInside[ii] = array_->get_size(ii)-1;
            }
            else
            {
                indInside[ii] = ind[ii];
            }
        }
        return (*array_)(indInside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( long long x )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(long long)array_->get_size(0)) ? array_->get_size(0)-1 : x );

        return (*array_)(x_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( long long x, long long y )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(long long)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(long long)array_->get_size(1)) ? array_->get_size(1)-1 : y );

        return (*array_)(x_inside, y_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( long long x, long long y, long long z )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(long long)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(long long)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(long long)array_->get_size(2)) ? array_->get_size(2)-1 : z );

        return (*array_)(x_inside, y_inside, z_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( long long x, long long y, long long z, long long s )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(long long)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(long long)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(long long)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(long long)array_->get_size(3)) ? array_->get_size(3)-1 : s );

        return (*array_)(x_inside, y_inside, z_inside, s_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(long long)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(long long)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(long long)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(long long)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(long long)array_->get_size(4)) ? array_->get_size(4)-1 : p );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(long long)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(long long)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(long long)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(long long)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(long long)array_->get_size(4)) ? array_->get_size(4)-1 : p );
        size_t r_inside = (r<0) ? 0 : ( (r>=(long long)array_->get_size(5)) ? array_->get_size(5)-1 : r );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(long long)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(long long)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(long long)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(long long)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(long long)array_->get_size(4)) ? array_->get_size(4)-1 : p );
        size_t r_inside = (r<0) ? 0 : ( (r>=(long long)array_->get_size(5)) ? array_->get_size(5)-1 : r );
        size_t a_inside = (a<0) ? 0 : ( (a>=(long long)array_->get_size(6)) ? array_->get_size(6)-1 : a );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(long long)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(long long)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(long long)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(long long)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(long long)array_->get_size(4)) ? array_->get_size(4)-1 : p );
        size_t r_inside = (r<0) ? 0 : ( (r>=(long long)array_->get_size(5)) ? array_->get_size(5)-1 : r );
        size_t a_inside = (a<0) ? 0 : ( (a>=(long long)array_->get_size(6)) ? array_->get_size(6)-1 : a );
        size_t q_inside = (q<0) ? 0 : ( (q>=(long long)array_->get_size(7)) ? array_->get_size(7)-1 : q );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q, long long u )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(long long)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(long long)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(long long)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(long long)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(long long)array_->get_size(4)) ? array_->get_size(4)-1 : p );
        size_t r_inside = (r<0) ? 0 : ( (r>=(long long)array_->get_size(5)) ? array_->get_size(5)-1 : r );
        size_t a_inside = (a<0) ? 0 : ( (a>=(long long)array_->get_size(6)) ? array_->get_size(6)-1 : a );
        size_t q_inside = (q<0) ? 0 : ( (q>=(long long)array_->get_size(7)) ? array_->get_size(7)-1 : q );
        size_t u_inside = (u<0) ? 0 : ( (u>=(long long)array_->get_size(8)) ? array_->get_size(8)-1 : u );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside, u_inside);
    }

    /// hoNDBoundaryHandlerPeriodic

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( const std::vector<long long>& ind )
    {
        unsigned int D = (unsigned int)array_->get_number_of_dimensions();
        std::vector<size_t> indInside(D);
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            if ( (ind[ii]<0) || (ind[ii]>=(long long)array_->get_size(ii)) )
            {
                indInside[ii] = this->mod(ind[ii], array_->get_size(ii));
            }
            else
            {
                indInside[ii] = ind[ii];
            }
        }
        return (*array_)(indInside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( long long x )
    {
        size_t x_inside = (x<0 || x>=(long long)array_->get_size(0)) ? (this->mod(x, (long long)array_->get_size(0))) : x;

        return (*array_)(x_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( long long x, long long y )
    {
        size_t x_inside = (x<0 || x>=(long long)array_->get_size(0)) ? (this->mod(x, (long long)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(long long)array_->get_size(1)) ? (this->mod(y, (long long)array_->get_size(1))) : y;

        return (*array_)(x_inside, y_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( long long x, long long y, long long z )
    {
        size_t x_inside = (x<0 || x>=(long long)array_->get_size(0)) ? (this->mod(x, (long long)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(long long)array_->get_size(1)) ? (this->mod(y, (long long)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(long long)array_->get_size(2)) ? (this->mod(z, (long long)array_->get_size(2))) : z;

        return (*array_)(x_inside, y_inside, z_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( long long x, long long y, long long z, long long s )
    {
        size_t x_inside = (x<0 || x>=(long long)array_->get_size(0)) ? (this->mod(x, (long long)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(long long)array_->get_size(1)) ? (this->mod(y, (long long)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(long long)array_->get_size(2)) ? (this->mod(z, (long long)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(long long)array_->get_size(3)) ? (this->mod(s, (long long)array_->get_size(3))) : s;

        return (*array_)(x_inside, y_inside, z_inside, s_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p )
    {
        size_t x_inside = (x<0 || x>=(long long)array_->get_size(0)) ? (this->mod(x, (long long)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(long long)array_->get_size(1)) ? (this->mod(y, (long long)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(long long)array_->get_size(2)) ? (this->mod(z, (long long)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(long long)array_->get_size(3)) ? (this->mod(s, (long long)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(long long)array_->get_size(4)) ? (this->mod(p, (long long)array_->get_size(4))) : p;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r )
    {
        size_t x_inside = (x<0 || x>=(long long)array_->get_size(0)) ? (this->mod(x, (long long)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(long long)array_->get_size(1)) ? (this->mod(y, (long long)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(long long)array_->get_size(2)) ? (this->mod(z, (long long)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(long long)array_->get_size(3)) ? (this->mod(s, (long long)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(long long)array_->get_size(4)) ? (this->mod(p, (long long)array_->get_size(4))) : p;
        size_t r_inside = (r<0 || r>=(long long)array_->get_size(5)) ? (this->mod(r, (long long)array_->get_size(5))) : r;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a )
    {
        size_t x_inside = (x<0 || x>=(long long)array_->get_size(0)) ? (this->mod(x, (long long)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(long long)array_->get_size(1)) ? (this->mod(y, (long long)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(long long)array_->get_size(2)) ? (this->mod(z, (long long)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(long long)array_->get_size(3)) ? (this->mod(s, (long long)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(long long)array_->get_size(4)) ? (this->mod(p, (long long)array_->get_size(4))) : p;
        size_t r_inside = (r<0 || r>=(long long)array_->get_size(5)) ? (this->mod(r, (long long)array_->get_size(5))) : r;
        size_t a_inside = (a<0 || a>=(long long)array_->get_size(6)) ? (this->mod(a, (long long)array_->get_size(6))) : a;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q )
    {
        size_t x_inside = (x<0 || x>=(long long)array_->get_size(0)) ? (this->mod(x, (long long)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(long long)array_->get_size(1)) ? (this->mod(y, (long long)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(long long)array_->get_size(2)) ? (this->mod(z, (long long)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(long long)array_->get_size(3)) ? (this->mod(s, (long long)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(long long)array_->get_size(4)) ? (this->mod(p, (long long)array_->get_size(4))) : p;
        size_t r_inside = (r<0 || r>=(long long)array_->get_size(5)) ? (this->mod(r, (long long)array_->get_size(5))) : r;
        size_t a_inside = (a<0 || a>=(long long)array_->get_size(6)) ? (this->mod(a, (long long)array_->get_size(6))) : a;
        size_t q_inside = (q<0 || q>=(long long)array_->get_size(7)) ? (this->mod(q, (long long)array_->get_size(7))) : q;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q, long long u )
    {
        size_t x_inside = (x<0 || x>=(long long)array_->get_size(0)) ? (this->mod(x, (long long)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(long long)array_->get_size(1)) ? (this->mod(y, (long long)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(long long)array_->get_size(2)) ? (this->mod(z, (long long)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(long long)array_->get_size(3)) ? (this->mod(s, (long long)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(long long)array_->get_size(4)) ? (this->mod(p, (long long)array_->get_size(4))) : p;
        size_t r_inside = (r<0 || r>=(long long)array_->get_size(5)) ? (this->mod(r, (long long)array_->get_size(5))) : r;
        size_t a_inside = (a<0 || a>=(long long)array_->get_size(6)) ? (this->mod(a, (long long)array_->get_size(6))) : a;
        size_t q_inside = (q<0 || q>=(long long)array_->get_size(7)) ? (this->mod(q, (long long)array_->get_size(7))) : q;
        size_t u_inside = (u<0 || u>=(long long)array_->get_size(8)) ? (this->mod(u, (long long)array_->get_size(8))) : u;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside, u_inside);
    }

    /// hoNDBoundaryHandlerMirror

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( const std::vector<long long>& ind )
    {
        unsigned int D = array_->get_number_of_dimensions();
        std::vector<size_t> indInside(D);
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            if ( ind[ii] < 0 )
            {
                indInside[ii] = -ind[ii];
            }
            else if ( ind[ii] >= (long long)array_->get_size(ii) )
            {
                indInside[ii] = 2*(long long)array_->get_size(ii) - ind[ii] -2;
            }
            else
            {
                indInside[ii] = ind[ii];
            }
        }
        return (*array_)(indInside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( long long x )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(long long)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );

        return (*array_)(x_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( long long x, long long y )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(long long)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(long long)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );

        return (*array_)(x_inside, y_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( long long x, long long y, long long z )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(long long)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(long long)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(long long)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );

        return (*array_)(x_inside, y_inside, z_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( long long x, long long y, long long z, long long s )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(long long)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(long long)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(long long)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(long long)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );

        return (*array_)(x_inside, y_inside, z_inside, s_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(long long)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(long long)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(long long)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(long long)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(long long)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(long long)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(long long)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(long long)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(long long)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(long long)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );
        size_t r_inside = (r<0) ? -r : ( (r>=(long long)array_->get_size(5)) ? (2*array_->get_size(5)-r-2) : r );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(long long)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(long long)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(long long)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(long long)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(long long)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );
        size_t r_inside = (r<0) ? -r : ( (r>=(long long)array_->get_size(5)) ? (2*array_->get_size(5)-r-2) : r );
        size_t a_inside = (a<0) ? -a : ( (a>=(long long)array_->get_size(6)) ? (2*array_->get_size(6)-a-2) : a );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(long long)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(long long)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(long long)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(long long)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(long long)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );
        size_t r_inside = (r<0) ? -r : ( (r>=(long long)array_->get_size(5)) ? (2*array_->get_size(5)-r-2) : r );
        size_t a_inside = (a<0) ? -a : ( (a>=(long long)array_->get_size(6)) ? (2*array_->get_size(6)-a-2) : a );
        size_t q_inside = (q<0) ? -q : ( (q>=(long long)array_->get_size(7)) ? (2*array_->get_size(7)-q-2) : q );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( long long x, long long y, long long z, long long s, long long p, long long r, long long a, long long q, long long u )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(long long)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(long long)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(long long)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(long long)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(long long)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );
        size_t r_inside = (r<0) ? -r : ( (r>=(long long)array_->get_size(5)) ? (2*array_->get_size(5)-r-2) : r );
        size_t a_inside = (a<0) ? -a : ( (a>=(long long)array_->get_size(6)) ? (2*array_->get_size(6)-a-2) : a );
        size_t q_inside = (q<0) ? -q : ( (q>=(long long)array_->get_size(7)) ? (2*array_->get_size(7)-q-2) : q );
        size_t u_inside = (u<0) ? -u : ( (u>=(long long)array_->get_size(8)) ? (2*array_->get_size(8)-u-2) : u );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside, u_inside);
    }
}
