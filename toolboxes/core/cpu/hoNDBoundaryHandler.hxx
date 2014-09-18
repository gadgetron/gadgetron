/** \file       hoNDBoundaryHandler.hxx
    \brief      N-dimensional boundary condition handler

                Designed to work with hoNDArray and hoNDImage

    \author     Hui Xue
*/

namespace Gadgetron
{
    /// hoNDBoundaryHandlerFixedValue

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( const std::vector<gt_index_type>& ind )
    {
        return (array_->point_in_range(ind)) ? (*array_)(ind) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( gt_index_type x )
    {
        return (array_->point_in_range(x)) ? (*array_)( size_t(x) ) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( gt_index_type x, gt_index_type y )
    {
        return (array_->point_in_range(x, y)) ? (*array_)( size_t(x), size_t(y) ) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z )
    {
        return (array_->point_in_range(x, y, z)) ? (*array_)( size_t(x), size_t(y), size_t(z) ) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s )
    {
        return (array_->point_in_range(x, y, z, s)) ? (*array_)( size_t(x), size_t(y), size_t(z), size_t(s) ) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p )
    {
        return (array_->point_in_range(x, y, z, s, p)) ? (*array_)( size_t(x), size_t(y), size_t(z), size_t(s), size_t(p) ) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r )
    {
        return (array_->point_in_range(x, y, z, s, p, r)) ? (*array_)( size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r) ) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a )
    {
        return (array_->point_in_range(x, y, z, s, p, r, a)) ? (*array_)( size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r), size_t(a) ) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a, gt_index_type q )
    {
        return (array_->point_in_range(x, y, z, s, p, r, a, q)) ? (*array_)( size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r), size_t(a), size_t(q) ) : value_;
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerFixedValue<ArrayType>::T hoNDBoundaryHandlerFixedValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a, gt_index_type q, gt_index_type u )
    {
        return (array_->point_in_range(x, y, z, s, p, r, a, q, u)) ? (*array_)( size_t(x), size_t(y), size_t(z), size_t(s), size_t(p), size_t(r), size_t(a), size_t(q), size_t(u) ) : value_;
    }

    /// hoNDBoundaryHandlerBorderValue

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( const std::vector<gt_index_type>& ind )
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
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( gt_index_type x )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(gt_index_type)array_->get_size(0)) ? array_->get_size(0)-1 : x );

        return (*array_)(x_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( gt_index_type x, gt_index_type y )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(gt_index_type)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(gt_index_type)array_->get_size(1)) ? array_->get_size(1)-1 : y );

        return (*array_)(x_inside, y_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(gt_index_type)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(gt_index_type)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(gt_index_type)array_->get_size(2)) ? array_->get_size(2)-1 : z );

        return (*array_)(x_inside, y_inside, z_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(gt_index_type)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(gt_index_type)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(gt_index_type)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(gt_index_type)array_->get_size(3)) ? array_->get_size(3)-1 : s );

        return (*array_)(x_inside, y_inside, z_inside, s_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(gt_index_type)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(gt_index_type)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(gt_index_type)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(gt_index_type)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(gt_index_type)array_->get_size(4)) ? array_->get_size(4)-1 : p );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(gt_index_type)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(gt_index_type)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(gt_index_type)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(gt_index_type)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(gt_index_type)array_->get_size(4)) ? array_->get_size(4)-1 : p );
        size_t r_inside = (r<0) ? 0 : ( (r>=(gt_index_type)array_->get_size(5)) ? array_->get_size(5)-1 : r );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(gt_index_type)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(gt_index_type)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(gt_index_type)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(gt_index_type)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(gt_index_type)array_->get_size(4)) ? array_->get_size(4)-1 : p );
        size_t r_inside = (r<0) ? 0 : ( (r>=(gt_index_type)array_->get_size(5)) ? array_->get_size(5)-1 : r );
        size_t a_inside = (a<0) ? 0 : ( (a>=(gt_index_type)array_->get_size(6)) ? array_->get_size(6)-1 : a );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a, gt_index_type q )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(gt_index_type)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(gt_index_type)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(gt_index_type)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(gt_index_type)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(gt_index_type)array_->get_size(4)) ? array_->get_size(4)-1 : p );
        size_t r_inside = (r<0) ? 0 : ( (r>=(gt_index_type)array_->get_size(5)) ? array_->get_size(5)-1 : r );
        size_t a_inside = (a<0) ? 0 : ( (a>=(gt_index_type)array_->get_size(6)) ? array_->get_size(6)-1 : a );
        size_t q_inside = (q<0) ? 0 : ( (q>=(gt_index_type)array_->get_size(7)) ? array_->get_size(7)-1 : q );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerBorderValue<ArrayType>::T hoNDBoundaryHandlerBorderValue<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a, gt_index_type q, gt_index_type u )
    {
        size_t x_inside = (x<0) ? 0 : ( (x>=(gt_index_type)array_->get_size(0)) ? array_->get_size(0)-1 : x );
        size_t y_inside = (y<0) ? 0 : ( (y>=(gt_index_type)array_->get_size(1)) ? array_->get_size(1)-1 : y );
        size_t z_inside = (z<0) ? 0 : ( (z>=(gt_index_type)array_->get_size(2)) ? array_->get_size(2)-1 : z );
        size_t s_inside = (s<0) ? 0 : ( (s>=(gt_index_type)array_->get_size(3)) ? array_->get_size(3)-1 : s );
        size_t p_inside = (p<0) ? 0 : ( (p>=(gt_index_type)array_->get_size(4)) ? array_->get_size(4)-1 : p );
        size_t r_inside = (r<0) ? 0 : ( (r>=(gt_index_type)array_->get_size(5)) ? array_->get_size(5)-1 : r );
        size_t a_inside = (a<0) ? 0 : ( (a>=(gt_index_type)array_->get_size(6)) ? array_->get_size(6)-1 : a );
        size_t q_inside = (q<0) ? 0 : ( (q>=(gt_index_type)array_->get_size(7)) ? array_->get_size(7)-1 : q );
        size_t u_inside = (u<0) ? 0 : ( (u>=(gt_index_type)array_->get_size(8)) ? array_->get_size(8)-1 : u );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside, u_inside);
    }

    /// hoNDBoundaryHandlerPeriodic

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( const std::vector<gt_index_type>& ind )
    {
        unsigned int D = (unsigned int)array_->get_number_of_dimensions();
        std::vector<size_t> indInside(D);
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            if ( (ind[ii]<0) || (ind[ii]>=(gt_index_type)array_->get_size(ii)) )
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
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( gt_index_type x )
    {
        size_t x_inside = (x<0 || x>=(gt_index_type)array_->get_size(0)) ? (this->mod(x, (gt_index_type)array_->get_size(0))) : x;

        return (*array_)(x_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( gt_index_type x, gt_index_type y )
    {
        size_t x_inside = (x<0 || x>=(gt_index_type)array_->get_size(0)) ? (this->mod(x, (gt_index_type)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(gt_index_type)array_->get_size(1)) ? (this->mod(y, (gt_index_type)array_->get_size(1))) : y;

        return (*array_)(x_inside, y_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z )
    {
        size_t x_inside = (x<0 || x>=(gt_index_type)array_->get_size(0)) ? (this->mod(x, (gt_index_type)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(gt_index_type)array_->get_size(1)) ? (this->mod(y, (gt_index_type)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(gt_index_type)array_->get_size(2)) ? (this->mod(z, (gt_index_type)array_->get_size(2))) : z;

        return (*array_)(x_inside, y_inside, z_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s )
    {
        size_t x_inside = (x<0 || x>=(gt_index_type)array_->get_size(0)) ? (this->mod(x, (gt_index_type)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(gt_index_type)array_->get_size(1)) ? (this->mod(y, (gt_index_type)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(gt_index_type)array_->get_size(2)) ? (this->mod(z, (gt_index_type)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(gt_index_type)array_->get_size(3)) ? (this->mod(s, (gt_index_type)array_->get_size(3))) : s;

        return (*array_)(x_inside, y_inside, z_inside, s_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p )
    {
        size_t x_inside = (x<0 || x>=(gt_index_type)array_->get_size(0)) ? (this->mod(x, (gt_index_type)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(gt_index_type)array_->get_size(1)) ? (this->mod(y, (gt_index_type)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(gt_index_type)array_->get_size(2)) ? (this->mod(z, (gt_index_type)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(gt_index_type)array_->get_size(3)) ? (this->mod(s, (gt_index_type)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(gt_index_type)array_->get_size(4)) ? (this->mod(p, (gt_index_type)array_->get_size(4))) : p;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r )
    {
        size_t x_inside = (x<0 || x>=(gt_index_type)array_->get_size(0)) ? (this->mod(x, (gt_index_type)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(gt_index_type)array_->get_size(1)) ? (this->mod(y, (gt_index_type)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(gt_index_type)array_->get_size(2)) ? (this->mod(z, (gt_index_type)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(gt_index_type)array_->get_size(3)) ? (this->mod(s, (gt_index_type)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(gt_index_type)array_->get_size(4)) ? (this->mod(p, (gt_index_type)array_->get_size(4))) : p;
        size_t r_inside = (r<0 || r>=(gt_index_type)array_->get_size(5)) ? (this->mod(r, (gt_index_type)array_->get_size(5))) : r;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a )
    {
        size_t x_inside = (x<0 || x>=(gt_index_type)array_->get_size(0)) ? (this->mod(x, (gt_index_type)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(gt_index_type)array_->get_size(1)) ? (this->mod(y, (gt_index_type)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(gt_index_type)array_->get_size(2)) ? (this->mod(z, (gt_index_type)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(gt_index_type)array_->get_size(3)) ? (this->mod(s, (gt_index_type)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(gt_index_type)array_->get_size(4)) ? (this->mod(p, (gt_index_type)array_->get_size(4))) : p;
        size_t r_inside = (r<0 || r>=(gt_index_type)array_->get_size(5)) ? (this->mod(r, (gt_index_type)array_->get_size(5))) : r;
        size_t a_inside = (a<0 || a>=(gt_index_type)array_->get_size(6)) ? (this->mod(a, (gt_index_type)array_->get_size(6))) : a;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a, gt_index_type q )
    {
        size_t x_inside = (x<0 || x>=(gt_index_type)array_->get_size(0)) ? (this->mod(x, (gt_index_type)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(gt_index_type)array_->get_size(1)) ? (this->mod(y, (gt_index_type)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(gt_index_type)array_->get_size(2)) ? (this->mod(z, (gt_index_type)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(gt_index_type)array_->get_size(3)) ? (this->mod(s, (gt_index_type)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(gt_index_type)array_->get_size(4)) ? (this->mod(p, (gt_index_type)array_->get_size(4))) : p;
        size_t r_inside = (r<0 || r>=(gt_index_type)array_->get_size(5)) ? (this->mod(r, (gt_index_type)array_->get_size(5))) : r;
        size_t a_inside = (a<0 || a>=(gt_index_type)array_->get_size(6)) ? (this->mod(a, (gt_index_type)array_->get_size(6))) : a;
        size_t q_inside = (q<0 || q>=(gt_index_type)array_->get_size(7)) ? (this->mod(q, (gt_index_type)array_->get_size(7))) : q;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerPeriodic<ArrayType>::T hoNDBoundaryHandlerPeriodic<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a, gt_index_type q, gt_index_type u )
    {
        size_t x_inside = (x<0 || x>=(gt_index_type)array_->get_size(0)) ? (this->mod(x, (gt_index_type)array_->get_size(0))) : x;
        size_t y_inside = (y<0 || y>=(gt_index_type)array_->get_size(1)) ? (this->mod(y, (gt_index_type)array_->get_size(1))) : y;
        size_t z_inside = (z<0 || z>=(gt_index_type)array_->get_size(2)) ? (this->mod(z, (gt_index_type)array_->get_size(2))) : z;
        size_t s_inside = (s<0 || s>=(gt_index_type)array_->get_size(3)) ? (this->mod(s, (gt_index_type)array_->get_size(3))) : s;
        size_t p_inside = (p<0 || p>=(gt_index_type)array_->get_size(4)) ? (this->mod(p, (gt_index_type)array_->get_size(4))) : p;
        size_t r_inside = (r<0 || r>=(gt_index_type)array_->get_size(5)) ? (this->mod(r, (gt_index_type)array_->get_size(5))) : r;
        size_t a_inside = (a<0 || a>=(gt_index_type)array_->get_size(6)) ? (this->mod(a, (gt_index_type)array_->get_size(6))) : a;
        size_t q_inside = (q<0 || q>=(gt_index_type)array_->get_size(7)) ? (this->mod(q, (gt_index_type)array_->get_size(7))) : q;
        size_t u_inside = (u<0 || u>=(gt_index_type)array_->get_size(8)) ? (this->mod(u, (gt_index_type)array_->get_size(8))) : u;

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside, u_inside);
    }

    /// hoNDBoundaryHandlerMirror

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( const std::vector<gt_index_type>& ind )
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
            else if ( ind[ii] >= (gt_index_type)array_->get_size(ii) )
            {
                indInside[ii] = 2*(gt_index_type)array_->get_size(ii) - ind[ii] -2;
            }
            else
            {
                indInside[ii] = ind[ii];
            }
        }
        return (*array_)(indInside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( gt_index_type x )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(gt_index_type)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );

        return (*array_)(x_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( gt_index_type x, gt_index_type y )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(gt_index_type)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(gt_index_type)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );

        return (*array_)(x_inside, y_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(gt_index_type)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(gt_index_type)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(gt_index_type)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );

        return (*array_)(x_inside, y_inside, z_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(gt_index_type)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(gt_index_type)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(gt_index_type)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(gt_index_type)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );

        return (*array_)(x_inside, y_inside, z_inside, s_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(gt_index_type)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(gt_index_type)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(gt_index_type)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(gt_index_type)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(gt_index_type)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(gt_index_type)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(gt_index_type)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(gt_index_type)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(gt_index_type)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(gt_index_type)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );
        size_t r_inside = (r<0) ? -r : ( (r>=(gt_index_type)array_->get_size(5)) ? (2*array_->get_size(5)-r-2) : r );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(gt_index_type)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(gt_index_type)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(gt_index_type)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(gt_index_type)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(gt_index_type)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );
        size_t r_inside = (r<0) ? -r : ( (r>=(gt_index_type)array_->get_size(5)) ? (2*array_->get_size(5)-r-2) : r );
        size_t a_inside = (a<0) ? -a : ( (a>=(gt_index_type)array_->get_size(6)) ? (2*array_->get_size(6)-a-2) : a );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a, gt_index_type q )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(gt_index_type)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(gt_index_type)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(gt_index_type)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(gt_index_type)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(gt_index_type)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );
        size_t r_inside = (r<0) ? -r : ( (r>=(gt_index_type)array_->get_size(5)) ? (2*array_->get_size(5)-r-2) : r );
        size_t a_inside = (a<0) ? -a : ( (a>=(gt_index_type)array_->get_size(6)) ? (2*array_->get_size(6)-a-2) : a );
        size_t q_inside = (q<0) ? -q : ( (q>=(gt_index_type)array_->get_size(7)) ? (2*array_->get_size(7)-q-2) : q );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside);
    }

    template <typename ArrayType> 
    inline typename hoNDBoundaryHandlerMirror<ArrayType>::T hoNDBoundaryHandlerMirror<ArrayType>::operator()( gt_index_type x, gt_index_type y, gt_index_type z, gt_index_type s, gt_index_type p, gt_index_type r, gt_index_type a, gt_index_type q, gt_index_type u )
    {
        size_t x_inside = (x<0) ? -x : ( (x>=(gt_index_type)array_->get_size(0)) ? (2*array_->get_size(0)-x-2) : x );
        size_t y_inside = (y<0) ? -y : ( (y>=(gt_index_type)array_->get_size(1)) ? (2*array_->get_size(1)-y-2) : y );
        size_t z_inside = (z<0) ? -z : ( (z>=(gt_index_type)array_->get_size(2)) ? (2*array_->get_size(2)-z-2) : z );
        size_t s_inside = (s<0) ? -s : ( (s>=(gt_index_type)array_->get_size(3)) ? (2*array_->get_size(3)-s-2) : s );
        size_t p_inside = (p<0) ? -p : ( (p>=(gt_index_type)array_->get_size(4)) ? (2*array_->get_size(4)-p-2) : p );
        size_t r_inside = (r<0) ? -r : ( (r>=(gt_index_type)array_->get_size(5)) ? (2*array_->get_size(5)-r-2) : r );
        size_t a_inside = (a<0) ? -a : ( (a>=(gt_index_type)array_->get_size(6)) ? (2*array_->get_size(6)-a-2) : a );
        size_t q_inside = (q<0) ? -q : ( (q>=(gt_index_type)array_->get_size(7)) ? (2*array_->get_size(7)-q-2) : q );
        size_t u_inside = (u<0) ? -u : ( (u>=(gt_index_type)array_->get_size(8)) ? (2*array_->get_size(8)-u-2) : u );

        return (*array_)(x_inside, y_inside, z_inside, s_inside, p_inside, r_inside, a_inside, q_inside, u_inside);
    }
}
