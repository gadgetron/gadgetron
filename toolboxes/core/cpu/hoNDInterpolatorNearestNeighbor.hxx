/** \file       hoNDInterpolatorNearestNeighbor.h
    \brief      N-dimensional nearest neighbor interpolator

                Designed to work with hoNDArray and hoNDImage

    \author     Hui Xue
*/

namespace Gadgetron
{
    /// hoNDInterpolatorNearestNeighbor

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( const coord_type* pos )
    {
        unsigned int D = array_->get_number_of_dimensions();
        std::vector<long long> ind(D);
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            ind[ii] = static_cast<long long>(pos[ii]+0.5);
        }

        return (*bh_)(ind);
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( const std::vector<coord_type>& pos )
    {
        std::vector<long long> ind(pos.size());
        unsigned int ii;
        unsigned int D = array_->get_number_of_dimensions();
        for ( ii=0; ii<D; ii++ )
        {
            ind[ii] = static_cast<long long>(pos[ii]+0.5);
        }

        return (*bh_)(ind);
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x )
    {
        return (*bh_)(static_cast<long long>(x+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y )
    {
        return (*bh_)(static_cast<long long>(x+0.5), static_cast<long long>(y+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z )
    {
        return (*bh_)(static_cast<long long>(x+0.5), static_cast<long long>(y+0.5), static_cast<long long>(z+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s )
    {
        return (*bh_)(static_cast<long long>(x+0.5), static_cast<long long>(y+0.5), static_cast<long long>(z+0.5), static_cast<long long>(s+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p )
    {
        return (*bh_)(static_cast<long long>(x+0.5), static_cast<long long>(y+0.5), static_cast<long long>(z+0.5), static_cast<long long>(s+0.5), static_cast<long long>(p+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r )
    {
        return (*bh_)(static_cast<long long>(x+0.5), static_cast<long long>(y+0.5), static_cast<long long>(z+0.5), static_cast<long long>(s+0.5), static_cast<long long>(p+0.5), static_cast<long long>(r+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a )
    {
        return (*bh_)(static_cast<long long>(x+0.5), static_cast<long long>(y+0.5), static_cast<long long>(z+0.5), static_cast<long long>(s+0.5), static_cast<long long>(p+0.5), static_cast<long long>(r+0.5), static_cast<long long>(a+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q )
    {
        return (*bh_)(static_cast<long long>(x+0.5), static_cast<long long>(y+0.5), static_cast<long long>(z+0.5), static_cast<long long>(s+0.5), static_cast<long long>(p+0.5), static_cast<long long>(r+0.5), static_cast<long long>(a+0.5), static_cast<long long>(q+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type u )
    {
        return (*bh_)(static_cast<long long>(x+0.5), static_cast<long long>(y+0.5), static_cast<long long>(z+0.5), static_cast<long long>(s+0.5), static_cast<long long>(p+0.5), static_cast<long long>(r+0.5), static_cast<long long>(a+0.5), static_cast<long long>(q+0.5), static_cast<long long>(u+0.5));
    }
}
