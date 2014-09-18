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
        std::vector<gt_index_type> ind(D);
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            ind[ii] = static_cast<gt_index_type>(pos[ii]+0.5);
        }

        return (*bh_)(ind);
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( const std::vector<coord_type>& pos )
    {
        std::vector<gt_index_type> ind(pos.size());
        unsigned int ii;
        unsigned int D = array_->get_number_of_dimensions();
        for ( ii=0; ii<D; ii++ )
        {
            ind[ii] = static_cast<gt_index_type>(pos[ii]+0.5);
        }

        return (*bh_)(ind);
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x )
    {
        return (*bh_)(static_cast<gt_index_type>(x+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y )
    {
        return (*bh_)(static_cast<gt_index_type>(x+0.5), static_cast<gt_index_type>(y+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z )
    {
        return (*bh_)(static_cast<gt_index_type>(x+0.5), static_cast<gt_index_type>(y+0.5), static_cast<gt_index_type>(z+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s )
    {
        return (*bh_)(static_cast<gt_index_type>(x+0.5), static_cast<gt_index_type>(y+0.5), static_cast<gt_index_type>(z+0.5), static_cast<gt_index_type>(s+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p )
    {
        return (*bh_)(static_cast<gt_index_type>(x+0.5), static_cast<gt_index_type>(y+0.5), static_cast<gt_index_type>(z+0.5), static_cast<gt_index_type>(s+0.5), static_cast<gt_index_type>(p+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r )
    {
        return (*bh_)(static_cast<gt_index_type>(x+0.5), static_cast<gt_index_type>(y+0.5), static_cast<gt_index_type>(z+0.5), static_cast<gt_index_type>(s+0.5), static_cast<gt_index_type>(p+0.5), static_cast<gt_index_type>(r+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a )
    {
        return (*bh_)(static_cast<gt_index_type>(x+0.5), static_cast<gt_index_type>(y+0.5), static_cast<gt_index_type>(z+0.5), static_cast<gt_index_type>(s+0.5), static_cast<gt_index_type>(p+0.5), static_cast<gt_index_type>(r+0.5), static_cast<gt_index_type>(a+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q )
    {
        return (*bh_)(static_cast<gt_index_type>(x+0.5), static_cast<gt_index_type>(y+0.5), static_cast<gt_index_type>(z+0.5), static_cast<gt_index_type>(s+0.5), static_cast<gt_index_type>(p+0.5), static_cast<gt_index_type>(r+0.5), static_cast<gt_index_type>(a+0.5), static_cast<gt_index_type>(q+0.5));
    }

    template <typename ArrayType> 
    inline typename hoNDInterpolatorNearestNeighbor<ArrayType>::T hoNDInterpolatorNearestNeighbor<ArrayType>::operator()( coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type u )
    {
        return (*bh_)(static_cast<gt_index_type>(x+0.5), static_cast<gt_index_type>(y+0.5), static_cast<gt_index_type>(z+0.5), static_cast<gt_index_type>(s+0.5), static_cast<gt_index_type>(p+0.5), static_cast<gt_index_type>(r+0.5), static_cast<gt_index_type>(a+0.5), static_cast<gt_index_type>(q+0.5), static_cast<gt_index_type>(u+0.5));
    }
}
