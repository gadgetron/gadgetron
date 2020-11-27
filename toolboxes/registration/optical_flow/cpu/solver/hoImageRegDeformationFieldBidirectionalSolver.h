/** \file   hoImageRegDeformationFieldBidirectionalSolver.h
    \brief  Implement the PDE solver for bidirecitonal deformation field non-linear image registration

            The PDE solver is a classical gradient descent method, derived from the calculus of variation:

            [1] Gerardo Hermosillo, Christophe Chefd'Hotel, Olivier Faugeras. Variational Methods for Multimodal Image Matching. 
            International Journal of Computer Vision. December 2002, Volume 50, Issue 3, pp 329-343.
            http://link.springer.com/article/10.1023%2FA%3A1020830525823

            [2] Gerardo Hermosillo. Variational Methods for Multimodal Image Matching. PhD Thesis, UNIVERSIT´E DE NICE - SOPHIA ANTIPOLIS. May 2002.
            http://webdocs.cs.ualberta.ca/~dana/readingMedIm/papers/hermosilloPhD.pdf

            [3] Christophe Chefd'Hotel, Gerardo Hermosillo, Olivier D. Faugeras: Flows of diffeomorphisms for multimodal image registration. ISBI 2002: 753-756.
            http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1029367&tag=1

            [4] Christophe Chefd'Hotel, Geometric Methods in Computer Vision and Image Processing : Contributions and Applications. PhD Thesis, April 2005.

            The code is based on the listed source code at page 185 - 187 in ref [2] and extended according to the ref [3] and [4].

            [5] Christoph Guetter, Hui Xue, Christophe Chefd'Hotel, Jens Guehring: Efficient symmetric and inverse-consistent deformable registration through interleaved optimization. ISBI 2011: 590-593.

    \author Hui Xue
*/

#ifndef hoImageRegDeformationFieldBidirectionalSolver_H_
#define hoImageRegDeformationFieldBidirectionalSolver_H_

#pragma once

#include "hoImageRegDeformationFieldSolver.h"

#ifdef max
#undef max
#endif // max

#ifdef min
#undef min
#endif // min

namespace Gadgetron {

    /// ValueType: image pixel value type
    /// CoordType: transformation data type
    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegDeformationFieldBidirectionalSolver : public hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType>
    {
    public:

        typedef hoImageRegDeformationFieldBidirectionalSolver<TargetType, SourceType, CoordType> Self;
        typedef hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType> BaseClass;

        typedef typename TargetType::value_type ValueType;
        enum { D = TargetType::NDIM };
        enum { DIn = TargetType::NDIM };
        enum { DOut = SourceType::NDIM };

        typedef hoNDImage<ValueType, 2> Target2DType;
        typedef Target2DType Source2DType;

        typedef hoNDImage<ValueType, 3> Target3DType;
        typedef Target2DType Source3DType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef CoordType coord_type;

        typedef typename BaseClass::InterpolatorType InterpolatorType;

        typedef hoImageRegDeformationField<CoordType, D> TransformationType;
        typedef typename TransformationType::input_point_type input_point_type;
        typedef typename TransformationType::output_point_type output_point_type;
        typedef typename TransformationType::jacobian_position_type jacobian_position_type;
        typedef typename TransformationType::DeformationFieldType DeformationFieldType;

        typedef typename BaseClass::ImageRegWarperType ImageRegWarperType;

        typedef typename BaseClass::ImageRegDissimilarityType ImageRegDissimilarityType;

        hoImageRegDeformationFieldBidirectionalSolver();
        virtual ~hoImageRegDeformationFieldBidirectionalSolver();

        void setTransform(TransformationType& transform) { transform_ = &transform; }
        void setTransformInverse(TransformationType& transform) { transform_inverse_ = &transform; }

        virtual bool initialize();

        virtual bool solve();

        virtual void print(std::ostream& os) const;

        void setDissimilarityInverse(ImageRegDissimilarityType& dissimilarity) { dissimilarity_inverse_ = &dissimilarity; }
        void setWarperInverse(ImageRegWarperType& warper) { warper_inverse_ = &warper; }
        void setInterpolatorInverse(InterpolatorType& interp) { interp_inverse_ = &interp; }

        virtual bool enforceInverseTransform(TransformationType* transform, TransformationType* transform_inverse, DeformationFieldType* deform_delta_inverse, unsigned int iter_num=10);

        /// number of iterations to improve the estimation of the inverse transform
        unsigned int inverse_deform_enforce_iter_;
        /// weight to update the estimation of the inverse transform, must be within [0 1]
        CoordType inverse_deform_enforce_weight_;

        using BaseClass::regularization_hilbert_strength_;
        using BaseClass::apply_in_FOV_constraint_;
        using BaseClass::apply_divergence_free_constraint_;
        using BaseClass::iter_num_;
        using BaseClass::max_iter_num_;
        using BaseClass::dissimilarity_thres_;
        using BaseClass::parameter_thres_;
        using BaseClass::div_num_;
        using BaseClass::step_size_para_;
        using BaseClass::step_size_div_para_;
        using BaseClass::verbose_;
        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

    protected:

        using BaseClass::transform_;

        using BaseClass::curr_dissimilarity_;
        using BaseClass::prev_dissimilarity_;

        using BaseClass::deform_delta_;
        using BaseClass::deform_updated_;
        using BaseClass::deform_norm_;
        using BaseClass::deform_norm_one_dim_;
        using BaseClass::gradient_warpped_;

        using BaseClass::target_;
        using BaseClass::source_;
        using BaseClass::warpped_;
        using BaseClass::bg_value_;
        using BaseClass::interp_;
        using BaseClass::warper_;
        using BaseClass::dissimilarity_;
        using BaseClass::use_world_coordinate_;
        using BaseClass::deform_delta_scale_factor_;

        /// for the inverse transformation

        SourceType warpped_inverse_;

        InterpolatorType* interp_inverse_;

        ImageRegWarperType* warper_inverse_;

        ImageRegDissimilarityType* dissimilarity_inverse_;

        TransformationType* transform_inverse_;

        ValueType curr_dissimilarity_inverse_;
        ValueType prev_dissimilarity_inverse_;

        DeformationFieldType deform_delta_inverse_[D];
        DeformationFieldType deform_updated_inverse_[D];

        DeformationFieldType deform_norm_inverse_;
        DeformationFieldType deform_norm_one_dim_inverse_;

        TargetType gradient_warpped_inverse_[D];

        coord_type deform_delta_scale_factor_inverse_[D];
    };

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegDeformationFieldBidirectionalSolver<TargetType, SourceType, CoordType>::
    hoImageRegDeformationFieldBidirectionalSolver() : BaseClass(), inverse_deform_enforce_iter_(10), inverse_deform_enforce_weight_(0.5)
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegDeformationFieldBidirectionalSolver<TargetType, SourceType, CoordType>::~hoImageRegDeformationFieldBidirectionalSolver()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegDeformationFieldBidirectionalSolver<TargetType, SourceType, CoordType>::initialize()
    {
        GADGET_CHECK_RETURN_FALSE(interp_inverse_!=NULL);
        GADGET_CHECK_RETURN_FALSE(warper_inverse_!=NULL);
        GADGET_CHECK_RETURN_FALSE(dissimilarity_inverse_!=NULL);
        GADGET_CHECK_RETURN_FALSE(transform_inverse_!=NULL);

        GADGET_CHECK_RETURN_FALSE(BaseClass::initialize());

        warper_inverse_->setInterpolator(*interp_inverse_);
        warper_inverse_->setBackgroundValue(bg_value_);

        dissimilarity_inverse_->setBackgroundValue(bg_value_);

        if ( !warpped_inverse_.dimensions_equal(*source_) )
        {
            warpped_inverse_ = *source_;
        }

        dissimilarity_inverse_->initialize(*source_);

        warper_inverse_->setTransformation(*transform_inverse_);

        std::vector<size_t> dim;
        source_->get_dimensions(dim);

        deform_norm_inverse_.copyImageInfo(*source_);
        deform_norm_one_dim_inverse_.copyImageInfo(*source_);

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            deform_delta_inverse_[ii].copyImageInfo(*source_);
            Gadgetron::clear(deform_delta_[ii]);

            deform_updated_inverse_[ii].copyImageInfo(*source_);
            Gadgetron::clear(deform_updated_[ii]);

            gradient_warpped_inverse_[ii].copyImageInfo(*source_);
        }

        deform_delta_scale_factor_inverse_[0] = 1;
        for ( ii=0; ii<D; ii++ )
        {
            deform_delta_scale_factor_inverse_[ii] = source_->get_pixel_size(0)/source_->get_pixel_size(ii);
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegDeformationFieldBidirectionalSolver<TargetType, SourceType, CoordType>::solve()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->initialize());

            prev_dissimilarity_ = std::numeric_limits<ValueType>::max();
            prev_dissimilarity_inverse_ = std::numeric_limits<ValueType>::max();

            unsigned int divTimes = 0;

            dissimilarity_->initialize(*target_);
            dissimilarity_inverse_->initialize(*source_);

            bool computeForwardTransform = false;
            bool stopIteration = false;

            for ( iter_num_=0; iter_num_<max_iter_num_; iter_num_++ )
            {
                if ( computeForwardTransform )
                {
                    GADGET_CHECK_RETURN_FALSE( this->solve_once(target_, source_, warpped_, iter_num_, max_iter_num_, 
                                                                divTimes, curr_dissimilarity_, prev_dissimilarity_, 
                                                                transform_, *warper_, *dissimilarity_,
                                                                stopIteration, 
                                                                gradient_warpped_, deform_delta_, 
                                                                deform_updated_, deform_norm_, deform_norm_one_dim_,
                                                                deform_delta_scale_factor_) );

                    if ( stopIteration ) break;

                    GADGET_CHECK_RETURN_FALSE(this->enforceInverseTransform(transform_, transform_inverse_, deform_delta_inverse_, 6));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE( this->solve_once(source_, target_, warpped_inverse_, iter_num_, max_iter_num_, 
                                                                divTimes, curr_dissimilarity_inverse_, prev_dissimilarity_inverse_, 
                                                                transform_inverse_, *warper_inverse_, *dissimilarity_inverse_,
                                                                stopIteration, 
                                                                gradient_warpped_inverse_, deform_delta_inverse_, 
                                                                deform_updated_inverse_, deform_norm_inverse_, deform_norm_one_dim_inverse_,
                                                                deform_delta_scale_factor_inverse_) );

                    if ( stopIteration ) break;

                    GADGET_CHECK_RETURN_FALSE(this->enforceInverseTransform(transform_inverse_, transform_, deform_delta_, 6));
                }

                computeForwardTransform = !computeForwardTransform;
            }

            GADGET_CHECK_RETURN_FALSE( this->enforceInverseTransform(transform_inverse_, transform_, deform_delta_, inverse_deform_enforce_iter_) );

            if ( verbose_ ) { GDEBUG_STREAM("----> Total iteration number : " << iter_num_); }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationFieldBidirectionalSolver<TargetType, SourceType, CoordType>::solve() ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegDeformationFieldBidirectionalSolver<TargetType, SourceType, CoordType>::
    enforceInverseTransform(TransformationType* transform, TransformationType* transform_inverse, DeformationFieldType* deform_delta, unsigned int iter_num)
    {
        try
        {

            std::vector<size_t> dim, dim_inverse;

            DeformationFieldType& deform = transform->getDeformationField(0);
            deform.get_dimensions(dim);

            DeformationFieldType& deform_inverse = transform_inverse->getDeformationField(0);
            deform_inverse.get_dimensions(dim_inverse);

            unsigned int iter_enforce;
            for ( iter_enforce=0; iter_enforce<iter_num; iter_enforce++ )
            {
                if ( use_world_coordinate_ )
                {
                    if ( D == 2 )
                    {
                        long long sx = (long long)dim_inverse[0];
                        long long sy = (long long)dim_inverse[1];

                        long long y;
                        // #pragma omp parallel default(none) private(y) shared(sx, sy, transform, transform_inverse, deform_delta, deform, deform_inverse) if(sx*sy>64*1024) num_threads(2)
                        {
                            CoordType ix, iy, px, py, px_inverse, py_inverse, dx, dy, dx_inverse, dy_inverse;
                            size_t offset;

                            // #pragma omp for 
                            for ( y=0; y<(long long)sy; y++ )
                            {
                                for ( size_t x=0; x<sx; x++ )
                                {
                                    transform_inverse->get(x, (size_t)y, dx_inverse, dy_inverse);

                                    deform_inverse.image_to_world(x, y, px_inverse, py_inverse);
                                    px = px_inverse + dx_inverse;
                                    py = py_inverse + dy_inverse;

                                    deform.world_to_image(px, py, ix, iy);

                                    transform->get(ix, iy, dx, dy);

                                    offset = x + y*sx;

                                    deform_delta[0](offset) = dx;
                                    deform_delta[1](offset) = dy;
                                }
                            }
                        }
                    }
                    else if ( D == 3 )
                    {
                        long long sx = (long long)dim_inverse[0];
                        long long sy = (long long)dim_inverse[1];
                        long long sz = (long long)dim_inverse[2];

                        long long z;
                        #pragma omp parallel default(none) private(z) shared(sx, sy, sz, transform, transform_inverse, deform_delta, deform, deform_inverse)
                        {
                            CoordType ix, iy, iz, px, py, pz, px_inverse, py_inverse, pz_inverse, dx, dy, dz, dx_inverse, dy_inverse, dz_inverse;

                            #pragma omp for 
                            for ( z=0; z<(long long)sz; z++ )
                            {
                                for ( size_t y=0; y<sy; y++ )
                                {
                                    size_t offset = z*sx*sy + y*sx;

                                    for ( size_t x=0; x<sx; x++ )
                                    {
                                        transform_inverse->get(x, y, (size_t)z, dx_inverse, dy_inverse, dz_inverse);

                                        deform_inverse.image_to_world(x, y, z, px_inverse, py_inverse, pz_inverse);
                                        px = px_inverse + dx_inverse;
                                        py = py_inverse + dy_inverse;
                                        pz = pz_inverse + dz_inverse;

                                        deform.world_to_image(px, py, pz, ix, iy, iz);

                                        transform->get(ix, iy, iz, dx, dy, dz);

                                        deform_delta[0](offset+x) = dx;
                                        deform_delta[1](offset+x) = dy;
                                        deform_delta[2](offset+x) = dz;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        size_t N = deform_inverse.get_number_of_elements();

                        long long n;
                        #pragma omp parallel default(none) private(n) shared(N, transform, transform_inverse, deform_delta, deform, deform_inverse)
                        {
                            size_t ind[D];
                            CoordType wind[D], wind_inverse[D], d_inverse[D], pt[D], d[D];

                            for ( n=0; n<(long long)N; n++ )
                            {
                                deform_inverse.calculate_index( (unsigned long long)(n), ind);
                                deform_inverse.image_to_world(ind, wind_inverse);

                                transform_inverse->get(ind, d_inverse);

                                unsigned int ii;
                                for ( ii=0; ii<D; ii++ ) pt[ii] = wind_inverse[ii] + d_inverse[ii];

                                deform.world_to_image(pt, wind);

                                transform->get(wind, d);
                                for ( ii=0; ii<D; ii++ )
                                {
                                    deform_delta[ii](n) = d[ii];
                                }
                            }
                        }
                    }
                }
                else
                {
                    if ( D == 2 )
                    {
                        long long sx = (long long)dim_inverse[0];
                        long long sy = (long long)dim_inverse[1];

                        long long y;
                        // #pragma omp parallel default(none) private(y) shared(sx, sy, transform, transform_inverse, deform_delta) if(sx*sy>64*1024) num_threads(2)
                        {
                            CoordType px, py, dx, dy, dx_inverse, dy_inverse;
                            size_t offset;

                            // #pragma omp for 
                            for ( y=0; y<(long long)sy; y++ )
                            {
                                for ( size_t x=0; x<sx; x++ )
                                {
                                    transform_inverse->get(x, (size_t)y, dx_inverse, dy_inverse);

                                    px = x + dx_inverse;
                                    py = y + dy_inverse;

                                    transform->get(px, py, dx, dy);

                                    offset = x + y*sx;

                                    deform_delta[0](offset) = dx;
                                    deform_delta[1](offset) = dy;
                                }
                            }
                        }
                    }
                    else if ( D == 3 )
                    {
                        long long sx = (long long)dim_inverse[0];
                        long long sy = (long long)dim_inverse[1];
                        long long sz = (long long)dim_inverse[2];

                        long long z;
                        #pragma omp parallel default(none) private(z) shared(sx, sy, sz, transform, transform_inverse, deform_delta)
                        {
                            CoordType px, py, pz, dx, dy, dz, dx_inverse, dy_inverse, dz_inverse;
                            size_t offset;

                            #pragma omp for 
                            for ( z=0; z<(long long)sz; z++ )
                            {
                                for ( size_t y=0; y<sy; y++ )
                                {
                                    offset = z*sx*sy + y*sx;

                                    for ( size_t x=0; x<sx; x++ )
                                    {
                                        transform_inverse->get(x, y, (size_t)z, dx_inverse, dy_inverse, dz_inverse);

                                        px = x + dx_inverse;
                                        py = y + dy_inverse;
                                        pz = z + dz_inverse;

                                        transform->get(px, py, pz, dx, dy, dz);

                                        deform_delta[0](offset+x) = dx;
                                        deform_delta[1](offset+x) = dy;
                                        deform_delta[2](offset+x) = dz;
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        size_t N = deform_inverse.get_number_of_elements();

                        long long n;
                        #pragma omp parallel default(none) private(n) shared(N, transform, transform_inverse, deform_delta, deform_inverse)
                        {
                            size_t ind[D];
                            CoordType d_inverse[D], pt[D], d[D];

                            for ( n=0; n<(long long)N; n++ )
                            {
                                deform_inverse.calculate_index( (unsigned long long)(n), ind);

                                transform_inverse->get(ind, d_inverse);

                                unsigned int ii;
                                for ( ii=0; ii<D; ii++ ) pt[ii] = ind[ii] + d_inverse[ii];

                                transform->get(pt, d);
                                for ( ii=0; ii<D; ii++ )
                                {
                                    deform_delta[ii](n) = d[ii];
                                }
                            }
                        }
                    }
                }

                unsigned int ii;
                for ( ii=0; ii<D; ii++ )
                {
                    DeformationFieldType& deform_inverse = transform_inverse->getDeformationField(ii);

                    Gadgetron::scal( CoordType(1-inverse_deform_enforce_weight_), deform_inverse);

                    Gadgetron::scal( CoordType(-1*inverse_deform_enforce_weight_), deform_delta[ii]);

                    Gadgetron::add(deform_delta[ii], deform_inverse, deform_inverse);
                }

                if ( apply_in_FOV_constraint_ )
                {
                    if ( !use_world_coordinate_ )
                    {
                        if ( D == 2 )
                        {
                            long long sx = (long long)dim_inverse[0];
                            long long sy = (long long)dim_inverse[1];

                            DeformationFieldType& dxInv = transform_inverse->getDeformationField(0);
                            DeformationFieldType& dyInv = transform_inverse->getDeformationField(1);

                            long long x, y;
                            // #pragma omp parallel for default(none) private(y, x) shared(sx, sy, dxInv, dyInv) if(sx*sy>64*1024) num_threads(2)
                            for ( y=0; y<sy; y++ )
                            {
                                for ( x=0; x<sx; x++ )
                                {
                                    size_t offset = x + y*sx;

                                    CoordType tx = x + dxInv(offset);
                                    CoordType ty = y + dyInv(offset);

                                    if ( tx < 0 )
                                    {
                                        dxInv(offset) = FLT_EPSILON - x;
                                    }
                                    else if (tx > sx-1 )
                                    {
                                        dxInv(offset) = sx-1-FLT_EPSILON - x;
                                    }

                                    if ( ty < 0 )
                                    {
                                        dyInv(offset) = FLT_EPSILON - y;
                                    }
                                    else if (ty > sy-1 )
                                    {
                                        dyInv(offset) = sy-1-FLT_EPSILON - y;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationFieldBidirectionalSolver<TargetType, SourceType, CoordType>::enforceInverseTransform(...) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegDeformationFieldBidirectionalSolver<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image registration non-parametric solver for pixel-wise bidirectional deformation field -------------" << endl;
        os << "Image dimension is : " << D << endl;
        os << "Image data type is : " << std::string(typeid(ValueType).name()) << std::endl;
        os << "Transformation data type is : " << std::string(typeid(CoordType).name()) << std::endl;
        os << "Use world coordinate is : " << use_world_coordinate_ << std::endl;
        os << "Maximal iteration number is : " << max_iter_num_ << std::endl;
        os << "Dissimilarity threshold is : " << dissimilarity_thres_ << std::endl;
        os << "Parameter threshold is : " << parameter_thres_ << std::endl;
        os << "Number of search size division is : " << div_num_ << std::endl;
        os << "Solver step size is : " << step_size_para_ << std::endl;
        os << "Step size division ratio is : " << step_size_div_para_ << std::endl;
        os << "Step size division ratio is : " << step_size_div_para_ << std::endl;
        os << "Number of iterations to improve the estimation of the inverse transform is : " << inverse_deform_enforce_iter_ << std::endl;
        os << "Weight to update the estimation of the inverse transform is : " << inverse_deform_enforce_weight_ << std::endl;
    }
}
#endif // hoImageRegDeformationFieldBidirectionalSolver_H_
