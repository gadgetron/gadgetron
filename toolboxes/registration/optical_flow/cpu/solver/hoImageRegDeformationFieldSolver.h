/** \file   hoImageRegDeformationFieldSolver.h
    \brief  Implement the PDE solver for deformation field non-linear image registration

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

    \author Hui Xue
*/

#ifndef hoImageRegDeformationFieldSolver_H_
#define hoImageRegDeformationFieldSolver_H_

#include "hoImageRegNonParametricSolver.h"
#include "hoImageRegDeformationField.h"

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
    class hoImageRegDeformationFieldSolver : public hoImageRegNonParametricSolver<TargetType, SourceType, CoordType>
    {
    public:

        typedef hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType> Self;
        typedef hoImageRegNonParametricSolver<TargetType, SourceType, CoordType> BaseClass;

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

        typedef hoNDImage< std::complex<CoordType>, D> DeformCplxType;

        typedef typename BaseClass::InterpolatorType InterpolatorType;

        typedef hoImageRegDeformationField<CoordType, D> TransformationType;
        typedef typename TransformationType::input_point_type input_point_type;
        typedef typename TransformationType::output_point_type output_point_type;
        typedef typename TransformationType::jacobian_position_type jacobian_position_type;
        typedef typename TransformationType::DeformationFieldType DeformationFieldType;

        typedef typename BaseClass::ImageRegWarperType ImageRegWarperType;

        typedef typename BaseClass::ImageRegDissimilarityType ImageRegDissimilarityType;

        hoImageRegDeformationFieldSolver();
        virtual ~hoImageRegDeformationFieldSolver();

        void setTransform(TransformationType& transform) { transform_ = &transform; }

        virtual bool initialize();

        virtual bool solve();

        /// perform one iteration of optimization
        virtual bool solve_once(TargetType* target, SourceType* source, TargetType& warped, 
                                unsigned int iter_num, unsigned int max_iter_num, 
                                unsigned int& divTimes, 
                                ValueType& curr_dissimilarity, ValueType& prev_dissimilarity, 
                                TransformationType* transform, ImageRegWarperType& warper, ImageRegDissimilarityType& dissimilarity,
                                bool& stopIteration, 
                                TargetType* gradient_warpped, DeformationFieldType* deform_delta, DeformationFieldType* deform_updated, 
                                DeformationFieldType& deform_norm , DeformationFieldType& deform_norm_one_dim,
                                CoordType* deform_delta_scale_factor);

        /// print function
        virtual void print(std::ostream& os) const;

        /// the regularization method in ref [3] is used
        /// in the unit of pixel
        ValueType regularization_hilbert_strength_[D];

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

        TransformationType* transform_;

        ValueType curr_dissimilarity_;
        ValueType prev_dissimilarity_;

        DeformationFieldType deform_delta_[D];
        DeformationFieldType deform_updated_[D];

        DeformationFieldType deform_norm_;
        DeformationFieldType deform_norm_one_dim_;

        TargetType gradient_warpped_[D];

        DeformCplxType deform_cplx_[D];
        DeformCplxType deform_fft_cplx_[D];
        DeformCplxType deform_fft_buf_cplx_[D];

        /// compensate for the non-isotropic pixel sizes
        coord_type deform_delta_scale_factor_[D];

        using BaseClass::target_;
        using BaseClass::source_;
        using BaseClass::warpped_;
        using BaseClass::bg_value_;
        using BaseClass::interp_;
        using BaseClass::warper_;
        using BaseClass::dissimilarity_;
        using BaseClass::use_world_coordinate_;
    };

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType>::
    hoImageRegDeformationFieldSolver() : BaseClass()
    {
        for ( unsigned int ii=0; ii<D; ii++ )
        {
            regularization_hilbert_strength_[ii] = 12;
            deform_delta_scale_factor_[ii] = 1;
        }
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType>::~hoImageRegDeformationFieldSolver()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType>::initialize()
    {
        GADGET_CHECK_RETURN_FALSE(BaseClass::initialize());
        warper_->setTransformation(*transform_);

        std::vector<size_t> dim;
        target_->get_dimensions(dim);

        deform_norm_.copyImageInfo(*target_);
        deform_norm_one_dim_.copyImageInfo(*target_);

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            deform_delta_[ii].copyImageInfo(*target_);
            Gadgetron::clear(deform_delta_[ii]);

            deform_updated_[ii].copyImageInfo(*target_);
            Gadgetron::clear(deform_updated_[ii]);

            gradient_warpped_[ii].copyImageInfo(*target_);
        }

        deform_delta_scale_factor_[0] = 1;
        for ( ii=0; ii<D; ii++ )
        {
            deform_delta_scale_factor_[ii] = target_->get_pixel_size(0)/target_->get_pixel_size(ii);
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType>::
    solve_once(TargetType* target, SourceType* source, TargetType& warped, 
                unsigned int iter_num, unsigned int max_iter_num, 
                unsigned int& divTimes, 
                ValueType& curr_dissimilarity, ValueType& prev_dissimilarity, 
                TransformationType* transform, ImageRegWarperType& warper, ImageRegDissimilarityType& dissimilarity,
                bool& stopIteration, 
                TargetType* gradient_warpped, DeformationFieldType* deform_delta, DeformationFieldType* deform_updated, 
                DeformationFieldType& deform_norm , DeformationFieldType& deform_norm_one_dim,
                CoordType* deform_delta_scale_factor)
    {
        try
        {
            unsigned int ii;

            long long sx = (long long)(target_->get_size(0));
            long long sy = (long long)(target_->get_size(1));
            long long sz = (long long)(target_->get_size(2));

            long long x, y, z;

            if ( !debugFolder_.empty() )
            {
                for ( ii=0; ii<D; ii++ )
                {
                    std::ostringstream ostr;
                    ostr << "DeformationFieldSolver_deformfield_" << ii;
                    const DeformationFieldType& def = transform->getDeformationField(ii);
                    gt_exporter_.export_image(def, debugFolder_+ostr.str());
                }
            }

            // warp the source

            if ( use_world_coordinate_ )
            {
                GADGET_CHECK_RETURN_FALSE(warper.warpWithDeformationFieldWorldCoordinate(*target, *source, warped));
            }
            else
            {
                GADGET_CHECK_RETURN_FALSE(warper.warp(*target, *source, use_world_coordinate_, warped));
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.export_image(warped, debugFolder_+"DeformationFieldSolver_warpped"); }

            // evaluate the dissimilarity and get the intensity comparison function
            GADGET_CHECK_RETURN_FALSE(dissimilarity.evaluateDeriv(warped));

            curr_dissimilarity = dissimilarity.getDissimilarity();
            if ( verbose_ ) { GDEBUG_STREAM("--> Iteration " << iter_num << " [out of " << max_iter_num << "] : \t" << curr_dissimilarity); }

            if ( prev_dissimilarity < curr_dissimilarity + dissimilarity_thres_ )
            {
                if ( ++divTimes > div_num_ )
                {
                    stopIteration = true;
                    return true;
                }

                step_size_para_ *= step_size_div_para_;

                if ( verbose_ ) { GDEBUG_STREAM("----> Parameter division " << divTimes << " [out of " << div_num_ << "] "); }
            }

            prev_dissimilarity = curr_dissimilarity;

            /// gradient is in the 1/pixel unit
            Gadgetron::gradient(warped, gradient_warpped);

            const TargetType& deriv = dissimilarity.getDeriv();

            size_t N = deriv.get_number_of_elements();
            const ValueType* pD = deriv.begin();

            for ( ii=0; ii<D; ii++ )
            {
                ValueType* pG = gradient_warpped[ii].begin();
                CoordType* pR = deform_delta[ii].begin();

                for (size_t n = 0; n < N; n++)
                {
                    pR[n] = pG[n] * pD[n];
                }
                // Gadgetron::multiply(gradient_warpped[ii], deriv, deform_delta[ii]);
            }

            if ( !debugFolder_.empty() )
            {
                gt_exporter_.export_image(deriv, debugFolder_+"DeformationFieldSolver_deriv");

                for ( ii=0; ii<D; ii++ )
                {
                    std::ostringstream ostr;
                    ostr << "DeformationFieldSolver_gradient_warpped_" << ii;

                    gt_exporter_.export_image(gradient_warpped[ii], debugFolder_+ostr.str());

                    std::ostringstream ostr2;
                    ostr2 << "DeformationFieldSolver_deform_delta_" << ii;

                    gt_exporter_.export_image(deform_delta[ii], debugFolder_+ostr2.str());
                }
            }

            /// compensate for non-isotropic pixel sizes
            for ( ii=0; ii<D; ii++ )
            {
                if ( std::abs(deform_delta_scale_factor[ii]-1) > FLT_EPSILON )
                {
                    Gadgetron::scal(deform_delta_scale_factor[ii], deform_delta[ii]);
                }
            }

            /// filter sigma is in the unit of pixel size
            for ( ii=0; ii<D; ii++ )
            {
                Gadgetron::filterGaussian(deform_delta[ii], regularization_hilbert_strength_);
            }

            if ( !debugFolder_.empty() )
            {
                for ( ii=0; ii<D; ii++ )
                {
                    std::ostringstream ostr;
                    ostr << "DeformationFieldSolver_deform_delta_filtered_" << ii;

                    gt_exporter_.export_image(deform_delta[ii], debugFolder_+ostr.str());
                }
            }

            // compute the max norm of hilbert derivative
            Gadgetron::clear(deform_norm);
            for ( ii=0; ii<D; ii++ )
            {
                Gadgetron::multiply(deform_delta[ii], deform_delta[ii], deform_norm_one_dim);
                Gadgetron::add(deform_norm_one_dim, deform_norm, deform_norm);
            }

            CoordType* pDeformNorm = deform_norm.begin();

            ValueType max_norm_deform_delta = pDeformNorm[0];
            // size_t max_ind;

            for ( ii=1; ii<sx*sy; ii++ )
            {
                if ( max_norm_deform_delta < pDeformNorm[ii] ) max_norm_deform_delta = pDeformNorm[ii];
            }

            // Gadgetron::maxAbsolute(deform_norm, max_norm_deform_delta, max_ind);

            ValueType PDE_time_integration_step_size = 0;
            if ( max_norm_deform_delta > 1e-5 )
            {
                PDE_time_integration_step_size = step_size_para_ / std::sqrt(max_norm_deform_delta);
            }

            if ( PDE_time_integration_step_size > 0 )
            {
                for ( ii=0; ii<D; ii++ )
                {
                    Gadgetron::scal( (CoordType)(PDE_time_integration_step_size), deform_delta[ii]);
                }

                if ( use_world_coordinate_ )
                {
                    // Note: the deform_delta is in the unit of pixel so far, need to convert it to the world coordinate

                    if ( D == 2 )
                    {
                        CoordType ix, iy, wx, wy, pX, pY, deltaWX, deltaWY;

                        // #pragma omp parallel for default(none) private(y, x, ix, iy, wx, wy, pX, pY, deltaWX, deltaWY) shared(sx, sy, target, deform_delta, deform_updated, transform) num_threads(2)
                        for ( y=0; y<sy; y++ )
                        {
                            for ( x=0; x<sx; x++ )
                            {
                                size_t offset = x + y*sx;

                                target->image_to_world( (size_t)x, (size_t)y, wx, wy);

                                CoordType deltaX = deform_delta[0](offset);
                                CoordType deltaY = deform_delta[1](offset);

                                // because the delta deformation is in the pixel size unit, it needs to be converted to world coordinate
                                target->image_to_world( deltaX, deltaY, deltaWX, deltaWY);

                                target->world_to_image(wx+deltaWX, wy+deltaWY, ix, iy);

                                transform->get(ix, iy, pX, pY);

                                deform_updated[0](offset) = deltaWX + pX;
                                deform_updated[1](offset) = deltaWY + pY;
                            }
                        }
                    }
                    else if ( D == 3 )
                    {
                        CoordType ix, iy, iz, wx, wy, wz, pX, pY, pZ, deltaWX, deltaWY, deltaWZ;

                        #pragma omp parallel for default(none) private(y, x, z, ix, iy, iz, wx, wy, wz, pX, pY, pZ, deltaWX, deltaWY, deltaWZ) shared(sx, sy, sz, target, deform_delta, deform_updated, transform)
                        for ( z=0; z<sz; z++ )
                        {
                            for ( y=0; y<sy; y++ )
                            {
                                for ( x=0; x<sx; x++ )
                                {
                                    size_t offset = x + y*sx + z*sx*sy;

                                    target->image_to_world( (size_t)x, (size_t)y, (size_t)z, wx, wy, wz);

                                    CoordType deltaX = deform_delta[0](offset);
                                    CoordType deltaY = deform_delta[1](offset);
                                    CoordType deltaZ = deform_delta[2](offset);

                                    target->image_to_world( deltaX, deltaY, deltaZ, deltaWX, deltaWY, deltaWZ);

                                    target->world_to_image(wx+deltaWX, wy+deltaWY, wz+deltaWZ, ix, iy, iz);

                                    transform->get(ix, iy, iz, pX, pY, pZ);

                                    deform_updated[0](offset) = deltaWX + pX;
                                    deform_updated[1](offset) = deltaWY + pY;
                                    deform_updated[2](offset) = deltaWZ + pZ;
                                }
                            }
                        }
                    }
                    else
                    {
                        size_t N = target_->get_number_of_elements();

                        long long n;

                        #pragma omp parallel default(none) private(n, ii) shared(N, target, deform_delta, deform_updated, transform)
                        {
                            size_t ind[D];
                            CoordType pos[D];
                            CoordType pDelta[D];
                            CoordType pDeltaWorld[D];
                            CoordType indDeform[D];
                            CoordType pDeform[D];

                            #pragma omp for 
                            for ( n=0; n<(long long)N; n++ )
                            {
                                deform_delta[0].calculate_index(n, ind);

                                target->image_to_world( ind, pos);

                                for ( ii=0; ii<D; ii++ )
                                {
                                    pDelta[ii] = deform_delta[ii](n);
                                }

                                target->image_to_world( pDelta, pDeltaWorld);

                                for ( ii=0; ii<D; ii++ )
                                {
                                    pDeltaWorld[ii] += pos[ii];
                                }

                                target->world_to_image(pDeltaWorld, indDeform);
                                transform->get(indDeform, pDeform);

                                for ( ii=0; ii<D; ii++ )
                                {
                                    deform_updated[ii](n) = pDeltaWorld[ii] + pDeform[ii];
                                }
                            }
                        }
                    }
                }
                else
                {
                    if ( D == 2 )
                    {
                        CoordType pX, pY;

                        // #pragma omp parallel for default(none) private(y, x, pX, pY) shared(sx, sy, deform_delta, deform_updated, transform) num_threads(2)
                        for ( y=0; y<sy; y++ )
                        {
                            for ( x=0; x<sx; x++ )
                            {
                                size_t offset = x + y*sx;

                                CoordType deltaX = deform_delta[0](offset);
                                CoordType deltaY = deform_delta[1](offset);

                                transform->get(x+deltaX, y+deltaY, pX, pY);

                                deform_updated[0](offset) = deltaX + pX;
                                deform_updated[1](offset) = deltaY + pY;
                            }
                        }
                    }
                    else if ( D == 3 )
                    {
                        CoordType pX, pY, pZ;

                        #pragma omp parallel for default(none) private(y, x, z, pX, pY, pZ) shared(sx, sy, sz, deform_delta, deform_updated, transform)
                        for ( z=0; z<sz; z++ )
                        {
                            for ( y=0; y<sy; y++ )
                            {
                                for ( x=0; x<sx; x++ )
                                {
                                    size_t offset = x + y*sx + z*sx*sy;

                                    CoordType deltaX = deform_delta[0](offset);
                                    CoordType deltaY = deform_delta[1](offset);
                                    CoordType deltaZ = deform_delta[2](offset);

                                    transform->get(x+deltaX, y+deltaY, z+deltaZ, pX, pY, pZ);

                                    deform_updated[0](offset) = deltaX + pX;
                                    deform_updated[1](offset) = deltaY + pY;
                                    deform_updated[2](offset) = deltaZ + pZ;
                                }
                            }
                        }
                    }
                    else
                    {
                        size_t N = target_->get_number_of_elements();

                        long long n;

                        #pragma omp parallel default(none) private(n, ii) shared(N, deform_delta, deform_updated, transform)
                        {
                            size_t ind[D];
                            CoordType pDelta[D];
                            CoordType indDeform[D];
                            CoordType pDeform[D];

                            #pragma omp for 
                            for ( n=0; n<(long long)N; n++ )
                            {
                                deform_delta[0].calculate_index(n, ind);

                                for ( ii=0; ii<D; ii++ )
                                {
                                    pDelta[ii] = deform_delta[ii](n);
                                    indDeform[ii] = ind[ii] + pDelta[ii];
                                }

                                transform->get(indDeform, pDeform);

                                for ( ii=0; ii<D; ii++ )
                                {
                                    deform_updated[ii](n) = pDelta[ii] + pDeform[ii];
                                }
                            }
                        }
                    }
                }

                if ( !debugFolder_.empty() )
                {
                    for ( ii=0; ii<D; ii++ )
                    {
                        std::ostringstream ostr;
                        ostr << "DeformationFieldSolver_deform_updated_" << ii;
                        gt_exporter_.export_image(deform_updated[ii], debugFolder_+ostr.str());
                    }
                }

                for ( ii=0; ii<D; ii++ )
                {
                    transform->setDeformationField(deform_updated[ii], ii);
                }
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType>::solve()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->initialize());

            prev_dissimilarity_ = std::numeric_limits<ValueType>::max();

            unsigned int divTimes = 0;

            dissimilarity_->initialize(*target_);

            if ( !debugFolder_.empty() )
            {
                gt_exporter_.export_image(*target_, debugFolder_+"DeformationFieldSolver_target");
                gt_exporter_.export_image(*source_, debugFolder_+"DeformationFieldSolver_source");
            }

            bool stopIteration = false;

            if ( verbose_ ) { GDEBUG_STREAM("--> DeformationFieldSolver ... "); }
            for ( iter_num_=0; iter_num_<max_iter_num_; iter_num_++ )
            {
                GADGET_CHECK_RETURN_FALSE( this->solve_once(target_, source_, warpped_, iter_num_, max_iter_num_, 
                                                            divTimes, curr_dissimilarity_, prev_dissimilarity_, 
                                                            transform_, *warper_, *dissimilarity_,
                                                            stopIteration, 
                                                            gradient_warpped_, deform_delta_, deform_updated_, 
                                                            deform_norm_ , deform_norm_one_dim_, deform_delta_scale_factor_) );

                if ( stopIteration ) break;
            }

            if ( verbose_ ) { GDEBUG_STREAM("----> Total iteration number : " << iter_num_); }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType>::solve() ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegDeformationFieldSolver<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image registration non-parametric solver for pixel-wise deformation field -------------" << endl;
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
    }
}
#endif // hoImageRegDeformationFieldSolver_H_
