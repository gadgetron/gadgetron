/** \file   hoImageRegDeformationFieldBidirectionalRegister.h
    \brief  Define the class to perform non-rigid image registration to estimate bi-directional variational deformation field
    \author Hui Xue
*/

#ifndef hoImageRegDeformationFieldBidirectionalRegister_H_
#define hoImageRegDeformationFieldBidirectionalRegister_H_

#pragma once

#include "hoImageRegDeformationFieldRegister.h"

namespace Gadgetron {

    template<typename TargetType, typename CoordType> 
    class hoImageRegDeformationFieldBidirectionalRegister : public hoImageRegDeformationFieldRegister<TargetType, CoordType>
    {
    public:

        typedef hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType> Self;
        typedef hoImageRegDeformationFieldRegister<TargetType, CoordType> BaseClass;
        typedef hoImageRegNonParametricRegister<TargetType, TargetType, CoordType> NonParametricRegisterClass;

        typedef typename TargetType::value_type ValueType;
        enum { D = TargetType::NDIM };
        enum { DIn = TargetType::NDIM };
        enum { DOut = TargetType::NDIM };

        typedef typename BaseClass::Target2DType Target2DType;
        typedef typename BaseClass::Source2DType Source2DType;

        typedef typename BaseClass::Target3DType Target3DType;
        typedef typename BaseClass::Source3DType Source3DType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        /// boundary handler and interpolator for target image
        typedef typename BaseClass::BoundaryHandlerTargetType BoundaryHandlerTargetType;
        typedef typename BaseClass::BoundaryHandlerTargetFixedValueType BoundaryHandlerTargetFixedValueType;
        typedef typename BaseClass::BoundaryHandlerTargetBorderValueType BoundaryHandlerTargetBorderValueType;
        typedef typename BaseClass::BoundaryHandlerTargetPeriodicType BoundaryHandlerTargetPeriodicType;
        typedef typename BaseClass::BoundaryHandlerTargetMirrorType BoundaryHandlerTargetMirrorType;

        typedef typename BaseClass::InterpTargetType InterpTargetType;
        typedef typename BaseClass::InterpTargetLinearType InterpTargetLinearType;
        typedef typename BaseClass::InterpTargetNearestNeighborType InterpTargetNearestNeighborType;
        typedef typename BaseClass::InterpTargetBSplineType InterpTargetBSplineType;

        /// boundary handler and interpolator for source image
        typedef typename BaseClass::BoundaryHandlerSourceType BoundaryHandlerSourceType;
        typedef typename BaseClass::BoundaryHandlerSourceFixedValueType BoundaryHandlerSourceFixedValueType;
        typedef typename BaseClass::BoundaryHandlerSourceBorderValueType BoundaryHandlerSourceBorderValueType;
        typedef typename BaseClass::BoundaryHandlerSourcePeriodicType BoundaryHandlerSourcePeriodicType;
        typedef typename BaseClass::BoundaryHandlerSourceMirrorType BoundaryHandlerSourceMirrorType;

        typedef typename BaseClass::InterpSourceType InterpSourceType;
        typedef typename BaseClass::InterpSourceLinearType InterpSourceLinearType;
        typedef typename BaseClass::InterpSourceNearestNeighborType InterpSourceNearestNeighborType;
        typedef typename BaseClass::InterpSourceBSplineType InterpSourceBSplineType;

        /// warper type
        typedef typename BaseClass::WarperType WarperType;

        /// image dissimilarity type
        typedef typename BaseClass::DissimilarityType DissimilarityType;

        /// transformation type
        typedef hoImageRegDeformationField<CoordType, D> TransformationType;
        typedef typename TransformationType::input_point_type input_point_type;
        typedef typename TransformationType::output_point_type output_point_type;
        typedef typename TransformationType::jacobian_position_type jacobian_position_type;
        typedef typename TransformationType::DeformationFieldType DeformationFieldType;
        typedef typename TransformationType::coord_type coord_type;

        /// solver type
        typedef hoImageRegDeformationFieldBidirectionalSolver<TargetType, TargetType, CoordType> SolverType;

        hoImageRegDeformationFieldBidirectionalRegister(unsigned int resolution_pyramid_levels=3, bool use_world_coordinates=false, ValueType bg_value=ValueType(0));
        virtual ~hoImageRegDeformationFieldBidirectionalRegister();

        /// initialize the registration
        /// should be called after all images and parameters of registration are set
        virtual bool initialize();

        /// perform the registration
        virtual bool performRegistration();

        virtual void printContent(std::ostream& os) const;
        virtual void print(std::ostream& os) const;

        /// parameters

        using BaseClass::use_world_coordinates_;
        using BaseClass::resolution_pyramid_divided_by_2_;
        using BaseClass::resolution_pyramid_levels_;
        using BaseClass::resolution_pyramid_downsample_ratio_;
        using BaseClass::resolution_pyramid_blurring_sigma_;
        using BaseClass::boundary_handler_type_warper_;
        using BaseClass::interp_type_warper_;
        using BaseClass::boundary_handler_type_pyramid_construction_;
        using BaseClass::interp_type_pyramid_construction_;
        using BaseClass::dissimilarity_type_;
        using BaseClass::apply_in_FOV_constraint_;
        using BaseClass::apply_divergence_free_constraint_;
        using BaseClass::solver_type_;
        using BaseClass::deform_field_bh_;
        using BaseClass::deform_field_interp_;

        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

        using BaseClass::max_iter_num_pyramid_level_;
        using BaseClass::dissimilarity_thres_pyramid_level_;
        using BaseClass::div_num_pyramid_level_;
        using BaseClass::step_size_para_pyramid_level_;
        using BaseClass::step_size_div_para_pyramid_level_;
        using BaseClass::regularization_hilbert_strength_world_coordinate_;
        using BaseClass::regularization_hilbert_strength_pyramid_level_;
        using BaseClass::verbose_;

        /// number of iterations to improve the estimation of the inverse transform
        std::vector<unsigned int> inverse_deform_enforce_iter_pyramid_level_;
        /// weight to update the estimation of the inverse transform, must be within [0 1]
        std::vector<CoordType> inverse_deform_enforce_weight_pyramid_level_;

        /// set the default parameters
        virtual bool setDefaultParameters(unsigned int resolution_pyramid_levels, bool use_world_coordinates);

        /// deformation field transformation, defined in the grid of target image
        using BaseClass::transform_;

        TransformationType* transform_inverse_;

        /// solver
        std::vector<SolverType> solver_pyramid_inverse_;

    protected:

        using BaseClass::target_;
        using BaseClass::source_;
        using BaseClass::bg_value_;
        using BaseClass::target_pyramid_;
        using BaseClass::source_pyramid_;
        using BaseClass::target_bh_warper_;
        using BaseClass::target_interp_warper_;
        using BaseClass::source_bh_warper_;
        using BaseClass::source_interp_warper_;
        using BaseClass::target_bh_pyramid_construction_;
        using BaseClass::target_interp_pyramid_construction_;
        using BaseClass::source_bh_pyramid_construction_;
        using BaseClass::source_interp_pyramid_construction_;
        using BaseClass::warper_pyramid_;
        using BaseClass::warper_pyramid_inverse_;
        using BaseClass::dissimilarity_pyramid_;
        using BaseClass::dissimilarity_pyramid_inverse_;
        using BaseClass::preset_transform_;
    };

    template<typename TargetType, typename CoordType> 
    hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType>::
    hoImageRegDeformationFieldBidirectionalRegister(unsigned int resolution_pyramid_levels, bool use_world_coordinates, ValueType bg_value) 
    : BaseClass(resolution_pyramid_levels, use_world_coordinates, bg_value)
    {
        inverse_deform_enforce_iter_pyramid_level_.clear();
        inverse_deform_enforce_iter_pyramid_level_.resize(resolution_pyramid_levels, 10);

        inverse_deform_enforce_weight_pyramid_level_.clear();
        inverse_deform_enforce_weight_pyramid_level_.resize(resolution_pyramid_levels, 0.5);
    }

    template<typename TargetType, typename CoordType> 
    hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType>::~hoImageRegDeformationFieldBidirectionalRegister()
    {
        if ( !preset_transform_ )
        {
            // delete transform_;
            delete transform_inverse_;
            transform_inverse_ = NULL;
        }
    }

    template<typename TargetType, typename CoordType> 
    bool hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType>::setDefaultParameters(unsigned int resolution_pyramid_levels, bool use_world_coordinates)
    {
        BaseClass::setDefaultParameters(resolution_pyramid_levels, use_world_coordinates);

        inverse_deform_enforce_iter_pyramid_level_.clear();
        inverse_deform_enforce_iter_pyramid_level_.resize(resolution_pyramid_levels, 10);

        inverse_deform_enforce_weight_pyramid_level_.clear();
        inverse_deform_enforce_weight_pyramid_level_.resize(resolution_pyramid_levels, 0.5);

        return true;
    }

    template<typename TargetType, typename CoordType> 
    bool hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType>::initialize()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(NonParametricRegisterClass::initialize());

            if ( transform_ == NULL )
            {
                std::vector<size_t> dim;

                target_->get_dimensions(dim);
                transform_ = new TransformationType(dim);

                source_->get_dimensions(dim);
                transform_inverse_ = new TransformationType(dim);

                preset_transform_ = false;
            }

            warper_pyramid_.resize(resolution_pyramid_levels_);
            warper_pyramid_inverse_.resize(resolution_pyramid_levels_);

            solver_pyramid_inverse_.resize(resolution_pyramid_levels_);

            unsigned int ii, jj;
            for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
            {
                warper_pyramid_[ii].setTransformation(*transform_);
                warper_pyramid_[ii].setInterpolator( *source_interp_warper_[ii] );
                warper_pyramid_[ii].setBackgroundValue(bg_value_);
                warper_pyramid_[ii].debugFolder_ = this->debugFolder_;

                warper_pyramid_inverse_[ii].setTransformation(*transform_inverse_);
                warper_pyramid_inverse_[ii].setInterpolator( *target_interp_warper_[ii] );
                warper_pyramid_inverse_[ii].setBackgroundValue(bg_value_);
                warper_pyramid_inverse_[ii].debugFolder_ = this->debugFolder_;

                solver_pyramid_inverse_[ii].setTransform(*transform_);
                solver_pyramid_inverse_[ii].setTransformInverse(*transform_inverse_);

                if ( regularization_hilbert_strength_world_coordinate_ )
                {
                    // world to pixel
                    std::vector<coord_type> pixelSize;
                    target_->get_pixel_size(pixelSize);

                    for ( jj=0; jj<D; jj++ )
                    {
                        solver_pyramid_inverse_[ii].regularization_hilbert_strength_[jj] = (regularization_hilbert_strength_pyramid_level_[ii][jj] / pixelSize[jj]);
                    }
                }
                else
                {
                    for ( jj=0; jj<D; jj++ )
                    {
                        solver_pyramid_inverse_[ii].regularization_hilbert_strength_[jj] = regularization_hilbert_strength_pyramid_level_[ii][jj];
                    }
                }

                solver_pyramid_inverse_[ii].max_iter_num_ = max_iter_num_pyramid_level_[ii];
                solver_pyramid_inverse_[ii].dissimilarity_thres_ = dissimilarity_thres_pyramid_level_[ii];
                solver_pyramid_inverse_[ii].div_num_ = div_num_pyramid_level_[ii];
                solver_pyramid_inverse_[ii].step_size_para_ = step_size_para_pyramid_level_[ii];
                solver_pyramid_inverse_[ii].step_size_div_para_ = step_size_div_para_pyramid_level_[ii];
                solver_pyramid_inverse_[ii].verbose_ = verbose_;
                solver_pyramid_inverse_[ii].debugFolder_ = this->debugFolder_;

                solver_pyramid_inverse_[ii].setTarget(target_pyramid_[ii]);
                solver_pyramid_inverse_[ii].setSource(source_pyramid_[ii]);

                solver_pyramid_inverse_[ii].setDissimilarity(*dissimilarity_pyramid_[ii]);
                solver_pyramid_inverse_[ii].setWarper(warper_pyramid_[ii]);
                solver_pyramid_inverse_[ii].setInterpolator(*source_interp_warper_[ii]);

                solver_pyramid_inverse_[ii].setDissimilarityInverse(*dissimilarity_pyramid_inverse_[ii]);
                solver_pyramid_inverse_[ii].setWarperInverse(warper_pyramid_inverse_[ii]);
                solver_pyramid_inverse_[ii].setInterpolatorInverse(*target_interp_warper_[ii]);

                solver_pyramid_inverse_[ii].setBackgroundValue(bg_value_);
                solver_pyramid_inverse_[ii].setUseWorldCoordinate(use_world_coordinates_);

                solver_pyramid_inverse_[ii].inverse_deform_enforce_iter_ = inverse_deform_enforce_iter_pyramid_level_[ii];
                solver_pyramid_inverse_[ii].inverse_deform_enforce_weight_ = inverse_deform_enforce_weight_pyramid_level_[ii];

                solver_pyramid_inverse_[ii].apply_in_FOV_constraint_ = apply_in_FOV_constraint_;
                solver_pyramid_inverse_[ii].apply_divergence_free_constraint_ = apply_divergence_free_constraint_;
            }

            // downsample the deformation field if necessary
            if ( !transform_->getDeformationField(0).dimensions_equal(target_pyramid_[resolution_pyramid_levels_-1]) )
            {
                std::vector<size_t> dim;
                target_pyramid_[resolution_pyramid_levels_-1].get_dimensions(dim);

                std::vector<size_t> dimInv;
                source_pyramid_[resolution_pyramid_levels_-1].get_dimensions(dimInv);

                for ( jj=0; jj<D; jj++ )
                {
                    DeformationFieldType& deField = transform_->getDeformationField(jj);
                    DeformationFieldType& deField_inverse = transform_inverse_->getDeformationField(jj);

                    if ( preset_transform_ )
                    {
                        // forward
                        DeformationFieldType deFieldResampled;

                        hoNDBoundaryHandlerBorderValue<DeformationFieldType> bhBorderValue(deField);
                        hoNDInterpolatorLinear<DeformationFieldType> interpLinear(deField, bhBorderValue);

                        GADGET_CHECK_RETURN_FALSE(Gadgetron::resampleImage(deField, interpLinear, dim, deFieldResampled));

                        deField = deFieldResampled;
                        deField.copyImageInfoWithoutImageSize(target_pyramid_[resolution_pyramid_levels_-1]);

                        // inverse
                        DeformationFieldType deFieldResampled_inverse;

                        bhBorderValue.setArray(deField_inverse);
                        interpLinear.setArray(deField_inverse);
                        GADGET_CHECK_RETURN_FALSE(Gadgetron::resampleImage(deField_inverse, interpLinear, dimInv, deFieldResampled_inverse));

                        deField_inverse = deFieldResampled_inverse;
                        deField_inverse.copyImageInfoWithoutImageSize(source_pyramid_[resolution_pyramid_levels_-1]);
                    }
                    else
                    {
                        deField.createFrom(target_pyramid_[resolution_pyramid_levels_-1]);
                        Gadgetron::clear(deField);

                        deField_inverse.createFrom(source_pyramid_[resolution_pyramid_levels_-1]);
                        Gadgetron::clear(deField_inverse);
                    }
                }
            }

            // create bh and interp for deformation fields
            if (deform_field_bh_ == NULL)
            {
                deform_field_bh_ = createBoundaryHandler<DeformationFieldType>(this->boundary_handler_type_pyramid_construction_);
            }

            if (deform_field_interp_ == NULL)
            {
                deform_field_interp_ = createInterpolator<DeformationFieldType, D>(this->interp_type_pyramid_construction_);
            }

            deform_field_interp_->setBoundaryHandler(*deform_field_bh_);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType>::initialize() ... ");
        }

        return true;
    }

    template<typename TargetType, typename CoordType> 
    bool hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType>::performRegistration()
    {
        try
        {
            // starting from the most coarse level

            int level;
            for ( level=(int)resolution_pyramid_levels_-1; level>=0; level-- )
            {
                // update the transform for multi-resolution pyramid
                transform_->update();
                transform_inverse_->update();

                GADGET_CHECK_RETURN_FALSE(solver_pyramid_inverse_[level].solve());

                if ( !debugFolder_.empty() )
                {
                    unsigned int jj;
                    for ( jj=0; jj<D; jj++ )
                    {
                        std::ostringstream ostr;
                        ostr << "deform_" << jj;

                        gt_exporter_.export_image(transform_->getDeformationField(jj), debugFolder_+ostr.str());

                        std::ostringstream ostr2;
                        ostr2 << "deform_inverse_" << jj;

                        gt_exporter_.export_image(transform_inverse_->getDeformationField(jj), debugFolder_+ostr2.str());
                    }
                }

                // expand the deformation field for next resolution level
                if ( level>0 )
                {
                    std::vector<float> ratio = resolution_pyramid_downsample_ratio_[level-1];

                    unsigned int jj;
                    bool downsampledBy2 = true;
                    for ( jj=0; jj<D; jj++ )
                    {
                        if ( std::abs(ratio[jj]-2.0f) > FLT_EPSILON )
                        {
                            downsampledBy2 = false;
                            break;
                        }
                    }

                    DeformationFieldType deformExpanded;
                    deformExpanded.createFrom(target_pyramid_[level-1]);
                    Gadgetron::clear(deformExpanded);

                    DeformationFieldType deformInverseExpanded;
                    deformInverseExpanded.createFrom(source_pyramid_[level-1]);
                    Gadgetron::clear(deformInverseExpanded);

                    if ( downsampledBy2 )
                    {
                        for ( jj=0; jj<D; jj++ )
                        {
                            DeformationFieldType& deform = transform_->getDeformationField(jj);
                            Gadgetron::expandImageBy2(deform, *deform_field_bh_, deformExpanded);

                            if ( !use_world_coordinates_ )
                            {
                                Gadgetron::scal(CoordType(2.0), deformExpanded); // the deformation vector should be doubled in length
                            }

                            deform = deformExpanded;

                            DeformationFieldType& deformInv = transform_inverse_->getDeformationField(jj);
                            Gadgetron::expandImageBy2(deformInv, *deform_field_bh_, deformInverseExpanded);

                            if ( !use_world_coordinates_ )
                            {
                                Gadgetron::scal(CoordType(2.0), deformInverseExpanded); // the deformation vector should be doubled in length
                            }

                            deformInv = deformInverseExpanded;
                        }
                    }
                    else
                    {
                        for ( jj=0; jj<D; jj++ )
                        {
                            DeformationFieldType& deform = transform_->getDeformationField(jj);
                            Gadgetron::upsampleImage(deform, *deform_field_interp_, deformExpanded, &ratio[0]);

                            if ( !use_world_coordinates_ )
                            {
                                Gadgetron::scal(CoordType(ratio[jj]), deformExpanded);
                            }

                            deform = deformExpanded;

                            DeformationFieldType& deformInv = transform_inverse_->getDeformationField(jj);
                            Gadgetron::upsampleImage(deformInv, *deform_field_interp_, deformInverseExpanded, &ratio[0]);

                            if ( !use_world_coordinates_ )
                            {
                                Gadgetron::scal(CoordType(ratio[jj]), deformInverseExpanded);
                            }

                            deformInv = deformInverseExpanded;
                        }
                    }
                }

                if ( !debugFolder_.empty() )
                {
                    unsigned int jj;
                    for ( jj=0; jj<D; jj++ )
                    {
                        std::ostringstream ostr;
                        ostr << "deformExpanded_" << jj;

                        gt_exporter_.export_image(transform_->getDeformationField(jj), debugFolder_+ostr.str());

                        std::ostringstream ostr2;
                        ostr2 << "deformExpanded_inverse_" << jj;

                        gt_exporter_.export_image(transform_inverse_->getDeformationField(jj), debugFolder_+ostr2.str());
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType>::performRegistration() ... ");
        }

        return true;
    }

    template<typename TargetType, typename CoordType> 
    void hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType>::printContent(std::ostream& os) const
    {
        using namespace std;
        BaseClass::printContent(os);

        unsigned int ii;

        os << "------------" << std::endl;
        os << "Number of iterations to improve the estimation of the inverse transform is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << inverse_deform_enforce_iter_pyramid_level_[ii] << std::endl;
        }

        os << "------------" << std::endl;
        os << "Weight to update the estimation of the inverse transform is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << inverse_deform_enforce_weight_pyramid_level_[ii] << std::endl;
        }
    }

    template<typename TargetType, typename CoordType> 
    void hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron non-parametric bi-directional deformation field image register -------------" << endl;
        this->printContent(os);
        os << "--------------------------------------------------------------------" << endl << ends;
    }
}
#endif // hoImageRegDeformationFieldBidirectionalRegister_H_
