/** \file   hoImageRegDeformationFieldRegister.h
    \brief  Define the class to perform non-rigid image registration to estimate variational deformation field
    \author Hui Xue
*/

#ifndef hoImageRegDeformationFieldRegister_H_
#define hoImageRegDeformationFieldRegister_H_

#pragma once

#include "hoImageRegNonParametricRegister.h"

namespace Gadgetron {

    template<typename TargetType, typename CoordType> 
    class hoImageRegDeformationFieldRegister : public hoImageRegNonParametricRegister<TargetType, TargetType, CoordType>
    {
    public:

        typedef hoImageRegDeformationFieldRegister<TargetType, CoordType> Self;
        typedef hoImageRegNonParametricRegister<TargetType, TargetType, CoordType> BaseClass;

        typedef typename TargetType::value_type ValueType;
        enum { D = TargetType::NDIM };
        enum { DIn = TargetType::NDIM };
        enum { DOut = TargetType::NDIM };

        typedef typename BaseClass::Target2DType Target2DType;
        typedef typename BaseClass::Source2DType Source2DType;

        typedef typename BaseClass::Target3DType Target3DType;
        typedef typename BaseClass::Source3DType Source3DType;

        typedef hoNDImage<CoordType, D> DeformFieldType;

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

        /// boundary handler and interpolator for deformation fields
        typedef hoNDBoundaryHandler<DeformFieldType> BoundaryHandlerDeformFieldType;
        typedef hoNDBoundaryHandlerFixedValue<DeformFieldType> BoundaryHandlerDeformFieldFixedValueType;
        typedef hoNDBoundaryHandlerBorderValue<DeformFieldType> BoundaryHandlerDeformFieldBorderValueType;
        typedef hoNDBoundaryHandlerPeriodic<DeformFieldType> BoundaryHandlerDeformFieldPeriodicType;
        typedef hoNDBoundaryHandlerMirror<DeformFieldType> BoundaryHandlerDeformFieldMirrorType;

        typedef hoNDInterpolator<DeformFieldType> InterpDeformFieldType;
        typedef hoNDInterpolatorLinear<DeformFieldType> InterpDeformFieldLinearType;
        typedef hoNDInterpolatorNearestNeighbor<DeformFieldType> InterpDeformFieldNearestNeighborType;
        typedef hoNDInterpolatorBSpline<DeformFieldType, D> InterpDeformFieldBSplineType;

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
        typedef hoImageRegDeformationFieldSolver<TargetType, TargetType, CoordType> SolverType;

        hoImageRegDeformationFieldRegister(unsigned int resolution_pyramid_levels=3, bool use_world_coordinates=false, ValueType bg_value=ValueType(0));
        virtual ~hoImageRegDeformationFieldRegister();

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
        using BaseClass::solver_type_;

        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

        /// number of iterations for every pyramid level
        std::vector<unsigned int> max_iter_num_pyramid_level_;
        /// threshold for dissimilarity for every pyramid level
        std::vector<ValueType> dissimilarity_thres_pyramid_level_;
        /// number of search size division for every pyramid level
        std::vector<unsigned int> div_num_pyramid_level_;
        /// solver step size for every pyramid level
        std::vector<ValueType> step_size_para_pyramid_level_;
        /// step size division ratio for every pyramid level
        std::vector<ValueType> step_size_div_para_pyramid_level_;
        /// regularization strength for every pyramid level
        /// if regularization_hilbert_strength_world_coordinate_=true, this strength is in the unit of world coordinate
        /// if regularization_hilbert_strength_world_coordinate_=false, this strength is in the unit of pixel
        bool regularization_hilbert_strength_world_coordinate_;
        std::vector< std::vector<ValueType> > regularization_hilbert_strength_pyramid_level_;

        /// in-FOV constraint
        bool apply_in_FOV_constraint_;

        /// divergence free constraint
        bool apply_divergence_free_constraint_;

        /// verbose mode
        bool verbose_;

        /// set the default parameters
        virtual bool setDefaultParameters(unsigned int resolution_pyramid_levels, bool use_world_coordinates);

        /// deformation field transformation, defined in the world coordinate of target image
        TransformationType* transform_;

        /// solver
        std::vector<SolverType> solver_pyramid_;

        /// boundary handler for deformation fields
        BoundaryHandlerDeformFieldType* deform_field_bh_;
        InterpDeformFieldType* deform_field_interp_;

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
        using BaseClass::dissimilarity_pyramid_;
        using BaseClass::warper_pyramid_inverse_;
        using BaseClass::dissimilarity_pyramid_inverse_;

        /// whether the transformation is preset or not
        /// the preset transformation can be used to pass in an initial deformation field
        bool preset_transform_;
    };

    template<typename TargetType, typename CoordType> 
    hoImageRegDeformationFieldRegister<TargetType, CoordType>::
    hoImageRegDeformationFieldRegister(unsigned int resolution_pyramid_levels, bool use_world_coordinates, ValueType bg_value) 
    : transform_(NULL), regularization_hilbert_strength_world_coordinate_(false), verbose_(false), deform_field_bh_(NULL), deform_field_interp_(NULL), preset_transform_(false), BaseClass(resolution_pyramid_levels, bg_value)
    {
        this->setDefaultParameters(resolution_pyramid_levels, use_world_coordinates);
    }

    template<typename TargetType, typename CoordType> 
    hoImageRegDeformationFieldRegister<TargetType, CoordType>::~hoImageRegDeformationFieldRegister()
    {
        if ( !preset_transform_ )
        {
            delete transform_;
            transform_ = NULL;
        }

        if(deform_field_bh_!=NULL)
        {
            delete deform_field_bh_;
            deform_field_bh_ = NULL;
        }

        if (deform_field_interp_ != NULL)
        {
            delete deform_field_interp_;
            deform_field_interp_ = NULL;
        }
    }

    template<typename TargetType, typename CoordType> 
    bool hoImageRegDeformationFieldRegister<TargetType, CoordType>::setDefaultParameters(unsigned int resolution_pyramid_levels, bool use_world_coordinates)
    {
        use_world_coordinates_ = use_world_coordinates;
        resolution_pyramid_levels_ = resolution_pyramid_levels;

        resolution_pyramid_downsample_ratio_.clear();
        resolution_pyramid_downsample_ratio_.resize(resolution_pyramid_levels_-1, std::vector<float>(D, 2.0) );

        resolution_pyramid_blurring_sigma_.clear();
        resolution_pyramid_blurring_sigma_.resize(resolution_pyramid_levels_, std::vector<float>(D, 0.0) );

        boundary_handler_type_warper_.clear();
        // boundary_handler_type_warper_.resize(resolution_pyramid_levels_, GT_BOUNDARY_CONDITION_FIXEDVALUE);
        boundary_handler_type_warper_.resize(resolution_pyramid_levels_, GT_BOUNDARY_CONDITION_BORDERVALUE);

        interp_type_warper_.clear();
        interp_type_warper_.resize(resolution_pyramid_levels_, GT_IMAGE_INTERPOLATOR_LINEAR);

        boundary_handler_type_pyramid_construction_ = GT_BOUNDARY_CONDITION_BORDERVALUE;
        interp_type_pyramid_construction_ = GT_IMAGE_INTERPOLATOR_LINEAR;

        dissimilarity_type_.clear();
        dissimilarity_type_.resize(resolution_pyramid_levels_, GT_IMAGE_DISSIMILARITY_LocalCCR);

        solver_type_.clear();
        solver_type_.resize(resolution_pyramid_levels_, GT_IMAGE_REG_SOLVER_PDE_TIME_INTEGRATION);

        max_iter_num_pyramid_level_.clear();
        max_iter_num_pyramid_level_.resize(resolution_pyramid_levels_, 32);

        dissimilarity_thres_pyramid_level_.clear();
        dissimilarity_thres_pyramid_level_.resize(resolution_pyramid_levels_, 1e-6);

        div_num_pyramid_level_.clear();
        div_num_pyramid_level_.resize(resolution_pyramid_levels_, 2);

        step_size_para_pyramid_level_.clear();
        step_size_para_pyramid_level_.resize(resolution_pyramid_levels_, 0.8);

        step_size_div_para_pyramid_level_.clear();
        step_size_div_para_pyramid_level_.resize(resolution_pyramid_levels_, 0.5);

        regularization_hilbert_strength_world_coordinate_ = false;

        regularization_hilbert_strength_pyramid_level_.clear();
        regularization_hilbert_strength_pyramid_level_.resize(resolution_pyramid_levels_);

        unsigned int ii;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            regularization_hilbert_strength_pyramid_level_[ii].resize(D, 12.0);
        }

        apply_in_FOV_constraint_ = false;
        apply_divergence_free_constraint_ = false;

        verbose_ = false;

        return true;
    }

    template<typename TargetType, typename CoordType> 
    bool hoImageRegDeformationFieldRegister<TargetType, CoordType>::initialize()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(BaseClass::initialize());

            if ( transform_ == NULL )
            {
                std::vector<size_t> dim;
                target_->get_dimensions(dim);
                transform_ = new TransformationType(dim);
                preset_transform_ = false;
            }

            warper_pyramid_.resize(resolution_pyramid_levels_);
            solver_pyramid_.resize(resolution_pyramid_levels_);

            unsigned int ii, jj;
            for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
            {
                warper_pyramid_[ii].setTransformation(*transform_);
                warper_pyramid_[ii].setInterpolator( *source_interp_warper_[ii] );
                warper_pyramid_[ii].setBackgroundValue(bg_value_);
                warper_pyramid_[ii].debugFolder_ = this->debugFolder_;

                solver_pyramid_[ii].setTransform(*transform_);

                if ( regularization_hilbert_strength_world_coordinate_ )
                {
                    // world to pixel
                    std::vector<coord_type> pixelSize;
                    target_->get_pixel_size(pixelSize);

                    for ( jj=0; jj<D; jj++ )
                    {
                        solver_pyramid_[ii].regularization_hilbert_strength_[jj] = (regularization_hilbert_strength_pyramid_level_[ii][jj] / pixelSize[jj]);
                    }
                }
                else
                {
                    for ( jj=0; jj<D; jj++ )
                    {
                        solver_pyramid_[ii].regularization_hilbert_strength_[jj] = regularization_hilbert_strength_pyramid_level_[ii][jj];
                    }
                }

                solver_pyramid_[ii].max_iter_num_ = max_iter_num_pyramid_level_[ii];
                solver_pyramid_[ii].dissimilarity_thres_ = dissimilarity_thres_pyramid_level_[ii];
                solver_pyramid_[ii].div_num_ = div_num_pyramid_level_[ii];
                solver_pyramid_[ii].step_size_para_ = step_size_para_pyramid_level_[ii];
                solver_pyramid_[ii].step_size_div_para_ = step_size_div_para_pyramid_level_[ii];
                solver_pyramid_[ii].verbose_ = verbose_;
                solver_pyramid_[ii].debugFolder_ = this->debugFolder_;

                solver_pyramid_[ii].setTarget(target_pyramid_[ii]);
                solver_pyramid_[ii].setSource(source_pyramid_[ii]);
                solver_pyramid_[ii].setDissimilarity(*dissimilarity_pyramid_[ii]);
                solver_pyramid_[ii].setWarper(warper_pyramid_[ii]);
                solver_pyramid_[ii].setInterpolator(*source_interp_warper_[ii]);
                solver_pyramid_[ii].setBackgroundValue(bg_value_);
                solver_pyramid_[ii].setUseWorldCoordinate(use_world_coordinates_);

                solver_pyramid_[ii].apply_in_FOV_constraint_ = apply_in_FOV_constraint_;
                solver_pyramid_[ii].apply_divergence_free_constraint_ = apply_divergence_free_constraint_;
            }

            // downsample the deformation field if necessary
            if ( !transform_->getDeformationField(0).dimensions_equal(target_pyramid_[resolution_pyramid_levels_-1]) )
            {
                std::vector<size_t> dim;
                target_pyramid_[resolution_pyramid_levels_-1].get_dimensions(dim);

                for ( jj=0; jj<D; jj++ )
                {
                    DeformationFieldType& deField = transform_->getDeformationField(jj);

                    if ( preset_transform_ )
                    {
                        DeformationFieldType deFieldResampled;

                        hoNDBoundaryHandlerBorderValue<DeformationFieldType> bhBorderValue(deField);
                        hoNDInterpolatorLinear<DeformationFieldType> interpLinear(deField, bhBorderValue);

                        GADGET_CHECK_RETURN_FALSE(Gadgetron::resampleImage(deField, interpLinear, dim, deFieldResampled));

                        deField = deFieldResampled;
                        deField.copyImageInfoWithoutImageSize(target_pyramid_[resolution_pyramid_levels_-1]);
                    }
                    else
                    {
                        deField.createFrom(target_pyramid_[resolution_pyramid_levels_-1]);
                        Gadgetron::clear(deField);
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
            GERROR_STREAM("Errors happened in hoImageRegDeformationFieldRegister<TargetType, CoordType>::initialize() ... ");
        }

        return true;
    }

    template<typename TargetType, typename CoordType> 
    bool hoImageRegDeformationFieldRegister<TargetType, CoordType>::performRegistration()
    {
        try
        {
            // starting from the most coarse level

            int level;
            for ( level=(int)resolution_pyramid_levels_-1; level>=0; level-- )
            {
                // update the transform for multi-resolution pyramid
                transform_->update();

                // GADGET_CHECK_RETURN_FALSE(solver_pyramid_[level].initialize());
                GADGET_CHECK_RETURN_FALSE(solver_pyramid_[level].solve());

                if ( !debugFolder_.empty() )
                {
                    unsigned int jj;
                    for ( jj=0; jj<D; jj++ )
                    {
                        std::ostringstream ostr;
                        ostr << "deform_" << jj;

                        gt_exporter_.export_image(transform_->getDeformationField(jj), debugFolder_+ostr.str());
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
                    // Gadgetron::clear(deformExpanded);
                    memset(deformExpanded.begin(), 0, deformExpanded.get_number_of_bytes());

                    if ( downsampledBy2 || resolution_pyramid_divided_by_2_ )
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
                        }
                    }

                    if ( !debugFolder_.empty() )
                    {
                        for ( jj=0; jj<D; jj++ )
                        {
                            std::ostringstream ostr;
                            ostr << "deformExpanded_" << jj;

                            gt_exporter_.export_image(transform_->getDeformationField(jj), debugFolder_+ostr.str());
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationFieldRegister<TargetType, CoordType>::performRegistration() ... ");
        }

        return true;
    }

    template<typename TargetType, typename CoordType> 
    void hoImageRegDeformationFieldRegister<TargetType, CoordType>::printContent(std::ostream& os) const
    {
        using namespace std;
        BaseClass::printContent(os);

        unsigned int ii, jj;

        os << "------------" << std::endl;
        os << "Maximal iteration number for every pyramid level is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << max_iter_num_pyramid_level_[ii] << std::endl;
        }

        os << "------------" << std::endl;
        os << "Threshold for dissimilarity for every pyramid level is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << dissimilarity_thres_pyramid_level_[ii] << std::endl;
        }

        os << "------------" << std::endl;
        os << "Number of search size division for every pyramid level is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << div_num_pyramid_level_[ii] << std::endl;
        }

        os << "------------" << std::endl;
        os << "Solver step size for every pyramid level is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << step_size_para_pyramid_level_[ii] << std::endl;
        }

        os << "------------" << std::endl;
        os << "Step size division ratio for every pyramid level is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << step_size_div_para_pyramid_level_[ii] << std::endl;
        }

        os << "------------" << std::endl;
        if ( regularization_hilbert_strength_world_coordinate_ )
        {
            os << "Regularization strength  is in the unit of physical metric, e.g. mm ... ";
        }
        else
        {
            os << "Regularization strength  is in the unit of image pixel size ... ";
        }

        os << "Regularization strength for every pyramid level is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - [ ";
            for( jj=0; jj<D; jj++ )
            {
                os << regularization_hilbert_strength_pyramid_level_[ii][jj] << " ";
            } 
            os << " ] " << std::endl;
        }
        os << "------------" << std::endl;
        os << "Apply in FOV constraint is : " << apply_in_FOV_constraint_ << std::endl;
        os << "Apply divergence free constraint is : " << apply_divergence_free_constraint_ << std::endl;

        os << "------------" << std::endl;
        os << "Verbose mode is : " << verbose_ << std::endl;
    }

    template<typename TargetType, typename CoordType> 
    void hoImageRegDeformationFieldRegister<TargetType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron non-parametric deformation field image register -------------" << endl;
        this->printContent(os);
        os << "--------------------------------------------------------------------" << endl << ends;
    }
}
#endif // hoImageRegDeformationFieldRegister_H_
