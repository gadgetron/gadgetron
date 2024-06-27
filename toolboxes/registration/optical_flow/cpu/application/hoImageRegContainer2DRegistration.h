/** \file   hoImageRegContainer2DRegistration.h
    \brief  Define the class to perform image registration over a 2D image container
    \author Hui Xue
*/

#ifndef hoImageRegContainer2DRegistration_H_
#define hoImageRegContainer2DRegistration_H_

#pragma once

#include <sstream>
#include "hoNDArray.h"
#include "hoNDImage.h"
#include "hoMRImage.h"
#include "hoNDInterpolator.h"
#include "hoNDBoundaryHandler.h"
#include "hoMatrix.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"

// transformation
#include "hoImageRegTransformation.h"
#include "hoImageRegParametricTransformation.h"
#include "hoImageRegDeformationField.h"

// warper
#include "hoImageRegWarper.h"

// solver
#include "hoImageRegDeformationFieldSolver.h"
#include "hoImageRegParametricSolver.h"
#include "hoImageRegDeformationFieldBidirectionalSolver.h"

// dissimilarity
#include "hoImageRegDissimilaritySSD.h"
#include "hoImageRegDissimilarityLocalCCR.h"
#include "hoImageRegDissimilarityMutualInformation.h"
#include "hoImageRegDissimilarityNormalizedMutualInformation.h"

// register
#include "hoImageRegDeformationFieldRegister.h"
#include "hoImageRegDeformationFieldBidirectionalRegister.h"

// container2D
#include "hoNDImageContainer2D.h"

namespace Gadgetron {

    template <typename ObjType> void printInfo(const ObjType& obj)
    {
        std::ostringstream outs;
        obj.print(outs);
        outs << std::ends;
        std::string msg(outs.str());
        GDEBUG_STREAM(msg.c_str());
    }

    enum GT_IMAGE_REG_CONTAINER_MODE
    {
        GT_IMAGE_REG_CONTAINER_PAIR_WISE,
        GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE,
        GT_IMAGE_REG_CONTAINER_PROGRESSIVE
    };

    inline std::string getImageRegContainerModeName(GT_IMAGE_REG_CONTAINER_MODE v)
    {
        std::string name;

        switch (v)
        {
            case GT_IMAGE_REG_CONTAINER_PAIR_WISE:
                name = "Pair-wise";
                break;

            case GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE:
                name = "FixedReference";
                break;

            case GT_IMAGE_REG_CONTAINER_PROGRESSIVE:
                name = "Progressive";
                break;

            default:
                GERROR_STREAM("Unrecognized image registration container mode type : " << v);
        }

        return name;
    }

    inline GT_IMAGE_REG_CONTAINER_MODE getImageRegContainerModeType(const std::string& name)
    {
        GT_IMAGE_REG_CONTAINER_MODE v;

        if ( name == "Pair-wise" )
        {
            v = GT_IMAGE_REG_CONTAINER_PAIR_WISE;
        }
        else if ( name == "FixedReference" )
        {
            v = GT_IMAGE_REG_CONTAINER_FIXED_REFERENCE;
        }
        else if ( name == "Progressive" )
        {
            v = GT_IMAGE_REG_CONTAINER_PROGRESSIVE;
        }
        else
        {
            GERROR_STREAM("Unrecognized image registration container mode name : " << name);
        }

        return v;
    }

    /// perform the image registration over an image container2D
    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegContainer2DRegistration
    {
    public:

        typedef hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType> Self;
        typedef typename TargetType::value_type ValueType;
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

        /// boundary handler and interpolator for target image
        typedef hoNDBoundaryHandler<TargetType> BoundaryHandlerTargetType;
        typedef hoNDBoundaryHandlerFixedValue<TargetType> BoundaryHandlerTargetFixedValueType;
        typedef hoNDBoundaryHandlerBorderValue<TargetType> BoundaryHandlerTargetBorderValueType;
        typedef hoNDBoundaryHandlerPeriodic<TargetType> BoundaryHandlerTargetPeriodicType;
        typedef hoNDBoundaryHandlerMirror<TargetType> BoundaryHandlerTargetMirrorType;

        typedef hoNDInterpolator<TargetType> InterpTargetType;
        typedef hoNDInterpolatorLinear<TargetType> InterpTargetLinearType;
        typedef hoNDInterpolatorNearestNeighbor<TargetType> InterpTargetNearestNeighborType;
        typedef hoNDInterpolatorBSpline<TargetType, DIn> InterpTargetBSplineType;

        /// boundary handler and interpolator for source image
        typedef hoNDBoundaryHandler<SourceType> BoundaryHandlerSourceType;
        typedef hoNDBoundaryHandlerFixedValue<SourceType> BoundaryHandlerSourceFixedValueType;
        typedef hoNDBoundaryHandlerBorderValue<SourceType> BoundaryHandlerSourceBorderValueType;
        typedef hoNDBoundaryHandlerPeriodic<SourceType> BoundaryHandlerSourcePeriodicType;
        typedef hoNDBoundaryHandlerMirror<SourceType> BoundaryHandlerSourceMirrorType;

        typedef hoNDInterpolator<SourceType> InterpSourceType;
        typedef hoNDInterpolatorLinear<SourceType> InterpSourceLinearType;
        typedef hoNDInterpolatorNearestNeighbor<SourceType> InterpSourceNearestNeighborType;
        typedef hoNDInterpolatorBSpline<SourceType, DIn> InterpSourceBSplineType;

        /// warper type
        typedef hoImageRegWarper<TargetType, SourceType, CoordType> WarperType;

        /// image dissimilarity type
        typedef hoImageRegDissimilarity<SourceType> DissimilarityType;

        /// transformation
        typedef hoImageRegParametricTransformation<CoordType, DIn, DOut> TransformationParametricType;

        typedef hoImageRegDeformationField<CoordType, DIn> TransformationDeformationFieldType;
        typedef typename TransformationDeformationFieldType::input_point_type input_point_type;
        typedef typename TransformationDeformationFieldType::output_point_type output_point_type;
        typedef typename TransformationDeformationFieldType::jacobian_position_type jacobian_position_type;
        typedef typename TransformationDeformationFieldType::DeformationFieldType DeformationFieldType;

        /// container
        typedef hoNDImageContainer2D<TargetType> TargetContinerType;
        typedef hoNDImageContainer2D<SourceType> SourceContinerType;
        typedef hoNDImageContainer2D<DeformationFieldType> DeformationFieldContinerType;

        hoImageRegContainer2DRegistration(unsigned int resolution_pyramid_levels=3, bool use_world_coordinates=false, ValueType bg_value=ValueType(0));
        virtual ~hoImageRegContainer2DRegistration();

        /// set the default parameters
        virtual bool setDefaultParameters(unsigned int resolution_pyramid_levels=3, bool use_world_coordinates=false);

        /// register two images
        /// transform or deform can contain the initial transformation or deformation
        /// if warped == NULL, warped images will not be computed
        virtual bool registerTwoImagesParametric(const TargetType& target, const SourceType& source, bool initial, TargetType* warped, TransformationParametricType& transform);
        virtual bool registerTwoImagesDeformationField(const TargetType& target, const SourceType& source, bool initial, TargetType* warped, DeformationFieldType** deform);
        virtual bool registerTwoImagesDeformationFieldBidirectional(const TargetType& target, const SourceType& source, bool initial, TargetType* warped, DeformationFieldType** deform, DeformationFieldType** deformInv);

        /// if warped is true, the warped images will be computed; if initial is true, the registration will be initialized by deformation_field_ and deformation_field_inverse_
        virtual bool registerOverContainer2DPairWise(TargetContinerType& targetContainer, SourceContinerType& sourceContainer, bool warped, bool initial = false);
        virtual bool registerOverContainer2DFixedReference(TargetContinerType& targetContainer, const std::vector<unsigned int>& referenceFrame, bool warped, bool initial = false);
        virtual bool registerOverContainer2DProgressive(TargetContinerType& targetContainer, const std::vector<unsigned int>& referenceFrame);

        /// warp image containers
        template <typename TargetType2, typename SourceType2> 
        bool warpContainer2D(const hoNDImageContainer2D< TargetType2 >& targetContainer, 
                             const hoNDImageContainer2D< SourceType2 >& sourceContainer, 
                             DeformationFieldContinerType deformation_field[], 
                             hoNDImageContainer2D< SourceType2 >& warppedContainer,
                             Gadgetron::GT_BOUNDARY_CONDITION bh=GT_BOUNDARY_CONDITION_FIXEDVALUE)
        {
            try
            {
                typedef typename TargetType2::value_type ValueType2;

                typedef TargetType2 ImageTargetType;
                typedef SourceType2 ImageSourceType;

                size_t R = sourceContainer.rows();
                std::vector<size_t> cols = sourceContainer.cols();

                GADGET_CHECK_RETURN_FALSE(targetContainer.dimensions_equal_container(sourceContainer));
                GADGET_CHECK_RETURN_FALSE(targetContainer.dimensions_equal_container(deformation_field[0]));

                if ( !targetContainer.dimensions_equal_container(warppedContainer) )
                {
                    GADGET_CHECK_RETURN_FALSE(warppedContainer.copyFrom(targetContainer));
                }

                if ( R == 1 )
                {
                    long long N = (long long)cols[0];

                    long long c;
                    #pragma omp parallel private(c) shared(N, targetContainer, sourceContainer, warppedContainer, deformation_field, bh) if ( DIn==2 )
                    {
                        hoImageRegDeformationField<CoordType, DIn> deformTransform;
                        hoNDBoundaryHandlerFixedValue< ImageSourceType > bhFixedValue;
                        hoNDBoundaryHandlerBorderValue< ImageSourceType > bhBorderValue;
                        hoNDBoundaryHandlerPeriodic< ImageSourceType > bhPeriodic;
                        hoNDBoundaryHandlerMirror< ImageSourceType > bhMirror;

                        hoNDInterpolatorBSpline<ImageSourceType, DIn> interpBSpline(5);

                        hoImageRegWarper<ImageTargetType, ImageSourceType, CoordType> warper;
                        warper.setBackgroundValue(bg_value_);
                        warper.setTransformation(deformTransform);
                        warper.setInterpolator(interpBSpline);

                        #pragma omp for 
                        for ( c=0; c<N; c++ )
                        {
                            const ImageTargetType& target = targetContainer(0, c);
                            ImageSourceType& source = const_cast<ImageSourceType&>(sourceContainer(0, c));
                            ImageTargetType& warpped = warppedContainer(0, c);

                            bhFixedValue.setArray( source );
                            interpBSpline.setArray( source );

                            if ( bh == GT_BOUNDARY_CONDITION_FIXEDVALUE )
                                interpBSpline.setBoundaryHandler(bhFixedValue);
                            else if ( bh == GT_BOUNDARY_CONDITION_BORDERVALUE )
                                interpBSpline.setBoundaryHandler(bhBorderValue);
                            else if ( bh == GT_BOUNDARY_CONDITION_PERIODIC )
                                interpBSpline.setBoundaryHandler(bhPeriodic);
                            else if ( bh == GT_BOUNDARY_CONDITION_MIRROR )
                                interpBSpline.setBoundaryHandler(bhMirror);
                            else
                                interpBSpline.setBoundaryHandler(bhFixedValue);

                            for ( unsigned int ii=0; ii<DIn; ii++ )
                            {
                                deformTransform.setDeformationField( deformation_field[ii](0, c), ii );
                            }

                            warper.warp(target, source, use_world_coordinates_, warpped);
                        }
                    }
                }
                else
                {

                    long long r, c;
                    #pragma omp parallel default(none) private(r, c) shared(targetContainer, sourceContainer, warppedContainer, deformation_field, R, cols, bh) if ( DIn==2 )
                    {
                        hoImageRegDeformationField<CoordType, DIn> deformTransform;
                        hoNDBoundaryHandlerFixedValue< ImageSourceType > bhFixedValue;
                        hoNDBoundaryHandlerBorderValue< ImageSourceType > bhBorderValue;
                        hoNDBoundaryHandlerPeriodic< ImageSourceType > bhPeriodic;
                        hoNDBoundaryHandlerMirror< ImageSourceType > bhMirror;

                        hoNDInterpolatorBSpline<ImageSourceType, DIn> interpBSpline(5);

                        hoImageRegWarper<ImageTargetType, ImageSourceType, CoordType> warper;
                        warper.setBackgroundValue(bg_value_);
                        warper.setTransformation(deformTransform);
                        warper.setInterpolator(interpBSpline);

                        #pragma omp for 
                        for ( r=0; r<(long long)R; r++ )
                        {
                            long long N = (long long)cols[r];
                            for ( c=0; c<N; c++ )
                            {
                                const ImageTargetType& target = targetContainer(r, c);
                                ImageSourceType& source = const_cast<ImageSourceType&>(sourceContainer(r, c));
                                ImageTargetType& warpped = warppedContainer(r, c);

                                bhFixedValue.setArray( source );
                                interpBSpline.setArray( source );

                                if ( bh == GT_BOUNDARY_CONDITION_FIXEDVALUE )
                                    interpBSpline.setBoundaryHandler(bhFixedValue);
                                else if ( bh == GT_BOUNDARY_CONDITION_BORDERVALUE )
                                    interpBSpline.setBoundaryHandler(bhBorderValue);
                                else if ( bh == GT_BOUNDARY_CONDITION_PERIODIC )
                                    interpBSpline.setBoundaryHandler(bhPeriodic);
                                else if ( bh == GT_BOUNDARY_CONDITION_MIRROR )
                                    interpBSpline.setBoundaryHandler(bhMirror);
                                else
                                    interpBSpline.setBoundaryHandler(bhFixedValue);

                                for ( unsigned int ii=0; ii<DIn; ii++ )
                                {
                                    deformTransform.setDeformationField( deformation_field[ii](r, c), ii );
                                }

                                warper.warp(target, source, use_world_coordinates_, warpped);
                            }
                        }
                    }
                }
            }
            catch(...)
            {
                GERROR_STREAM("Errors happened in hoImageRegContainer2DRegistration<...>::warpContainer2D(...) ... ");
                return false;
            }

            return true;
        }

        /// print the class information
        virtual void print(std::ostream& os) const;

        // ----------------------------------
        // parameters
        // ----------------------------------

        /// mode for registration over the container
        GT_IMAGE_REG_CONTAINER_MODE container_reg_mode_;

        /// mode for transformation
        GT_IMAGE_REG_TRANSFORMATION container_reg_transformation_;

        /// back ground values, used to mark regions in the target image which will not be warped
        ValueType bg_value_;

        /// whether to perform world coordinate registration
        bool use_world_coordinates_;

        /// number of resolution pyramid levels
        unsigned int resolution_pyramid_levels_;

        /// number of iterations for every pyramid level
        std::vector<unsigned int> max_iter_num_pyramid_level_;

        /// dissimilarity
        GT_IMAGE_DISSIMILARITY dissimilarity_type_;

        /// threshold for dissimilarity for every pyramid level
        std::vector<ValueType> dissimilarity_thres_pyramid_level_;

        /// number of search size division for every pyramid level
        std::vector<unsigned int> div_num_pyramid_level_;

        /// parameters for dissimilarity measures, for every paramid level
        /// LocalCCR
        std::vector<std::vector<ValueType> > dissimilarity_LocalCCR_sigmaArg_;

        /// Histogram based
        /// Mutual information
        std::vector<ValueType> dissimilarity_MI_betaArg_;

        /// regularization strength for every pyramid level
        /// if regularization_hilbert_strength_world_coordinate_=true, this strength is in the unit of world coordinate
        /// if regularization_hilbert_strength_world_coordinate_=false, this strength is in the unit of pixel
        bool regularization_hilbert_strength_world_coordinate_;
        std::vector< std::vector<ValueType> > regularization_hilbert_strength_pyramid_level_;

        /// boundary handler type
        std::vector<GT_BOUNDARY_CONDITION> boundary_handler_type_warper_;
        std::vector<GT_IMAGE_INTERPOLATOR> interp_type_warper_;

        /// number of iterations to improve the estimation of the inverse transform
        std::vector<unsigned int> inverse_deform_enforce_iter_pyramid_level_;
        /// weight to update the estimation of the inverse transform, must be within [0 1]
        std::vector<CoordType> inverse_deform_enforce_weight_pyramid_level_;

        /// in-FOV constraint
        bool apply_in_FOV_constraint_;

        /// divergence free constraint
        bool apply_divergence_free_constraint_;

        /// verbose mode
        bool verbose_;

        // ----------------------------------
        // debug and timing
        // ----------------------------------
        // clock for timing
        Gadgetron::GadgetronTimer gt_timer1_;
        Gadgetron::GadgetronTimer gt_timer2_;
        Gadgetron::GadgetronTimer gt_timer3_;

        bool performTiming_;

        // exporter
        Gadgetron::ImageIOAnalyze gt_exporter_;

        // debug folder
        std::string debugFolder_;

        // ----------------------------------
        // registration results
        // ----------------------------------

        /// warpped images
        TargetContinerType warped_container_;

        /// for parametric registration
        std::vector< std::vector<TransformationParametricType*> > parametric_tranformation_;

        /// deformation field registration
        DeformationFieldContinerType deformation_field_[DIn];
        DeformationFieldContinerType deformation_field_inverse_[DIn];

    protected:

        bool initialize(const TargetContinerType& targetContainer, bool warped);

    };

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::
    hoImageRegContainer2DRegistration(unsigned int resolution_pyramid_levels, bool use_world_coordinates, ValueType bg_value) 
    : bg_value_(bg_value), use_world_coordinates_(use_world_coordinates), resolution_pyramid_levels_(resolution_pyramid_levels), performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);

        GADGET_CHECK_THROW(this->setDefaultParameters(resolution_pyramid_levels, use_world_coordinates));
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::
    ~hoImageRegContainer2DRegistration()
    {
        if ( !parametric_tranformation_.empty() )
        {
            size_t r, c;
            for ( r=0; r<parametric_tranformation_.size(); r++ )
            {
                if ( !parametric_tranformation_[r].empty() )
                {
                    for ( c=0; c<parametric_tranformation_[r].size(); c++ )
                    {
                        if ( parametric_tranformation_[r][c] != NULL )
                        {
                            delete parametric_tranformation_[r][c];
                            parametric_tranformation_[r][c] = NULL;
                        }
                    }
                }
            }
        }
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::setDefaultParameters(unsigned int resolution_pyramid_levels, bool use_world_coordinates)
    {
        unsigned int ii;

        use_world_coordinates_ = use_world_coordinates;
        resolution_pyramid_levels_ = resolution_pyramid_levels;

        container_reg_mode_ = GT_IMAGE_REG_CONTAINER_PAIR_WISE;
        container_reg_transformation_ = GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD;

        max_iter_num_pyramid_level_.clear();
        max_iter_num_pyramid_level_.resize(resolution_pyramid_levels_, 32);
        max_iter_num_pyramid_level_[0] = 16;

        dissimilarity_type_ = GT_IMAGE_DISSIMILARITY_LocalCCR;

        dissimilarity_thres_pyramid_level_.clear();
        dissimilarity_thres_pyramid_level_.resize(resolution_pyramid_levels_, (ValueType)(1e-5) );

        div_num_pyramid_level_.clear();
        div_num_pyramid_level_.resize(resolution_pyramid_levels_, 2);

        dissimilarity_LocalCCR_sigmaArg_.clear();
        dissimilarity_LocalCCR_sigmaArg_.resize(resolution_pyramid_levels_);
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            dissimilarity_LocalCCR_sigmaArg_[ii].resize(DIn, 2.0);
        }

        dissimilarity_MI_betaArg_.clear();
        dissimilarity_MI_betaArg_.resize(resolution_pyramid_levels_, 2);

        regularization_hilbert_strength_world_coordinate_ = false;
        regularization_hilbert_strength_pyramid_level_.resize(resolution_pyramid_levels_);
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            regularization_hilbert_strength_pyramid_level_[ii].resize(DIn, 12.0);
        }

        boundary_handler_type_warper_.clear();
        boundary_handler_type_warper_.resize(resolution_pyramid_levels_, GT_BOUNDARY_CONDITION_BORDERVALUE);

        interp_type_warper_.clear();
        interp_type_warper_.resize(resolution_pyramid_levels_, GT_IMAGE_INTERPOLATOR_LINEAR);

        inverse_deform_enforce_iter_pyramid_level_.clear();
        inverse_deform_enforce_iter_pyramid_level_.resize(resolution_pyramid_levels_, 10);

        inverse_deform_enforce_weight_pyramid_level_.clear();
        inverse_deform_enforce_weight_pyramid_level_.resize(resolution_pyramid_levels_, 0.5);

        apply_in_FOV_constraint_ = false;
        apply_divergence_free_constraint_ = false;

        verbose_ = false;

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::
    registerTwoImagesParametric(const TargetType& target, const SourceType& source, bool initial, TargetType* warped, TransformationParametricType& transform)
    {
        try
        {
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::registerTwoImagesParametric(...) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::
    registerTwoImagesDeformationField(const TargetType& target, const SourceType& source, bool initial, TargetType* warped, DeformationFieldType** deform)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(deform!=NULL);

            hoImageRegDeformationFieldRegister<TargetType, CoordType> reg(resolution_pyramid_levels_, use_world_coordinates_, bg_value_);

            if ( !debugFolder_.empty() )
            {
                reg.debugFolder_ = debugFolder_;
            }

            GADGET_CHECK_RETURN_FALSE(reg.setDefaultParameters(resolution_pyramid_levels_, use_world_coordinates_));

            reg.max_iter_num_pyramid_level_ = max_iter_num_pyramid_level_;
            reg.div_num_pyramid_level_ = div_num_pyramid_level_;
            reg.dissimilarity_MI_betaArg_ = dissimilarity_MI_betaArg_;
            reg.regularization_hilbert_strength_world_coordinate_ = regularization_hilbert_strength_world_coordinate_;
            reg.regularization_hilbert_strength_pyramid_level_ = regularization_hilbert_strength_pyramid_level_;
            reg.dissimilarity_LocalCCR_sigmaArg_ = dissimilarity_LocalCCR_sigmaArg_;
            reg.boundary_handler_type_warper_ = boundary_handler_type_warper_;
            reg.interp_type_warper_ = interp_type_warper_;
            reg.apply_in_FOV_constraint_ = apply_in_FOV_constraint_;
            reg.apply_divergence_free_constraint_ = apply_divergence_free_constraint_;
            reg.verbose_ = verbose_;

            reg.dissimilarity_type_.clear();
            reg.dissimilarity_type_.resize(resolution_pyramid_levels_, dissimilarity_type_);

            reg.setTarget( const_cast<TargetType&>(target) );
            reg.setSource( const_cast<TargetType&>(source) );

            if ( verbose_ )
            {
                std::ostringstream outs;
                reg.print(outs);
                GDEBUG_STREAM(outs.str());
            }

            GADGET_CHECK_RETURN_FALSE(reg.initialize());

            unsigned int d;

            if ( target.dimensions_equal( *(deform[0]) ) )
            {
                if ( initial )
                {
                    for ( d=0; d<DIn; d++ )
                    {
                        reg.transform_->setDeformationField( *(deform[d]), d);
                    }
                }
            }
            else
            {
                for ( d=0; d<DIn; d++ )
                {
                    deform[d]->copyImageInfo(target);
                    Gadgetron::clear( *(deform[d]) );
                }
            }

            GADGET_CHECK_RETURN_FALSE(reg.performRegistration());

            for ( d=0; d<DIn; d++ )
            {
                *(deform[d]) = reg.transform_->getDeformationField(d);
            }

            if ( warped != NULL )
            {
                /// bspline warp
                hoNDBoundaryHandlerFixedValue<SourceType> bhFixedValue;
                bhFixedValue.setArray( const_cast<SourceType&>(source) );

                hoNDInterpolatorBSpline<SourceType, DIn> interpBSpline(5);
                interpBSpline.setArray( const_cast<SourceType&>(source) );
                interpBSpline.setBoundaryHandler(bhFixedValue);

                hoImageRegWarper<TargetType, SourceType, CoordType> warper;
                warper.setBackgroundValue(bg_value_);
                warper.setTransformation(*reg.transform_);
                warper.setInterpolator(interpBSpline);

                warper.warp(target, source, use_world_coordinates_, *warped);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::registerTwoImagesDeformationField(...) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::
    registerTwoImagesDeformationFieldBidirectional(const TargetType& target, const SourceType& source, bool initial, TargetType* warped, DeformationFieldType** deform, DeformationFieldType** deformInv)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(deform!=NULL);
            GADGET_CHECK_RETURN_FALSE(deformInv!=NULL);

            hoImageRegDeformationFieldBidirectionalRegister<TargetType, CoordType> reg(resolution_pyramid_levels_, use_world_coordinates_, bg_value_);

            if ( !debugFolder_.empty() )
            {
                reg.debugFolder_ = debugFolder_;
            }

            GADGET_CHECK_RETURN_FALSE(reg.setDefaultParameters(resolution_pyramid_levels_, use_world_coordinates_));

            reg.max_iter_num_pyramid_level_ = max_iter_num_pyramid_level_;
            reg.div_num_pyramid_level_ = div_num_pyramid_level_;
            reg.dissimilarity_MI_betaArg_ = dissimilarity_MI_betaArg_;
            reg.regularization_hilbert_strength_world_coordinate_ = regularization_hilbert_strength_world_coordinate_;
            reg.regularization_hilbert_strength_pyramid_level_ = regularization_hilbert_strength_pyramid_level_;
            reg.dissimilarity_LocalCCR_sigmaArg_ = dissimilarity_LocalCCR_sigmaArg_;
            reg.boundary_handler_type_warper_ = boundary_handler_type_warper_;
            reg.interp_type_warper_ = interp_type_warper_;
            reg.inverse_deform_enforce_iter_pyramid_level_ = inverse_deform_enforce_iter_pyramid_level_;
            reg.inverse_deform_enforce_weight_pyramid_level_ = inverse_deform_enforce_weight_pyramid_level_;
            reg.apply_in_FOV_constraint_ = apply_in_FOV_constraint_;
            reg.apply_divergence_free_constraint_ = apply_divergence_free_constraint_;

            reg.verbose_ = verbose_;

            reg.dissimilarity_type_.clear();
            reg.dissimilarity_type_.resize(resolution_pyramid_levels_, dissimilarity_type_);

            reg.setTarget( const_cast<TargetType&>(target) );
            reg.setSource( const_cast<SourceType&>(source) );

            if ( verbose_ )
            {
                Gadgetron::printInfo(reg);
            }

            GADGET_CHECK_RETURN_FALSE(reg.initialize());

            unsigned int d;

            if ( target.dimensions_equal( *(deform[0]) ) )
            {
                if ( initial )
                {
                    for ( d=0; d<DIn; d++ )
                    {
                        reg.transform_->setDeformationField( *(deform[d]), d);
                        reg.transform_inverse_->setDeformationField( *(deformInv[d]), d);
                    }
                }
            }
            else
            {
                for ( d=0; d<DIn; d++ )
                {
                    deform[d]->copyImageInfo(target);
                    Gadgetron::clear( *(deform[d]) );
                    deformInv[d]->copyImageInfo(target);
                    Gadgetron::clear( *(deformInv[d]) );
                }
            }

            GADGET_CHECK_RETURN_FALSE(reg.performRegistration());

            for ( d=0; d<DIn; d++ )
            {
                *(deform[d]) = reg.transform_->getDeformationField(d);
                *(deformInv[d]) = reg.transform_inverse_->getDeformationField(d);
            }

            if ( warped != NULL )
            {
                /// bspline warp
                // hoNDBoundaryHandlerFixedValue<SourceType> bhFixedValue;
                hoNDBoundaryHandlerBorderValue<SourceType> bhFixedValue;
                bhFixedValue.setArray(const_cast<SourceType&>(source));

                hoNDInterpolatorBSpline<SourceType, DIn> interpBSpline(5);
                interpBSpline.setArray(const_cast<SourceType&>(source));
                interpBSpline.setBoundaryHandler(bhFixedValue);

                hoImageRegWarper<TargetType, SourceType, CoordType> warper;
                warper.setBackgroundValue(bg_value_);
                warper.setTransformation(*reg.transform_);
                warper.setInterpolator(interpBSpline);

                GADGET_CHECK_RETURN_FALSE(warper.warp(target, source, use_world_coordinates_, *warped));
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::registerTwoImagesDeformationFieldBidirectional(...) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::
    initialize(const TargetContinerType& targetContainer, bool warped)
    {
        try
        {
            if ( warped )
            {
                GADGET_CHECK_RETURN_FALSE(warped_container_.copyFrom(targetContainer));
            }

            std::vector<size_t> col = targetContainer.cols();

            unsigned int ii;

            if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD )
            {
                for ( ii=0; ii<DIn; ii++ )
                {
                    GADGET_CHECK_RETURN_FALSE(deformation_field_[ii].create(col));
                    GADGET_CHECK_RETURN_FALSE(deformation_field_[ii].fillWithZeros());
                }
            }
            else if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL )
            {
                for ( ii=0; ii<DIn; ii++ )
                {
                    GADGET_CHECK_RETURN_FALSE(deformation_field_[ii].create(col));
                    GADGET_CHECK_RETURN_FALSE(deformation_field_[ii].fillWithZeros());

                    GADGET_CHECK_RETURN_FALSE(deformation_field_inverse_[ii].create(col));
                    GADGET_CHECK_RETURN_FALSE(deformation_field_inverse_[ii].fillWithZeros());
                }
            }
            else if ( container_reg_transformation_==GT_IMAGE_REG_TRANSFORMATION_RIGID 
                        || container_reg_transformation_==GT_IMAGE_REG_TRANSFORMATION_AFFINE )
            {
                GDEBUG_STREAM("To be implemented ...");
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::initialize(const TargetContinerType& targetContainer) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::
    registerOverContainer2DPairWise(TargetContinerType& targetContainer, SourceContinerType& sourceContainer, bool warped, bool initial)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->initialize(targetContainer, warped));

            std::vector<TargetType*> targetImages;
            targetContainer.get_all_images(targetImages);

            std::vector<SourceType*> sourceImages;
            sourceContainer.get_all_images(sourceImages);

            long long numOfImages = targetImages.size();

            GADGET_CHECK_RETURN_FALSE(numOfImages==sourceImages.size());

            std::vector<SourceType*> warpedImages(numOfImages, NULL);
            if ( warped )
            {
                warped_container_.get_all_images(warpedImages);
            }

            GDEBUG_STREAM("registerOverContainer2DPairWise - threading ... ");

            int numOfThreads = 1;

#ifdef USE_OMP
            int numOfProcs = omp_get_num_procs();
            int nested = omp_get_nested();
            GDEBUG_STREAM("registerOverContainer2DPairWise - nested openMP is " << nested);
            /*if (numOfImages < numOfProcs - 1)
            {
                omp_set_nested(1);
                GDEBUG_STREAM("registerOverContainer2DPairWise - nested openMP on ... ");
            }
            else
            {
                omp_set_nested(0);
                GDEBUG_STREAM("registerOverContainer2DPairWise - nested openMP off ... ");
            }*/

            numOfThreads = (numOfImages>numOfProcs) ? numOfProcs : numOfImages;
#endif // USE_OMP

            unsigned int ii;
            long long n;

            if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD )
            {
                std::vector< std::vector<DeformationFieldType*> > deform(DIn);

                for ( ii=0; ii<DIn; ii++ )
                {
                    deformation_field_[ii].get_all_images(deform[ii]);
                }

                #pragma omp parallel default(none) private(n, ii) shared(numOfImages, initial, targetImages, sourceImages, deform, warpedImages) num_threads(numOfThreads)
                {
                    DeformationFieldType* deformCurr[DIn];

                    #pragma omp for 
                    for ( n=0; n<numOfImages; n++ )
                    {
                        TargetType& target = *(targetImages[n]);
                        SourceType& source = *(sourceImages[n]);

                        if ( &target == &source )
                        {
                            for ( ii=0; ii<DIn; ii++ )
                            {
                                deform[ii][n]->create(*target.get_dimensions());
                                Gadgetron::clear( *deform[ii][n] );
                            }
                        }
                        else
                        {
                            for ( ii=0; ii<DIn; ii++ )
                            {
                                deformCurr[ii] = deform[ii][n];
                            }

                            registerTwoImagesDeformationField(target, source, initial, warpedImages[n], deformCurr);
                        }
                    }
                }
            }
            else if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL )
            {
                std::vector< std::vector<DeformationFieldType*> > deform(DIn);
                std::vector< std::vector<DeformationFieldType*> > deformInv(DIn);

                for ( ii=0; ii<DIn; ii++ )
                {
                    deformation_field_[ii].get_all_images(deform[ii]);
                    deformation_field_inverse_[ii].get_all_images(deformInv[ii]);
                }

                #pragma omp parallel default(none) private(n, ii) shared(numOfImages, initial, targetImages, sourceImages, deform, deformInv, warpedImages) num_threads(numOfThreads)
                {
                    DeformationFieldType* deformCurr[DIn];
                    DeformationFieldType* deformInvCurr[DIn];

                    #pragma omp for 
                    for ( n=0; n<numOfImages; n++ )
                    {
                        TargetType& target = *(targetImages[n]);
                        SourceType& source = *(sourceImages[n]);

                        if ( &target == &source )
                        {
                            for ( ii=0; ii<DIn; ii++ )
                            {
                                deform[ii][n]->create(*target.get_dimensions());
                                Gadgetron::clear( *deform[ii][n] );

                                deformInv[ii][n]->create(*source.get_dimensions());
                                Gadgetron::clear( *deformInv[ii][n] );
                            }
                        }
                        else
                        {
                            for ( ii=0; ii<DIn; ii++ )
                            {
                                deformCurr[ii] = deform[ii][n];
                                deformInvCurr[ii] = deformInv[ii][n];
                            }

                            registerTwoImagesDeformationFieldBidirectional(target, source, initial, warpedImages[n], deformCurr, deformInvCurr);
                        }
                    }
                }
            }
            else if ( container_reg_transformation_==GT_IMAGE_REG_TRANSFORMATION_RIGID 
                        || container_reg_transformation_==GT_IMAGE_REG_TRANSFORMATION_AFFINE )
            {
                GDEBUG_STREAM("To be implemented ...");
            }

//#ifdef USE_OMP
//            omp_set_nested(nested);
//#endif // USE_OMP
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::registerOverContainer2DPairWise(...) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::
    registerOverContainer2DFixedReference(TargetContinerType& imageContainer, const std::vector<unsigned int>& referenceFrame, bool warped, bool initial)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->initialize(imageContainer, warped));

            size_t row = imageContainer.rows();
            std::vector<size_t> col = imageContainer.cols();

            GADGET_CHECK_RETURN_FALSE(referenceFrame.size() == col.size());

            std::vector<SourceType*> sourceImages;
            imageContainer.get_all_images(sourceImages);

            long long numOfImages = (long long)sourceImages.size();

            // warped images
            std::vector<SourceType*> warpedImages(numOfImages, NULL);
            if ( warped )
            {
                warped_container_.get_all_images(warpedImages);
            }

            unsigned int ii;
            long long n;
            size_t r, c;

            // fill in the reference frames
            std::vector<TargetType*> targetImages(numOfImages, NULL);

            size_t ind=0;
            for ( r=0; r<row; r++ )
            {
                TargetType& ref = const_cast<TargetType&>(imageContainer(r, referenceFrame[r]));

                for ( c=0; c<col[r]; c++ )
                {
                    targetImages[ind] = &ref;
                    ind++;
                }
            }

            GADGET_CHECK_RETURN_FALSE(numOfImages==targetImages.size());

            int numOfThreads = 1;

#ifdef USE_OMP
            int numOfProcs = omp_get_num_procs();
            int nested = omp_get_nested();
            GDEBUG_STREAM("registerOverContainer2DFixedReference - nested openMP is " << nested);
            //if (numOfImages < numOfProcs - 1)
            //{
            //    omp_set_nested(1);
            //    GDEBUG_STREAM("registerOverContainer2DFixedReference - nested openMP on ... ");
            //}
            //else
            //{
            //    omp_set_nested(0);
            //    GDEBUG_STREAM("registerOverContainer2DFixedReference - nested openMP off ... ");
            //}

            numOfThreads = (numOfImages>numOfProcs) ? numOfProcs : numOfImages;
#endif // USE_OMP

            if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD )
            {
                std::vector< std::vector<DeformationFieldType*> > deform(DIn);

                for ( ii=0; ii<DIn; ii++ )
                {
                    deformation_field_[ii].get_all_images(deform[ii]);
                }

                #pragma omp parallel default(none) private(n, ii) shared(numOfImages, initial, targetImages, sourceImages, deform, warpedImages) num_threads(numOfThreads)
                {
                    DeformationFieldType* deformCurr[DIn];

                    #pragma omp for 
                    for ( n=0; n<numOfImages; n++ )
                    {
                        if ( targetImages[n] == sourceImages[n] )
                        {
                            if ( warpedImages[n] != NULL )
                            {
                                *(warpedImages[n]) = *(targetImages[n]);
                            }

                            for ( ii=0; ii<DIn; ii++ )
                            {
                                deform[ii][n]->create(*targetImages[n]->get_dimensions());
                                Gadgetron::clear(*deform[ii][n]);
                            }

                            continue;
                        }

                        TargetType& target = *(targetImages[n]);
                        SourceType& source = *(sourceImages[n]);

                        for ( ii=0; ii<DIn; ii++ )
                        {
                            deformCurr[ii] = deform[ii][n];
                        }

                        registerTwoImagesDeformationField(target, source, initial, warpedImages[n], deformCurr);
                    }
                }
            }
            else if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL )
            {
                std::vector< std::vector<DeformationFieldType*> > deform(DIn);
                std::vector< std::vector<DeformationFieldType*> > deformInv(DIn);

                for ( ii=0; ii<DIn; ii++ )
                {
                    deformation_field_[ii].get_all_images(deform[ii]);
                    deformation_field_inverse_[ii].get_all_images(deformInv[ii]);
                }

                #pragma omp parallel default(none) private(n, ii) shared(numOfImages, initial, targetImages, sourceImages, deform, deformInv, warpedImages) num_threads(numOfThreads)
                {
                    DeformationFieldType* deformCurr[DIn];
                    DeformationFieldType* deformInvCurr[DIn];

                    #pragma omp for 
                    for ( n=0; n<numOfImages; n++ )
                    {
                        if ( targetImages[n] == sourceImages[n] )
                        {
                            if ( warpedImages[n] != NULL )
                            {
                                *(warpedImages[n]) = *(targetImages[n]);
                            }

                            for ( ii=0; ii<DIn; ii++ )
                            {
                                deform[ii][n]->create(*targetImages[n]->get_dimensions());
                                Gadgetron::clear(*deform[ii][n]);

                                deformInv[ii][n]->create(*targetImages[n]->get_dimensions());
                                Gadgetron::clear(*deformInv[ii][n]);
                            }

                            continue;
                        }

                        TargetType& target = *(targetImages[n]);
                        SourceType& source = *(sourceImages[n]);

                        for ( ii=0; ii<DIn; ii++ )
                        {
                            deformCurr[ii] = deform[ii][n];
                            deformInvCurr[ii] = deformInv[ii][n];
                        }

                        registerTwoImagesDeformationFieldBidirectional(target, source, initial, warpedImages[n], deformCurr, deformInvCurr);
                    }
                }
            }
            else if ( container_reg_transformation_==GT_IMAGE_REG_TRANSFORMATION_RIGID 
                        || container_reg_transformation_==GT_IMAGE_REG_TRANSFORMATION_AFFINE )
            {
                GDEBUG_STREAM("To be implemented ...");
            }

//#ifdef USE_OMP
//            omp_set_nested(nested);
//#endif // USE_OMP
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::registerOverContainer2DFixedReference(...) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::
    registerOverContainer2DProgressive(TargetContinerType& imageContainer, const std::vector<unsigned int>& referenceFrame)
    {
        try
        {
            bool warped = true;
            GADGET_CHECK_RETURN_FALSE(this->initialize(imageContainer, warped));

            long long row = (long long)imageContainer.rows();
            std::vector<size_t> col = imageContainer.cols();

            GADGET_CHECK_RETURN_FALSE(referenceFrame.size() == col.size());

            unsigned int ii;
            long long n;
            long long r, c;

            // for every row, two registration tasks can be formatted

            long long numOfTasks = (long long)(2*row);
            GDEBUG_STREAM("hoImageRegContainer2DRegistration<...>::registerOverContainer2DProgressive(...), numOfTasks : " << numOfTasks);

            std::vector< std::vector<TargetType*> > regImages(numOfTasks);
            std::vector< std::vector<TargetType*> > warpedImages(numOfTasks);

            std::vector< std::vector< std::vector<DeformationFieldType*> > > deform(DIn);
            std::vector< std::vector< std::vector<DeformationFieldType*> > > deformInv(DIn);

            for ( ii=0; ii<DIn; ii++ )
            {
                deform[ii].resize(numOfTasks);
                deformInv[ii].resize(numOfTasks);
            }

            for ( r=0; r<row; r++ )
            {
                unsigned int refFrame = referenceFrame[r];

                regImages[2*r].resize(col[r]-refFrame);
                regImages[2*r+1].resize(1+refFrame);

                warpedImages[2*r].resize(col[r]-refFrame);
                warpedImages[2*r+1].resize(1+refFrame);

                // copy over the reference frame
                warped_container_(r, refFrame) = imageContainer(r, refFrame);

                if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD )
                {
                    for ( ii=0; ii<DIn; ii++ )
                    {
                        deformation_field_[ii](r, refFrame).create(*imageContainer(r, refFrame).get_dimensions());
                        Gadgetron::clear(deformation_field_[ii](r, refFrame));
                    }
                }

                if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL )
                {
                    for ( ii=0; ii<DIn; ii++ )
                    {
                        deformation_field_[ii](r, refFrame).create(*imageContainer(r, refFrame).get_dimensions());
                        Gadgetron::clear(deformation_field_[ii](r, refFrame));

                        deformation_field_inverse_[ii](r, refFrame).create(*imageContainer(r, refFrame).get_dimensions());
                        Gadgetron::clear(deformation_field_inverse_[ii](r, refFrame));
                    }
                }

                // task one
                for ( c=refFrame; c<(long long)col[r]; c++ )
                {
                    regImages[2*r][c-refFrame] = &(imageContainer(r, c));
                    warpedImages[2*r][c-refFrame] = &(warped_container_(r, c));
                }

                // task two
                for ( c=refFrame; c>=0; c-- )
                {
                    regImages[2*r+1][refFrame-c] = &(imageContainer(r, c));
                    warpedImages[2*r+1][refFrame-c] = &(warped_container_(r, c));
                }

                for ( ii=0; ii<DIn; ii++ )
                {
                    if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD )
                    {
                        deform[ii][2*r].resize(col[r]-refFrame);
                        deform[ii][2*r+1].resize(1+refFrame);

                        // task one
                        for ( c=refFrame; c<(long long)col[r]; c++ )
                        {
                            deform[ii][2*r][c-refFrame] = &(deformation_field_[ii](r, c));
                        }

                        // task two
                        for ( c=refFrame; c>=0; c-- )
                        {
                            deform[ii][2*r+1][refFrame-c] = &(deformation_field_[ii](r, c));
                        }
                    }

                    if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL )
                    {
                        deform[ii][2*r].resize(col[r]-refFrame);
                        deform[ii][2*r+1].resize(1+refFrame);

                        deformInv[ii][2*r].resize(col[r]-refFrame);
                        deformInv[ii][2*r+1].resize(1+refFrame);

                        // task one
                        for ( c=refFrame; c<(long long)col[r]; c++ )
                        {
                            deform[ii][2*r][c-refFrame] = &(deformation_field_[ii](r, c));
                            deformInv[ii][2*r][c-refFrame] = &(deformation_field_inverse_[ii](r, c));
                        }

                        // task two
                        for ( c=refFrame; c>=0; c-- )
                        {
                            deform[ii][2*r+1][refFrame-c] = &(deformation_field_[ii](r, c));
                            deformInv[ii][2*r+1][refFrame-c] = &(deformation_field_inverse_[ii](r, c));
                        }
                    }
                }
            }

            if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD )
            {
                bool initial = false;

                #pragma omp parallel default(none) private(n, ii) shared(numOfTasks, initial, regImages, warpedImages, deform)
                {
                    DeformationFieldType* deformCurr[DIn];

                    #pragma omp for 
                    for ( n=0; n<numOfTasks; n++ )
                    {
                        size_t numOfImages = regImages[n].size();

                        // no need to copy the refrence frame to warped

                        size_t k;
                        for ( k=1; k<numOfImages; k++ )
                        {
                            TargetType& target = *(warpedImages[n][k-1]);
                            SourceType& source = *(regImages[n][k]);

                            for ( ii=0; ii<DIn; ii++ )
                            {
                                deformCurr[ii] = deform[ii][n][k];
                            }

                            registerTwoImagesDeformationField(target, source, initial, warpedImages[n][k], deformCurr);
                        }
                    }
                }
            }
            else if ( container_reg_transformation_ == GT_IMAGE_REG_TRANSFORMATION_DEFORMATION_FIELD_BIDIRECTIONAL )
            {
                bool initial = false;

                #pragma omp parallel default(none) private(n, ii) shared(numOfTasks, initial, regImages, warpedImages, deform, deformInv)
                {
                    DeformationFieldType* deformCurr[DIn];
                    DeformationFieldType* deformInvCurr[DIn];

                    #pragma omp for 
                    for ( n=0; n<numOfTasks; n++ )
                    {
                        size_t numOfImages = regImages[n].size();

                        size_t k;
                        for ( k=1; k<numOfImages; k++ )
                        {
                            TargetType& target = *(warpedImages[n][k-1]);
                            SourceType& source = *(regImages[n][k]);

                            for ( ii=0; ii<DIn; ii++ )
                            {
                                deformCurr[ii] = deform[ii][n][k];
                                deformInvCurr[ii] = deformInv[ii][n][k];
                            }

                            registerTwoImagesDeformationFieldBidirectional(target, source, initial, warpedImages[n][k], deformCurr, deformInvCurr);
                        }
                    }
                }
            }
            else if ( container_reg_transformation_==GT_IMAGE_REG_TRANSFORMATION_RIGID 
                        || container_reg_transformation_==GT_IMAGE_REG_TRANSFORMATION_AFFINE )
            {
                GDEBUG_STREAM("To be implemented ...");
            }
        }
        catch(...)
        {
            GERROR_STREAM("Error happened in hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::registerOverContainer2DProgressive(...) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegContainer2DRegistration<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;

        unsigned int ii, jj;

        os << "--------------Gagdgetron image registration container 2D -------------" << endl;

        os << "Input dimension is : " << DIn << endl;
        os << "Output dimension is : " << DOut << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Image data type is : " << elemTypeName << std::endl;

        elemTypeName = std::string(typeid(CoordType).name());
        os << "Transformation coordinate data type is : " << elemTypeName << std::endl;

        os << "Whether to apply in_FOV constraint : " << apply_in_FOV_constraint_ << std::endl;
        os << "Whether to apply divergence free constraint : " << apply_divergence_free_constraint_ << std::endl;
        os << "Whether to perform world coordinate registration is : " << use_world_coordinates_ << std::endl;
        os << "Number of resolution pyramid levels is : " << resolution_pyramid_levels_ << std::endl;

        os << "------------" << std::endl;
        os << "Number of iterations is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << max_iter_num_pyramid_level_[ii] << std::endl;
        }

        os << "------------" << std::endl;
        os << "Image dissimilarity is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << getDissimilarityName(dissimilarity_type_) << std::endl;
        }

        os << "------------" << std::endl;
        os << "Threshold for dissimilarity is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << dissimilarity_thres_pyramid_level_[ii] << std::endl;
        }

        os << "------------" << std::endl;
        os << "Number of search size division is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << div_num_pyramid_level_[ii] << std::endl;
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
            for( jj=0; jj<DIn; jj++ )
            {
                os << regularization_hilbert_strength_pyramid_level_[ii][jj] << " ";
            } 
            os << " ] " << std::endl;
        }

        os << "------------" << std::endl;
        os << "Boundary handler and interpolator type for warper is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << getBoundaryHandlerName(boundary_handler_type_warper_[ii]) 
                << " - " << getInterpolatorName(interp_type_warper_[ii]) << std::endl;
        }

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
        os << "------------" << std::endl;
    }
}

#endif // hoImageRegContainer2DRegistration_H_
