/** \file   hoImageRegParametricRegister.h
    \brief  Define the class to perform parametric image registration in gadgetron
            By default, the multi-level multi-step parametric solver is used
    \author Hui Xue
*/

#ifndef hoImageRegParametricRegister_H_
#define hoImageRegParametricRegister_H_

#pragma once

#include "hoImageRegRegister.h"

namespace Gadgetron {

    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegParametricRegister : public hoImageRegRegister<TargetType, SourceType, CoordType>
    {
    public:

        typedef hoImageRegParametricRegister<TargetType, SourceType, CoordType> Self;
        typedef hoImageRegRegister<TargetType, SourceType, CoordType> BaseClass;

        typedef typename TargetType::value_type ValueType;
        enum { DIn = TargetType::NDIM };
        enum { DOut = SourceType::NDIM };

        typedef typename BaseClass::Target2DType Target2DType;
        typedef typename BaseClass::Source2DType Source2DType;

        typedef typename BaseClass::Target3DType Target3DType;
        typedef typename BaseClass::Source3DType Source3DType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef CoordType coord_type;

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
        typedef hoImageRegParametricTransformation<CoordType, DIn, DOut> TransformationType;
        typedef typename TransformationType::input_point_type input_point_type;
        typedef typename TransformationType::output_point_type output_point_type;
        typedef typename TransformationType::jacobian_position_type jacobian_position_type;

        /// solver type
        typedef hoImageRegParametricSolver<TargetType, SourceType, CoordType> SolverType;

        hoImageRegParametricRegister(unsigned int resolution_pyramid_levels=3, bool use_world_coordinates=true, ValueType bg_value=ValueType(0));
        virtual ~hoImageRegParametricRegister();

        /// initialize the registration
        /// should be called after all images and parameters of registration are set
        virtual bool initialize();

        /// create parametric solver
        SolverType* createParametricSolver(GT_IMAGE_REG_SOLVER v, unsigned int level);

        /// set the default parameters
        virtual bool setDefaultParameters(unsigned int resolution_pyramid_levels, bool use_world_coordinates);

        /// perform the registration
        virtual bool performRegistration();

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

        using BaseClass::dissimilarity_LocalCCR_sigmaArg_;
        using BaseClass::dissimilarity_hist_num_bin_target_;
        using BaseClass::dissimilarity_hist_num_bin_warpped_;
        using BaseClass::dissimilarity_hist_pv_interpolation_;
        using BaseClass::dissimilarity_hist_step_size_ignore_pixel_;

        using BaseClass::dissimilarity_MI_betaArg_;

        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

        /// verbose mode
        bool verbose_;

        /// deformation field transformation, defined in the world coordinate of target image
        TransformationType* transform_;

        /// solver
        std::vector<SolverType*> solver_pyramid_;

        /// for solver of every pyramid level

        /// maximal number of iterations
        std::vector<unsigned int> max_iter_num_pyramid_level_;

        /// threshold for minimal dissimilarity changes
        ValueType dissimilarity_thres_;

        /// threshold for minimal parameter changes
        ValueType parameter_thres_;

        /// number of search division
        std::vector<unsigned int> div_num_pyramid_level_;

        /// step size for every parameter
        std::vector< std::vector<ValueType> > step_size_para_pyramid_level_;

        /// step size division ratio
        /// step_size_para_ = step_size_para_ .* step_size_div_para_ to reduce search step size
        std::vector< std::vector<ValueType> > step_size_div_para_pyramid_level_;

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

        bool preset_transform_;
    };

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegParametricRegister<TargetType, SourceType, CoordType>::
    hoImageRegParametricRegister(unsigned int resolution_pyramid_levels, bool use_world_coordinates, ValueType bg_value) : BaseClass(resolution_pyramid_levels, bg_value), verbose_(false), preset_transform_(false)
    {
        GADGET_CHECK_THROW(this->setDefaultParameters(resolution_pyramid_levels, use_world_coordinates));
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegParametricRegister<TargetType, SourceType, CoordType>::~hoImageRegParametricRegister()
    {

    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegParametricRegister<TargetType, SourceType, CoordType>::setDefaultParameters(unsigned int resolution_pyramid_levels, bool use_world_coordinates)
    {
        use_world_coordinates_ = use_world_coordinates;
        resolution_pyramid_levels_ = resolution_pyramid_levels;

        resolution_pyramid_downsample_ratio_.clear();
        resolution_pyramid_downsample_ratio_.resize(resolution_pyramid_levels_-1, std::vector<float>(std::max( (int)DIn, (int)DOut), 2.0) );

        resolution_pyramid_blurring_sigma_.clear();
        resolution_pyramid_blurring_sigma_.resize(resolution_pyramid_levels_, std::vector<float>(std::max( (int)DIn, (int)DOut), 0.0) );

        boundary_handler_type_warper_.clear();
        boundary_handler_type_warper_.resize(resolution_pyramid_levels_, GT_BOUNDARY_CONDITION_FIXEDVALUE);

        interp_type_warper_.clear();
        interp_type_warper_.resize(resolution_pyramid_levels_, GT_IMAGE_INTERPOLATOR_LINEAR);

        boundary_handler_type_pyramid_construction_ = GT_BOUNDARY_CONDITION_BORDERVALUE;
        interp_type_pyramid_construction_ = GT_IMAGE_INTERPOLATOR_LINEAR;

        dissimilarity_type_.clear();
        dissimilarity_type_.resize(resolution_pyramid_levels_, GT_IMAGE_DISSIMILARITY_NMI);

        solver_type_.clear();
        solver_type_.resize(resolution_pyramid_levels_, GT_IMAGE_REG_SOLVER_DOWNHILL);

        max_iter_num_pyramid_level_.clear();
        max_iter_num_pyramid_level_.resize(resolution_pyramid_levels_, 100);

        dissimilarity_thres_ = 1e-6;

        div_num_pyramid_level_.clear();
        div_num_pyramid_level_.resize(resolution_pyramid_levels_, 5);

        step_size_para_pyramid_level_.clear();
        step_size_para_pyramid_level_.resize(resolution_pyramid_levels_);

        step_size_div_para_pyramid_level_.clear();
        step_size_div_para_pyramid_level_.resize(resolution_pyramid_levels_);

        size_t maxParaNum = 4096;

        step_size_para_pyramid_level_[resolution_pyramid_levels_-1].resize(maxParaNum, 3.2);
        step_size_div_para_pyramid_level_[resolution_pyramid_levels_-1].resize(maxParaNum, 0.5);

        int ii;
        unsigned int jj;
        for ( ii=(int)resolution_pyramid_levels_-2; ii>=0; ii-- )
        {
            step_size_div_para_pyramid_level_[ii].resize(maxParaNum, 0.5);
            step_size_para_pyramid_level_[ii].resize(maxParaNum);

            for ( jj=0; jj<maxParaNum; jj++ )
            {
                step_size_para_pyramid_level_[ii][jj] = step_size_div_para_pyramid_level_[ii][jj]*step_size_para_pyramid_level_[ii+1][jj];
            }
        }

        verbose_ = false;

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegParametricRegister<TargetType, SourceType, CoordType>::initialize()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(BaseClass::initialize());

            GADGET_CHECK_RETURN_FALSE( transform_ != NULL );

            warper_pyramid_.resize(resolution_pyramid_levels_);
            solver_pyramid_.resize(resolution_pyramid_levels_);

            unsigned int ii, jj;
            for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
            {
                warper_pyramid_[ii].setTransformation(*transform_);
                warper_pyramid_[ii].setInterpolator( *source_interp_warper_[ii] );
                warper_pyramid_[ii].setBackgroundValue(bg_value_);
                warper_pyramid_[ii].debugFolder_ = this->debugFolder_;

                solver_pyramid_[ii] = this->createParametricSolver(solver_type_[ii], ii);

                solver_pyramid_[ii]->setTransform(*transform_);

                solver_pyramid_[ii]->verbose_ = verbose_;
                solver_pyramid_[ii]->debugFolder_ = this->debugFolder_;

                solver_pyramid_[ii]->setTarget(target_pyramid_[ii]);
                solver_pyramid_[ii]->setSource(source_pyramid_[ii]);
                solver_pyramid_[ii]->setDissimilarity(*dissimilarity_pyramid_[ii]);
                solver_pyramid_[ii]->setWarper(warper_pyramid_[ii]);
                solver_pyramid_[ii]->setInterpolator(*source_interp_warper_[ii]);
                solver_pyramid_[ii]->setBackgroundValue(bg_value_);
                solver_pyramid_[ii]->setUseWorldCoordinate(use_world_coordinates_);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegDeformationFieldRegister<ValueType, CoordType, D>::initialize() ... ");
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegParametricRegister<TargetType, SourceType, CoordType>::performRegistration()
    {
        try
        {
            // starting from the most coarse level

            if ( verbose_ )
            {
                GDEBUG_STREAM("Initial transformation : ");
                transform_->print(std::cout);
            }

            int level;
            for ( level=(int)resolution_pyramid_levels_-1; level>=0; level-- )
            {
                // GADGET_CHECK_RETURN_FALSE(solver_pyramid_[level].initialize());
                GADGET_CHECK_RETURN_FALSE(solver_pyramid_[level]->solve());

                if ( verbose_ )
                {
                    GDEBUG_STREAM("Transformation for level " << level << " : ");
                    transform_->printTransform(std::cout);
                }

                // adjust transformation for the next resolution level
                if ( level>0 )
                {
                    if ( !use_world_coordinates_ )
                    {
                        hoMatrix<ValueType> lowResI2W, highResI2W;
                        source_pyramid_[level].image_to_world_matrix(lowResI2W);
                        source_pyramid_[level-1].image_to_world_matrix(highResI2W);

                        GADGET_CHECK_RETURN_FALSE(transform_->adjustForResolutionPyramid(lowResI2W, highResI2W));
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegParametricRegister<TargetType, SourceType, CoordType>::performRegistration() ... ");
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegParametricSolver<TargetType, SourceType, CoordType>* hoImageRegParametricRegister<TargetType, SourceType, CoordType>::createParametricSolver(GT_IMAGE_REG_SOLVER v, unsigned int level)
    {
        SolverType* res = NULL;

        unsigned int ii;

        switch (v)
        {
            case GT_IMAGE_REG_SOLVER_DOWNHILL:
                res = new hoImageRegParametricDownHillSolver<TargetType, SourceType, CoordType>();
                break;

            case GT_IMAGE_REG_SOLVER_GRADIENT_DESCENT:
                res = new hoImageRegParametricGradientDescentSolver<TargetType, SourceType, CoordType>();
                break;

            default:
                GERROR_STREAM("Unrecognized parametric solver type : " << v);
        }

        res->max_iter_num_ = max_iter_num_pyramid_level_[level];
        res->dissimilarity_thres_ = dissimilarity_thres_;
        res->parameter_thres_ = parameter_thres_;
        res->div_num_ = div_num_pyramid_level_[level];
        res->step_size_para_ = step_size_para_pyramid_level_[level];
        res->step_size_div_para_ = step_size_div_para_pyramid_level_[level];

        return res;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegParametricRegister<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron parametric image register -------------" << endl;
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
                << dissimilarity_thres_ << std::endl;
        }

        os << "------------" << std::endl;
        os << "Number of search size division for every pyramid level is : " << std::endl;
        for ( ii=0; ii<resolution_pyramid_levels_; ii++ )
        {
            os << " Level " << ii << " - " 
                << div_num_pyramid_level_[ii] << std::endl;
        }
        os << "--------------------------------------------------------------------" << endl << ends;
    }
}
#endif // hoImageRegParametricRegister_H_
