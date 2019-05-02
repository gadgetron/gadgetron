/** \file   hoImageRegParametricSolver.h
    \brief  Define the base class of image registration solver for parametric image transformation
    \author Hui Xue
*/

#ifndef hoImageRegParametricSolver_H_
#define hoImageRegParametricSolver_H_

#pragma once

#include "hoImageRegSolver.h"

namespace Gadgetron {

    /// ValueType: image pixel value type
    /// CoordType: transformation data type
    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegParametricSolver : public hoImageRegSolver<TargetType, SourceType, CoordType>
    {
    public:

        typedef hoImageRegParametricSolver<TargetType, SourceType, CoordType> Self;
        typedef hoImageRegSolver<TargetType, SourceType, CoordType> BaseClass;

        typedef typename TargetType::value_type ValueType;
        enum { DIn = TargetType::NDIM };
        enum { DOut = SourceType::NDIM };

        typedef hoNDImage<ValueType, 2> Target2DType;
        typedef Target2DType Source2DType;

        typedef hoNDImage<ValueType, 3> Target3DType;
        typedef Target2DType Source3DType;

        typedef typename BaseClass::InterpolatorType InterpolatorType;

        typedef hoImageRegParametricTransformation<CoordType, DIn, DOut> TransformationType;
        typedef typename TransformationType::ParaStatus ParaStatusType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef CoordType coord_type;

        typedef typename TransformationType::input_point_type input_point_type;
        typedef typename TransformationType::output_point_type output_point_type;

        typedef typename TransformationType::jacobian_parameter_type jacobian_parameter_type;
        typedef typename TransformationType::jacobian_position_type jacobian_position_type;

        typedef typename BaseClass::ImageRegWarperType ImageRegWarperType;

        typedef typename BaseClass::ImageRegDissimilarityType ImageRegDissimilarityType;

        hoImageRegParametricSolver();
        virtual ~hoImageRegParametricSolver();

        void setTransform(TransformationType& transform) { transform_ = &transform; }

        virtual bool initialize();

        /// solve the minimization and find the optimal transformation
        virtual bool solve();

        /// perform one iteration of optimization
        virtual ValueType solver_once(ValueType curr_dissimilarity) = 0;

        /// compute the derivatives of dissimilarity measures to the transformation parameters
        /// if the analytic derivative is hard to compute, the central difference derivative is computed
        /// deriv_step_size is the step size used to compute central difference derivatives
        virtual bool evaluateDeriv(TransformationType* transform, ImageRegDissimilarityType* dissimilarity, const std::vector<ValueType>& deriv_step_size, std::vector<ValueType>& deriv);

        virtual void print(std::ostream& os) const;

        /// number of performed iterations
        unsigned int iter_num_;

        /// maximal number of iterations
        unsigned int max_iter_num_;

        /// threshold for minimal dissimilarity changes
        ValueType dissimilarity_thres_;

        /// threshold for minimal parameter changes
        ValueType parameter_thres_;

        /// number of search division
        unsigned int div_num_;

        /// step size for every parameters used in optimization
        /// depending on the optimization algorithm, this variable may not be used
        std::vector<ValueType> step_size_para_;
        /// step size division ratio
        /// step_size_para_ = step_size_para_ .* step_size_div_para_ to reduce search step size
        std::vector<ValueType> step_size_div_para_;

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
    hoImageRegParametricSolver<TargetType, SourceType, CoordType>::hoImageRegParametricSolver() 
        : BaseClass(), dissimilarity_thres_(1e-8), parameter_thres_(1e-8)
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegParametricSolver<TargetType, SourceType, CoordType>::~hoImageRegParametricSolver()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegParametricSolver<TargetType, SourceType, CoordType>::initialize()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(target_!=NULL);
            GADGET_CHECK_RETURN_FALSE(source_!=NULL);
            GADGET_CHECK_RETURN_FALSE(interp_!=NULL);
            GADGET_CHECK_RETURN_FALSE(warper_!=NULL);
            GADGET_CHECK_RETURN_FALSE(dissimilarity_!=NULL);
            GADGET_CHECK_RETURN_FALSE(transform_!=NULL);

            warper_->setTransformation(*transform_);
            warper_->setInterpolator(*interp_);
            warper_->setBackgroundValue(bg_value_);

            dissimilarity_->setBackgroundValue(bg_value_);

            if ( !warpped_.dimensions_equal(*target_) )
            {
                warpped_ = *target_;
            }

            dissimilarity_->initialize(*target_);

            if ( step_size_para_.size() != transform_->get_number_of_parameters() )
            {
                step_size_para_.resize(transform_->get_number_of_parameters(), 1.0);
            }

            if ( step_size_div_para_.size() != transform_->get_number_of_parameters() )
            {
                step_size_div_para_.resize(transform_->get_number_of_parameters(), 0.5);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegParametricSolver<TargetType, SourceType, CoordType>::initialize() ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegParametricSolver<TargetType, SourceType, CoordType>::solve()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->initialize());

            size_t numOfPara = transform_->get_number_of_parameters();

            GADGET_CHECK_RETURN_FALSE(warper_->warp(*target_, *source_, use_world_coordinate_, warpped_));
            curr_dissimilarity_ = dissimilarity_->evaluate(warpped_);

            if ( verbose_ ) { GDEBUG_STREAM("----> Initial image dissimilarity : " << curr_dissimilarity_); }

            unsigned int totalIterNum = 0;

            unsigned int div;
            for ( div=0; div<div_num_; div++ )
            {
                if ( verbose_ ) { GDEBUG_STREAM("----> Parameter division " << div << " [out of " << div_num_ << "] "); }

                for ( iter_num_=0; iter_num_<max_iter_num_; iter_num_++ )
                {
                    if ( verbose_ ) { GDEBUG_STREAM("--> Iteration " << iter_num_ << " [out of " << max_iter_num_ << "] : \t" << curr_dissimilarity_); }

                    prev_dissimilarity_ = curr_dissimilarity_;

                    curr_dissimilarity_ = this->solver_once(prev_dissimilarity_);

                    // if the dissimilarity stops decreasing
                    if ( prev_dissimilarity_ < curr_dissimilarity_ + dissimilarity_thres_ )
                    {
                        break;
                    }
                }
                if ( verbose_ ) { transform_->printTransform(std::cout); }

                totalIterNum += iter_num_;

                // reduce the step size
                size_t p;
                for ( p=0; p<numOfPara; p++ )
                {
                    step_size_para_[p] *= step_size_div_para_[p];
                }
            }

            if ( verbose_ ) { GDEBUG_STREAM("----> Total iteration number : " << totalIterNum); }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegParametricSolver<TargetType, SourceType, CoordType>::solve() ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegParametricSolver<TargetType, SourceType, CoordType>::
    evaluateDeriv(TransformationType* transform, ImageRegDissimilarityType* dissimilarity, const std::vector<ValueType>& deriv_step_size, std::vector<ValueType>& deriv)
    {
        try
        {
            bool has_analytic_deriv = false;

            /// for some transformation and dissimilarity combination, the analytical derivative is easier to be computed

            /// implement the central difference numerical derivative
            if ( !has_analytic_deriv )
            {
                size_t numOfPara = transform_->get_number_of_parameters();
                size_t i;

                deriv.resize(numOfPara, 0);

                ValueType currPara(0), positiveValue(0), negativeValue(0), normDeriv(0);
                for ( i=0; i<numOfPara; i++ )
                {
                    if ( transform_->get_para_status(i) == TransformationType::Active )
                    {
                        currPara = transform_->get_parameter(i);

                        // positive
                        transform_->set_parameter(i, currPara + deriv_step_size[i]);

                        GADGET_CHECK_RETURN_FALSE(warper_->warp(*target_, *source_, use_world_coordinate_, warpped_));
                        positiveValue = dissimilarity_->evaluate(warpped_);

                        // negative
                        transform_->set_parameter(i, currPara - deriv_step_size[i]);

                        GADGET_CHECK_RETURN_FALSE(warper_->warp(*target_, *source_, use_world_coordinate_, warpped_));
                        negativeValue = dissimilarity_->evaluate(warpped_);

                        deriv[i] = (positiveValue - negativeValue)/(2*deriv_step_size[i]);
                        normDeriv += deriv[i]*deriv[i];

                        transform_->set_parameter(i, currPara);
                    }
                }

                if ( normDeriv > 0 )
                {
                    ValueType distDeriv=std::sqrt(normDeriv);

                    for ( i=0; i<numOfPara; i++ )
                    {
                        deriv[i] /= distDeriv;
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in evaluateDeriv(TransformationType* transform, ImageRegDissimilarityType* dissimilarity, const std::vector<ValueType>& deriv_step_size, std::vector<ValueType>& deriv) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegParametricSolver<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron parametric image registration solver -------------" << endl;
        os << "Target image dimension is : " << DIn << endl;
        os << "Source image dimension is : " << DOut << endl;
        os << "Image data type is : " << std::string(typeid(ValueType).name()) << std::endl;
        os << "Transformation data type is : " << std::string(typeid(CoordType).name()) << std::endl;
        os << "Use world coordinate is : " << use_world_coordinate_ << std::endl;
        os << "Maximal iteration number is : " << max_iter_num_ << std::endl;
        os << "Dissimilarity threshold is : " << dissimilarity_thres_ << std::endl;
        os << "Parameter threshold is : " << parameter_thres_ << std::endl;
        os << "Number of search division is : " << div_num_ << std::endl;

        os << "Step size for every parameters used in optimization is : [ ";
        unsigned int ii;
        for ( ii=0; ii<step_size_para_.size(); ii++ )
        {
            os << step_size_para_[ii] << " ";
        }
        os << " ] " << endl;

        os << "Step size division ratio is : [ ";
        for ( ii=0; ii<step_size_div_para_.size(); ii++ )
        {
            os << step_size_div_para_[ii] << " ";
        }
        os << " ] " << endl;
    }
}
#endif // hoImageRegParametricSolver_H_
