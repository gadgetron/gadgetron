/** \file   hoImageRegParametricDownHillSolver.h
    \brief  Define the class of simple down-hill solver for parametric image transformation
    \author Hui Xue
*/

#ifndef hoImageRegParametricDownHillSolver_H_
#define hoImageRegParametricDownHillSolver_H_

#pragma once

#include "hoImageRegParametricSolver.h"

namespace Gadgetron {

    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegParametricDownHillSolver : public hoImageRegParametricSolver<TargetType, SourceType, CoordType>
    {
    public:

        typedef hoImageRegParametricDownHillSolver<TargetType, SourceType, CoordType> Self;
        typedef hoImageRegParametricSolver<TargetType, SourceType, CoordType> BaseClass;

        typedef typename TargetType::value_type ValueType;
        enum { DIn = TargetType::NDIM };
        enum { DOut = SourceType::NDIM };

        typedef hoNDImage<ValueType, 2> Target2DType;
        typedef Target2DType Source2DType;

        typedef hoNDImage<ValueType, 3> Target3DType;
        typedef Target2DType Source3DType;

        typedef typename BaseClass::InterpolatorType InterpolatorType;

        typedef typename BaseClass::TransformationType TransformationType;
        typedef typename BaseClass::ParaStatusType ParaStatusType;

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

        hoImageRegParametricDownHillSolver();
        virtual ~hoImageRegParametricDownHillSolver();

        /// perform one iteration of optimization
        virtual ValueType solver_once(ValueType curr_dissimilarity);

        virtual void print(std::ostream& os) const;

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
    hoImageRegParametricDownHillSolver<TargetType, SourceType, CoordType>::hoImageRegParametricDownHillSolver() : BaseClass()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegParametricDownHillSolver<TargetType, SourceType, CoordType>::~hoImageRegParametricDownHillSolver()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    typename hoImageRegParametricDownHillSolver<TargetType, SourceType, CoordType>::ValueType hoImageRegParametricDownHillSolver<TargetType, SourceType, CoordType>::solver_once(ValueType curr_dissimilarity)
    {
        ValueType prevValue = curr_dissimilarity;
        ValueType currValue;

        size_t optimizedDimIndex = 0;
        ValueType positiveStepFlag = 0;

        size_t numOfPara = transform_->get_number_of_parameters();
        size_t i;

        ValueType currPara(0);
        for ( i=0; i<numOfPara; i++ )
        {
            if ( transform_->get_para_status(i) == TransformationType::Active)
            {
                currPara = transform_->get_parameter(i);

                // positive
                transform_->set_parameter(i, currPara + step_size_para_[i]);

                GADGET_CHECK_RETURN_FALSE(warper_->warp(*target_, *source_, use_world_coordinate_, warpped_));
                currValue = dissimilarity_->evaluate(warpped_);

                if ( currValue < curr_dissimilarity + dissimilarity_thres_ )
                {
                    curr_dissimilarity = currValue;
                    optimizedDimIndex =  i;
                    positiveStepFlag = 1;
                }

                // negative
                transform_->set_parameter(i, currPara - step_size_para_[i]);

                GADGET_CHECK_RETURN_FALSE(warper_->warp(*target_, *source_, use_world_coordinate_, warpped_));
                currValue = dissimilarity_->evaluate(warpped_);

                if ( currValue < curr_dissimilarity + dissimilarity_thres_ )
                {
                    curr_dissimilarity = currValue;
                    optimizedDimIndex =  i;
                    positiveStepFlag = -1;
                }

                transform_->set_parameter(i, currPara);
            }
        }

        if ( curr_dissimilarity < prevValue )
        {
            currPara = transform_->get_parameter(optimizedDimIndex);
            transform_->set_parameter(optimizedDimIndex, currPara+positiveStepFlag*step_size_para_[optimizedDimIndex]);
        }

        return curr_dissimilarity;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegParametricDownHillSolver<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "-------------- Gagdgetron DownHill image registration solver -------------" << endl;
        BaseClass::print(os);
    }
}
#endif // hoImageRegParametricDownHillSolver_H_
