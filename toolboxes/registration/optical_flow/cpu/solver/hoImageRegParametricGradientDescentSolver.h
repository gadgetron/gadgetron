/** \file   hoImageRegParametricGradientDescentSolver.h
    \brief  Define the class of simple gradient descent solver for parametric image transformation, no linear search is performed in this solver
    \author Hui Xue
*/

#ifndef hoImageRegParametricGradientDescentSolver_H_
#define hoImageRegParametricGradientDescentSolver_H_

#pragma once

#include "hoImageRegParametricSolver.h"

namespace Gadgetron {

    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegParametricGradientDescentSolver : public hoImageRegParametricSolver<TargetType, SourceType, CoordType>
    {
    public:

        typedef hoImageRegParametricGradientDescentSolver<TargetType, SourceType, CoordType> Self;
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

        hoImageRegParametricGradientDescentSolver();
        virtual ~hoImageRegParametricGradientDescentSolver();

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
    hoImageRegParametricGradientDescentSolver<TargetType, SourceType, CoordType>::hoImageRegParametricGradientDescentSolver() : BaseClass()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegParametricGradientDescentSolver<TargetType, SourceType, CoordType>::~hoImageRegParametricGradientDescentSolver()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    typename hoImageRegParametricGradientDescentSolver<TargetType, SourceType, CoordType>::ValueType hoImageRegParametricGradientDescentSolver<TargetType, SourceType, CoordType>::solver_once(ValueType curr_dissimilarity)
    {
        std::vector<ValueType> deriv;
        GADGET_CHECK_RETURN_FALSE(this->evaluateDeriv(transform_, dissimilarity_, step_size_para_, deriv));

        size_t numOfPara = transform_->get_number_of_parameters();
        size_t i;
        ValueType prevValue, currPara(0);

        while ( 1 )
        {
            prevValue = curr_dissimilarity;

            for ( i=0; i<numOfPara; i++ )
            {
                currPara = transform_->get_parameter(i);
                transform_->set_parameter( i, currPara-step_size_para_[i]*deriv[i] );
            }

            GADGET_CHECK_RETURN_FALSE(warper_->warp(*target_, *source_, use_world_coordinate_, warpped_));
            curr_dissimilarity = dissimilarity_->evaluate(warpped_);

            if ( curr_dissimilarity > prevValue + dissimilarity_thres_ )
            {
                break;
            }
        }

        // rewind
        for ( i=0; i<numOfPara; i++ )
        {
            currPara = transform_->get_parameter(i);
            transform_->set_parameter( i, currPara+step_size_para_[i]*deriv[i] );
        }

        return prevValue;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegParametricGradientDescentSolver<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "-------------- Gagdgetron Gradient descent image registration solver -------------" << endl;
        BaseClass::print(os);
    }
}
#endif // hoImageRegParametricGradientDescentSolver_H_
