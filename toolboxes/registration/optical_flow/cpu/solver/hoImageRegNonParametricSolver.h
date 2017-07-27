/** \file   hoImageRegNonParametricSolver.h
    \brief  Define the base class of image registration solver for non-parametric image transformation
    \author Hui Xue
*/

#ifndef hoImageRegNonParametricSolver_H_
#define hoImageRegNonParametricSolver_H_

#pragma once

#include "hoImageRegSolver.h"

namespace Gadgetron {

    /// ValueType: image pixel value type
    /// CoordType: transformation data type
    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegNonParametricSolver : public hoImageRegSolver<TargetType, SourceType, CoordType>
    {
    public:

        typedef hoImageRegNonParametricSolver<TargetType, SourceType, CoordType> Self;
        typedef hoImageRegSolver<TargetType, SourceType, CoordType> BaseClass;

        typedef typename TargetType::value_type ValueType;
        enum { DIn = TargetType::NDIM };
        enum { DOut = SourceType::NDIM };

        typedef hoNDImage<ValueType, 2> Target2DType;
        typedef Target2DType Source2DType;

        typedef hoNDImage<ValueType, 3> Target3DType;
        typedef Target2DType Source3DType;

        typedef typename BaseClass::InterpolatorType InterpolatorType;

        typedef hoImageRegNonParametricTransformation<CoordType, DIn, DOut> TransformationType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef CoordType coord_type;

        typedef typename TransformationType::input_point_type input_point_type;
        typedef typename TransformationType::output_point_type output_point_type;
        typedef typename TransformationType::jacobian_position_type jacobian_position_type;

        typedef typename BaseClass::ImageRegWarperType ImageRegWarperType;

        typedef typename BaseClass::ImageRegDissimilarityType ImageRegDissimilarityType;

        hoImageRegNonParametricSolver();
        virtual ~hoImageRegNonParametricSolver();

        virtual bool initialize();

        /// solve the minimization and find the optimal transformation
        virtual bool solve() = 0;

        virtual void print(std::ostream& os) const;

        /// number of performed iterations
        unsigned int iter_num_;

        /// maximal number of iterations
        unsigned int max_iter_num_;

        /// threshold for minimal dissimilarity changes
        ValueType dissimilarity_thres_;

        /// threshold for minimal parameter changes
        ValueType parameter_thres_;

        /// number of search size division
        unsigned int div_num_;

        /// solver step size
        ValueType step_size_para_;
        /// step size division ratio
        ValueType step_size_div_para_;

        using BaseClass::verbose_;
        using BaseClass::gt_timer1_;
        using BaseClass::gt_timer2_;
        using BaseClass::gt_timer3_;
        using BaseClass::performTiming_;
        using BaseClass::gt_exporter_;
        using BaseClass::debugFolder_;

    protected:

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
    hoImageRegNonParametricSolver<TargetType, SourceType, CoordType>::hoImageRegNonParametricSolver() 
        : BaseClass(), dissimilarity_thres_(0), parameter_thres_( (ValueType)1e-8 ), div_num_(3), step_size_para_( (ValueType)0.8 ), step_size_div_para_( (ValueType)0.5 )
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegNonParametricSolver<TargetType, SourceType, CoordType>::~hoImageRegNonParametricSolver()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegNonParametricSolver<TargetType, SourceType, CoordType>::initialize()
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(target_!=NULL);
            GADGET_CHECK_RETURN_FALSE(source_!=NULL);
            GADGET_CHECK_RETURN_FALSE(interp_!=NULL);
            GADGET_CHECK_RETURN_FALSE(warper_!=NULL);
            GADGET_CHECK_RETURN_FALSE(dissimilarity_!=NULL);

            warper_->setInterpolator(*interp_);
            warper_->setBackgroundValue(bg_value_);

            dissimilarity_->setBackgroundValue(bg_value_);

            if ( !warpped_.dimensions_equal(*target_) )
            {
                warpped_ = *target_;
            }

            dissimilarity_->initialize(*target_);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegNonParametricSolver<TargetType, SourceType, CoordType>::initialize() ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegNonParametricSolver<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image registration non-parametric solver -------------" << endl;
        os << "Target image dimension is : " << DIn << endl;
        os << "Source image dimension is : " << DOut << endl;
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
#endif // hoImageRegNonParametricSolver_H_
