/** \file   hoImageRegSolver.h
    \brief  Define the base class of image registration solver for gadgetron

            The solver takes in the image warper, similarity, target and source images, and solves
            for an optimal image transformation.

    \author Hui Xue
*/

#ifndef hoImageRegSolver_H_
#define hoImageRegSolver_H_

#pragma once

#include "hoNDArray.h"
#include "hoNDImage.h"
#include "hoNDInterpolator.h"
#include "hoNDBoundaryHandler.h"
#include "hoMatrix.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"
#include "hoImageRegTransformation.h"
#include "hoImageRegWarper.h"
#include "hoImageRegDissimilarity.h"

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

namespace Gadgetron {

    // define the solver type
    enum GT_IMAGE_REG_SOLVER
    {
        GT_IMAGE_REG_SOLVER_DOWNHILL,
        GT_IMAGE_REG_SOLVER_GRADIENT_DESCENT,
        GT_IMAGE_REG_SOLVER_PDE_TIME_INTEGRATION,
        GT_IMAGE_REG_SOLVER_PDE_TIME_INTEGRATION_INV
    };

    inline std::string getImageRegSolverName(GT_IMAGE_REG_SOLVER v)
    {
        std::string name;

        switch (v)
        {
            case GT_IMAGE_REG_SOLVER_DOWNHILL:
                name = "DownHill";
                break;

            case GT_IMAGE_REG_SOLVER_GRADIENT_DESCENT:
                name = "GradientDescent";
                break;

            case GT_IMAGE_REG_SOLVER_PDE_TIME_INTEGRATION:
                name = "PDE_Time_Integration";
                break;

            case GT_IMAGE_REG_SOLVER_PDE_TIME_INTEGRATION_INV:
                name = "PDE_Time_Integration_Inv";
                break;

            default:
                GERROR_STREAM("Unrecognized image registration solver type : " << v);
        }

        return name;
    }

    inline GT_IMAGE_REG_SOLVER getImageRegSolverType(const std::string& name)
    {
        GT_IMAGE_REG_SOLVER v;

        if ( name == "DownHill" )
        {
            v = GT_IMAGE_REG_SOLVER_DOWNHILL;
        }
        else if ( name == "GradientDescent" )
        {
            v = GT_IMAGE_REG_SOLVER_GRADIENT_DESCENT;
        }
        else if ( name == "PDE_Time_Integration" )
        {
            v = GT_IMAGE_REG_SOLVER_PDE_TIME_INTEGRATION;
        }
        else if ( name == "PDE_Time_Integration_Inv" )
        {
            v = GT_IMAGE_REG_SOLVER_PDE_TIME_INTEGRATION_INV;
        }
        else
        {
            GERROR_STREAM("Unrecognized image registration solver name : " << name);
        }

        return v;
    }

    /// ValueType: image pixel value type
    /// CoordType: transformation data type
    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegSolver
    {
    public:

        typedef hoImageRegSolver<TargetType, SourceType, CoordType> Self;

        typedef typename TargetType::value_type ValueType;
        enum { DIn = TargetType::NDIM };
        enum { DOut = SourceType::NDIM };

        typedef hoNDImage<ValueType, 2> Target2DType;
        typedef Target2DType Source2DType;

        typedef hoNDImage<ValueType, 3> Target3DType;
        typedef Target2DType Source3DType;

        typedef hoNDInterpolator<SourceType> InterpolatorType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef CoordType coord_type;

        typedef hoImageRegWarper<TargetType, SourceType, CoordType> ImageRegWarperType;

        typedef hoImageRegDissimilarity<SourceType> ImageRegDissimilarityType;

        hoImageRegSolver();
        virtual ~hoImageRegSolver();

        void setTarget(TargetType& target) { target_ = &target; }
        void setSource(SourceType& source) { source_ = &source; }

        void setDissimilarity(ImageRegDissimilarityType& dissimilarity) { dissimilarity_ = &dissimilarity; }
        void setWarper(ImageRegWarperType& warper) { warper_ = &warper; }
        void setInterpolator(InterpolatorType& interp) { interp_ = &interp; }
        void setBackgroundValue(ValueType bg_value) { bg_value_ = bg_value; }

        void setUseWorldCoordinate(bool use_world_coordinate) { use_world_coordinate_ = use_world_coordinate; }

        virtual bool solve() = 0;

        virtual void print(std::ostream& os) const;

        /// if true, print out more intermediate information
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

    protected:

        TargetType* target_;
        SourceType* source_;

        // warped image
        TargetType warpped_;

        ValueType bg_value_;

        InterpolatorType* interp_;

        ImageRegWarperType* warper_;

        ImageRegDissimilarityType* dissimilarity_;

        /// whether to perform registration using the world coordinates
        bool use_world_coordinate_;
    };

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegSolver<TargetType, SourceType, CoordType>::hoImageRegSolver() 
        : target_(NULL), source_(NULL), bg_value_(0), interp_(NULL), warper_(NULL), dissimilarity_(NULL), verbose_(false), use_world_coordinate_(true), performTiming_(false)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegSolver<TargetType, SourceType, CoordType>::~hoImageRegSolver()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegSolver<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image registration solver -------------" << endl;
        os << "Target image dimension is : " << DIn << endl;
        os << "Source image dimension is : " << DOut << endl;
        os << "Image data type is : " << std::string(typeid(ValueType).name()) << std::endl;
        os << "Transformation data type is : " << std::string(typeid(CoordType).name()) << std::endl;
        os << "Use world coordinate is : " << use_world_coordinate_ << std::endl;
        os << "verbose flag is : " << verbose_ << std::endl;
    }
}
#endif // hoImageRegSolver_H_
