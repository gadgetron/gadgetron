/** \file   hoImageRegNonParametricRegister.h
    \brief  Define the class to perform non-parametric image registration in gadgetron
    \author Hui Xue
*/

#ifndef hoImageRegNonParametricRegister_H_
#define hoImageRegNonParametricRegister_H_

#pragma once

#include "hoImageRegRegister.h"

namespace Gadgetron {

    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegNonParametricRegister : public hoImageRegRegister<TargetType, SourceType, CoordType>
    {
    public:

        typedef hoImageRegNonParametricRegister<TargetType, SourceType, CoordType> Self;
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

        hoImageRegNonParametricRegister(unsigned int resolution_pyramid_levels=3, ValueType bg_value=ValueType(0));
        virtual ~hoImageRegNonParametricRegister();

        /// initialize the registration
        /// should be called after all images and parameters of registration are set
        virtual bool initialize() { return BaseClass::initialize(); }

        /// perform the registration
        virtual bool performRegistration() = 0;

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
    };

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegNonParametricRegister<TargetType, SourceType, CoordType>::
    hoImageRegNonParametricRegister(unsigned int resolution_pyramid_levels, ValueType bg_value) : BaseClass(resolution_pyramid_levels, bg_value)
    {

    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegNonParametricRegister<TargetType, SourceType, CoordType>::~hoImageRegNonParametricRegister()
    {

    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegNonParametricRegister<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron non-parametric image register -------------" << endl;
        BaseClass::printContent(os);
        os << "--------------------------------------------------------------------" << endl << ends;
    }
}
#endif // hoImageRegNonParametricRegister_H_
