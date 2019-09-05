/** \file   hoImageRegWarper.h
    \brief  Define the class to perform image warpping using the geometric transformation in gadgetron registration
    \author Hui Xue
*/

#ifndef hoImageRegWarper_H_
#define hoImageRegWarper_H_

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
#include "hoImageRegDeformationField.h"

#include "GadgetronTimer.h"
#include "ImageIOAnalyze.h"

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

namespace Gadgetron {

    /// warp the source image to the grid of target image under a transformation
    /// both image domain warpping and world coordinate warpping is implemented
    /// for the image domain warpping, the pixels are in the coordinate of image grid
    /// input and output can have different dimensions
    /// input has DIn dimension and output has DOut dimension
    template<typename TargetType, typename SourceType, typename CoordType> 
    class hoImageRegWarper
    {
    public:

        typedef hoImageRegWarper<TargetType, SourceType, CoordType> Self;

        typedef typename TargetType::value_type ValueType;
        enum { DIn = TargetType::NDIM };
        enum { DOut = SourceType::NDIM };

        typedef hoNDImage<ValueType, 2> Target2DType;
        typedef Target2DType Source2DType;

        typedef hoNDImage<ValueType, 3> Target3DType;
        typedef Target2DType Source3DType;

        typedef hoNDInterpolator<SourceType> InterpolatorType;

        typedef hoImageRegTransformation<CoordType, DIn, DOut> TransformationType;
        typedef hoImageRegDeformationField<CoordType, DIn> DeformTransformationType;

        typedef ValueType T;
        typedef ValueType element_type;
        typedef ValueType value_type;

        typedef CoordType coord_type;

        typedef typename TransformationType::input_point_type input_point_type;
        typedef typename TransformationType::output_point_type output_point_type;

        typedef typename TransformationType::jacobian_parameter_type jacobian_parameter_type;
        typedef typename TransformationType::jacobian_position_type jacobian_position_type;

        hoImageRegWarper(ValueType bg_values = 0);
        virtual ~hoImageRegWarper();

        void setTransformation(TransformationType& transform);
        void setInterpolator(InterpolatorType& interp);
        void setBackgroundValue(ValueType bg_value);

        virtual bool warp(const TargetType& target, const SourceType& source, bool useWorldCoordinate, TargetType& warped);
        //virtual bool warp(const Target2DType& target, const Source2DType& source, bool useWorldCoordinate, Target2DType& warped);
        //virtual bool warp(const Target3DType& target, const Source3DType& source, bool useWorldCoordinate, Target3DType& warped);

        /// warp at the target image grid using the DeformationField transformation
        /// the DeformationField takes in the target pixel indexes and returns the transformed position in the world coordinates
        /// the deformation field grid should be the same as the target images
        virtual bool warpWithDeformationFieldWorldCoordinate(const TargetType& target, const SourceType& source, TargetType& warped);

        virtual void print(std::ostream& os) const;

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

        TransformationType* transform_;
        InterpolatorType* interp_;

        /// back ground values, used to mark regions in the target image which will not be warped
        ValueType bg_value_;
    };

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegWarper<TargetType, SourceType, CoordType>::hoImageRegWarper(ValueType bg_value) : transform_(NULL), interp_(NULL), performTiming_(false), bg_value_(bg_value)
    {
        gt_timer1_.set_timing_in_destruction(false);
        gt_timer2_.set_timing_in_destruction(false);
        gt_timer3_.set_timing_in_destruction(false);
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    hoImageRegWarper<TargetType, SourceType, CoordType>::~hoImageRegWarper()
    {
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    inline void hoImageRegWarper<TargetType, SourceType, CoordType>::setTransformation(TransformationType& transform)
    {
        transform_ = &transform;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    inline void hoImageRegWarper<TargetType, SourceType, CoordType>::setInterpolator(InterpolatorType& interp)
    {
        interp_ = &interp;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    inline void hoImageRegWarper<TargetType, SourceType, CoordType>::setBackgroundValue(ValueType bg_value)
    {
        bg_value_ = bg_value;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegWarper<TargetType, SourceType, CoordType>::
    warp(const TargetType& target, const SourceType& source, bool useWorldCoordinate, TargetType& warped)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(transform_!=NULL);

            if ( useWorldCoordinate )
            {
                // if the transformation is the deformation filed, special version of warp should be called
                DeformTransformationType* transformDeformField = dynamic_cast<DeformTransformationType*>(transform_);
                if( transformDeformField != NULL )
                {
                    return this->warpWithDeformationFieldWorldCoordinate(target, source, warped);
                }
            }

            GADGET_DEBUG_CHECK_RETURN_FALSE(interp_!=NULL);
            interp_->setArray( const_cast<SourceType&>(source) );

            warped = target;

            if ( DIn==2 && DOut==2 )
            {
                size_t sx = target.get_size(0);
                size_t sy = target.get_size(1);

                long long y;

                if ( useWorldCoordinate )
                {
                    // #pragma omp parallel private(y) shared(sx, sy, target, source, warped) num_threads(2)
                    {
                        typename TargetType::coord_type px, py, px_source, py_source, ix_source, iy_source;

                        // #pragma omp for 
                        for ( y=0; y<(long long)sy; y++ )
                        {
                            for ( size_t x=0; x<sx; x++ )
                            {
                                size_t offset = x + y*sx;

                                if ( target( offset ) != bg_value_ )
                                {
                                    // target to world
                                    target.image_to_world(x, size_t(y), px, py);

                                    // transform the point
                                    transform_->transform(px, py, px_source, py_source);

                                    // world to source
                                    source.world_to_image(px_source, py_source, ix_source, iy_source);

                                    // interpolate the source
                                    warped( offset ) = (*interp_)(ix_source, iy_source);
                                }
                            }
                        }
                    }
                }
                else
                {
                    // #pragma omp parallel private(y) shared(sx, sy, target, source, warped) num_threads(2)
                    {
                        typename TargetType::coord_type ix_source, iy_source;

                        // #pragma omp for 
                        for ( y=0; y<(long long)sy; y++ )
                        {
                            for ( size_t x=0; x<sx; x++ )
                            {
                                size_t offset = x + y*sx;

                                if ( target( offset ) != bg_value_ )
                                {
                                    // transform the point
                                    transform_->transform(x, size_t(y), ix_source, iy_source);

                                    // interpolate the source
                                    warped( offset ) = (*interp_)(ix_source, iy_source);
                                }
                            }
                        }
                    }
                }
            }
            else if ( DIn==3 && DOut==3 )
            {
                size_t sx = target.get_size(0);
                size_t sy = target.get_size(1);
                size_t sz = target.get_size(2);

                long long z;

                if ( useWorldCoordinate )
                {
                    #pragma omp parallel private(z) shared(sx, sy, sz, target, source, warped)
                    {
                        typename TargetType::coord_type px, py, pz, px_source, py_source, pz_source, ix_source, iy_source, iz_source;

                        #pragma omp for 
                        for ( z=0; z<(long long)sz; z++ )
                        {
                            for ( size_t y=0; y<sy; y++ )
                            {
                                size_t offset = y*sx + z*sx*sy;

                                for ( size_t x=0; x<sx; x++ )
                                {
                                    if ( target( x+offset ) != bg_value_ )
                                    {
                                        // target to world
                                        target.image_to_world(x, y, size_t(z), px, py, pz);

                                        // transform the point
                                        transform_->transform(px, py, pz, px_source, py_source, pz_source);

                                        // world to source
                                        source.world_to_image(px_source, py_source, pz_source, ix_source, iy_source, iz_source);

                                        // interpolate the source
                                        warped( x+offset ) = (*interp_)(ix_source, iy_source, iz_source);
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    #pragma omp parallel private(z) shared(sx, sy, sz, target, source, warped)
                    {
                        typename TargetType::coord_type ix_source, iy_source, iz_source;

                        #pragma omp for 
                        for ( z=0; z<(long long)sz; z++ )
                        {
                            for ( size_t y=0; y<sy; y++ )
                            {
                                size_t offset = y*sx + z*sx*sy;

                                for ( size_t x=0; x<sx; x++ )
                                {
                                    if ( target( x+offset ) != bg_value_ )
                                    {
                                        // transform the point
                                        transform_->transform(x, y, size_t(z), ix_source, iy_source, iz_source);

                                        // interpolate the source
                                        warped( x+offset ) = (*interp_)(ix_source, iy_source, iz_source);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                size_t numOfPixels = target.get_number_of_elements();

                long long n;

                if ( useWorldCoordinate )
                {
                    #pragma omp parallel private(n) shared(numOfPixels, target, source, warped)
                    {
                        size_t ind_target[DIn];
                        typename TargetType::coord_type pt_target[DIn];
                        typename TargetType::coord_type pt_source[DOut];
                        typename TargetType::coord_type ind_source[DOut];

                        #pragma omp for 
                        for ( n=0; n<(long long)numOfPixels; n++ )
                        {
                            if ( target( size_t(n) ) != bg_value_ )
                            {
                                // target to world
                                target.calculate_index( size_t(n), ind_target );

                                target.image_to_world(ind_target, pt_target);

                                // transform the point
                                transform_->transform(pt_target, pt_source);

                                // world to source
                                source.world_to_image(pt_source, ind_source);

                                // interpolate the source
                                warped( size_t(n) ) = (*interp_)(ind_source);
                            }
                        }
                    }
                }
                else
                {
                    #pragma omp parallel private(n) shared(numOfPixels, target, source, warped)
                    {
                        typename TargetType::coord_type pt_target[DIn];
                        typename TargetType::coord_type pt_source[DOut];

                        #pragma omp for 
                        for ( n=0; n<(long long)numOfPixels; n++ )
                        {
                            if ( target( size_t(n) ) != bg_value_ )
                            {
                                target.calculate_index( size_t(n), pt_target );

                                // transform the point
                                this->transform_->transform(pt_target, pt_source);

                                // interpolate the source
                                warped( size_t(n) ) = (*interp_)(pt_source);
                            }
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegWarper<TargetType, SourceType, CoordType>::\
                                    warp(const TargetType& target, const SourceType& source, bool useWorldCoordinate, TargetType& warped) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    bool hoImageRegWarper<TargetType, SourceType, CoordType>::
    warpWithDeformationFieldWorldCoordinate(const TargetType& target, const SourceType& source, TargetType& warped)
    {
        try
        {
            GADGET_DEBUG_CHECK_RETURN_FALSE(DIn==DOut);
            GADGET_DEBUG_CHECK_RETURN_FALSE(transform_!=NULL);

            DeformTransformationType* transformDeformField = dynamic_cast<DeformTransformationType*>(transform_);
            GADGET_DEBUG_CHECK_RETURN_FALSE(transformDeformField!=NULL);

            GADGET_DEBUG_CHECK_RETURN_FALSE(interp_!=NULL);
            interp_->setArray( const_cast<SourceType&>(source) );

            warped = target;

            if ( DIn==2 && DOut==2 )
            {
                size_t sx = target.get_size(0);
                size_t sy = target.get_size(1);

                long long y;

                // #pragma omp parallel private(y) shared(sx, sy, target, source, warped) num_threads(2)
                {
                    coord_type px, py, dx, dy, ix_source, iy_source;

                    // #pragma omp for 
                    for ( y=0; y<(long long)sy; y++ )
                    {
                        for ( size_t x=0; x<sx; x++ )
                        {
                            size_t offset = x + y*sx;

                            if ( target( offset ) != bg_value_ )
                            {
                                // target to world
                                target.image_to_world(x, size_t(y), px, py);

                                // transform the point
                                transformDeformField->get(x, size_t(y), dx, dy);

                                // world to source
                                source.world_to_image(px+dx, py+dy, ix_source, iy_source);

                                // interpolate the source
                                warped( offset ) = (*interp_)(ix_source, iy_source);
                            }
                        }
                    }
                }
            }
            else if ( DIn==3 && DOut==3 )
            {
                size_t sx = target.get_size(0);
                size_t sy = target.get_size(1);
                size_t sz = target.get_size(2);

                long long z;

                #pragma omp parallel private(z) shared(sx, sy, sz, target, source, warped)
                {
                    coord_type px, py, pz, dx, dy, dz, ix_source, iy_source, iz_source;

                    #pragma omp for 
                    for ( z=0; z<(long long)sz; z++ )
                    {
                        for ( size_t y=0; y<sy; y++ )
                        {
                            size_t offset = y*sx + z*sx*sy;

                            for ( size_t x=0; x<sx; x++ )
                            {
                                if ( target( x+offset ) != bg_value_ )
                                {
                                    // target to world
                                    target.image_to_world(x, y, size_t(z), px, py, pz);

                                    // transform the point
                                    transformDeformField->get(x, y, size_t(z), dx, dy, dz);

                                    // world to source
                                    source.world_to_image(px+dx, py+dy, pz+dz, ix_source, iy_source, iz_source);

                                    // interpolate the source
                                    warped( x+offset ) = (*interp_)(ix_source, iy_source, iz_source);
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                size_t numOfPixels = target.get_number_of_elements();

                long long n;

                #pragma omp parallel private(n) shared(numOfPixels, target, source, warped)
                {
                    size_t ind_target[DIn];
                    coord_type pt_target[DIn];
                    coord_type pt_source[DOut];
                    coord_type ind_source[DOut];

                    unsigned int ii;

                    #pragma omp for 
                    for ( n=0; n<(long long)numOfPixels; n++ )
                    {
                        if ( target( size_t(n) ) != bg_value_ )
                        {
                            // target to world
                            target.calculate_index( size_t(n), ind_target );

                            target.image_to_world(ind_target, pt_target);

                            // transform the point
                            transformDeformField->get(ind_target, pt_source);

                            for ( ii=0; ii<DIn; ii++ )
                            {
                                pt_source[ii] += pt_target[ii];
                            }

                            // world to source
                            source.world_to_image(pt_source, ind_source);

                            // interpolate the source
                            warped( size_t(n) ) = (*interp_)(ind_source);
                        }
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoImageRegWarper<TargetType, SourceType, CoordType>::\
                                    warpWithDeformationFieldWorldCoordinate(const TargetType& target, const SourceType& source, TargetType& warped) ... ");
            return false;
        }

        return true;
    }

    template<typename TargetType, typename SourceType, typename CoordType> 
    void hoImageRegWarper<TargetType, SourceType, CoordType>::print(std::ostream& os) const
    {
        using namespace std;
        os << "--------------Gagdgetron image warper -------------" << endl;
        os << "Input dimension is : " << DIn << endl;
        os << "Output dimension is : " << DOut << endl;

        std::string elemTypeName = std::string(typeid(ValueType).name());
        os << "Image data type is : " << elemTypeName << std::endl;

        elemTypeName = std::string(typeid(CoordType).name());
        os << "Transformation coordinate data type is : " << elemTypeName << std::endl;
    }
}
#endif // hoImageRegWarper_H_
